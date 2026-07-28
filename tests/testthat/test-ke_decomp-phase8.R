.phase8_test_data <- function() {
  n_strata <- 6L
  psu_per_stratum <- 4L
  observations_per_psu <- 5L
  n <- n_strata * psu_per_stratum * observations_per_psu
  source_id <- seq_len(n)
  data.frame(
    health = 0.27 +
      0.045 * sin(source_id / 8) +
      0.0025 * source_id +
      0.02 * as.numeric(source_id %% 4 == 0),
    ses = rep(seq_len(30), each = 4),
    income = 190 + 6 * source_id +
      9 * cos(source_id / 7),
    x1 = sin(source_id / 6) + source_id / 100,
    x2 = as.numeric(source_id %% 4 == 0),
    weight = 0.8 +
      0.1 * rep(
        seq_len(n_strata),
        each = psu_per_stratum * observations_per_psu
      ) +
      0.08 * rep(
        rep(
          seq_len(psu_per_stratum),
          each = observations_per_psu
        ),
        n_strata
      ),
    strata = rep(
      seq_len(n_strata),
      each = psu_per_stratum * observations_per_psu
    ),
    psu = rep(
      rep(
        seq_len(psu_per_stratum),
        each = observations_per_psu
      ),
      n_strata
    )
  )
}

.phase8_manual_design_vcov <- function(fit) {
  matrix <- fit$raw$model$matrix
  residual <- fit$raw$model$residuals
  weights <- unname(fit$raw$design$weights)
  strata <- unname(fit$raw$design$strata)
  psu <- unname(fit$raw$design$psu)
  bread <- solve(crossprod(matrix, weights * matrix))
  score <- matrix * (weights * residual)
  nested_psu <- interaction(
    strata,
    psu,
    drop = TRUE,
    lex.order = TRUE
  )
  psu_score <- rowsum(score, nested_psu, reorder = FALSE)
  psu_strata <- vapply(
    split(strata, nested_psu),
    `[`,
    character(1),
    1L
  )
  meat <- matrix(
    0,
    nrow = ncol(matrix),
    ncol = ncol(matrix)
  )
  for (stratum in unique(psu_strata)) {
    stratum_score <- psu_score[
      psu_strata == stratum,
      ,
      drop = FALSE
    ]
    psu_count <- nrow(stratum_score)
    centered <- sweep(
      stratum_score,
      2,
      colMeans(stratum_score),
      "-"
    )
    meat <- meat +
      psu_count / (psu_count - 1) *
      crossprod(centered)
  }
  covariance <- bread %*% meat %*% bread
  dimnames(covariance) <- dimnames(fit$raw$model$vcov)
  covariance
}

test_that("KE phase 8 weight-only survey matches conditional linearization", {
  data <- .phase8_test_data()
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    use_svy = FALSE,
    weight_var = "weight",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expect_true(fit$raw$settings$use_svy)
  expect_equal(fit$summary_stats$N, nrow(data))
  expect_equal(fit$summary_stats$Strata, 1)
  expect_equal(fit$summary_stats$PSUs, nrow(data))
  expect_equal(fit$summary_stats$DF, nrow(data) - 1)
  expect_equal(
    fit$raw$model$vcov,
    .phase8_manual_design_vcov(fit),
    tolerance = 2e-14
  )
  expect_match(
    fit$raw$model$vcov_type,
    "survey design covariance conditional"
  )
  expect_identical(
    fit$raw$vcov$scope,
    "conditional_on_empirically_constructed_target"
  )
  expect_false(
    fit$raw$vcov$includes_generated_target_uncertainty
  )
  expect_true(is.na(fit$results_overall$Std_Error))
  expect_null(fit$raw$replication)
})

test_that("KE phase 8 complex survey uses nested PSU design degrees", {
  data <- .phase8_test_data()
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "income",
    index_type = "level",
    correction = "standard",
    weight_var = "weight",
    strata_var = "strata",
    psu_var = "psu",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expected_psu <- nrow(unique(data[c("strata", "psu")]))
  expected_strata <- length(unique(data$strata))
  expected_df <- expected_psu - expected_strata
  expect_true(fit$raw$settings$use_svy)
  expect_equal(fit$summary_stats$Strata, expected_strata)
  expect_equal(fit$summary_stats$PSUs, expected_psu)
  expect_equal(fit$summary_stats$DF, expected_df)
  expect_equal(
    fit$raw$design$degrees_of_freedom,
    expected_df
  )
  expect_equal(
    fit$raw$model$vcov,
    .phase8_manual_design_vcov(fit),
    tolerance = 2e-14
  )

  standard_errors <- sqrt(diag(fit$raw$model$vcov))
  statistics <- fit$raw$model$coefficients / standard_errors
  p_values <- 2 * stats::pt(
    abs(statistics),
    df = expected_df,
    lower.tail = FALSE
  )
  critical_value <- stats::qt(0.975, df = expected_df)
  expect_equal(
    fit$results_detailed$Std_Error,
    unname(standard_errors),
    tolerance = 2e-14
  )
  expect_equal(
    fit$results_detailed$Statistic,
    unname(statistics),
    tolerance = 2e-13
  )
  expect_equal(
    fit$results_detailed$P_Value,
    unname(p_values),
    tolerance = 2e-14
  )
  expect_equal(
    fit$results_detailed$Conf_Low,
    unname(
      fit$raw$model$coefficients -
        critical_value * standard_errors
    ),
    tolerance = 2e-14
  )
  expect_equal(
    fit$results_detailed$Conf_High,
    unname(
      fit$raw$model$coefficients +
        critical_value * standard_errors
    ),
    tolerance = 2e-14
  )

  x1_row <- fit$results_anova[
    fit$results_anova$Variable == "x1",
    ,
    drop = FALSE
  ]
  expected_f <- statistics[["x1"]]^2
  expect_equal(
    x1_row$F_Statistic,
    expected_f,
    tolerance = 2e-13
  )
  expect_equal(
    x1_row$Prob_F,
    stats::pf(
      expected_f,
      1,
      expected_df,
      lower.tail = FALSE
    ),
    tolerance = 2e-14
  )
})

test_that("KE phase 8 strata and PSU activate unit-weight survey mode", {
  data <- .phase8_test_data()
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    use_svy = FALSE,
    weight_var = NULL,
    strata_var = "strata",
    psu_var = "psu",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expect_true(fit$raw$settings$use_svy)
  expect_equal(
    unname(fit$raw$design$weights),
    rep(1, nrow(data)),
    tolerance = 0
  )
  expect_equal(fit$summary_stats$Pop, nrow(data))
  expect_equal(fit$summary_stats$Strata, 6)
  expect_equal(fit$summary_stats$PSUs, 24)
  expect_equal(fit$summary_stats$DF, 18)
  expect_equal(
    fit$raw$model$vcov,
    .phase8_manual_design_vcov(fit),
    tolerance = 2e-14
  )
})

test_that("KE phase 8 is invariant to global survey-weight scaling", {
  data <- .phase8_test_data()
  data$weight_scaled <- 100 * data$weight
  arguments <- list(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "erreygers",
    dep_min = 0,
    dep_max = 1,
    strata_var = "strata",
    psu_var = "psu",
    vce_method = "linearized",
    quiet = TRUE
  )
  original <- suppressWarnings(do.call(
    ke_decomp,
    c(arguments, list(weight_var = "weight"))
  ))
  scaled <- suppressWarnings(do.call(
    ke_decomp,
    c(arguments, list(weight_var = "weight_scaled"))
  ))

  expect_equal(
    original$results_overall$Estimate,
    scaled$results_overall$Estimate,
    tolerance = 2e-14
  )
  expect_equal(
    original$raw$target$position,
    scaled$raw$target$position,
    tolerance = 2e-14
  )
  expect_equal(
    original$raw$model$coefficients,
    scaled$raw$model$coefficients,
    tolerance = 2e-14
  )
  expect_equal(
    original$raw$model$vcov,
    scaled$raw$model$vcov,
    tolerance = 2e-14
  )
  expect_equal(
    original$results_detailed$P_Value,
    scaled$results_detailed$P_Value,
    tolerance = 2e-14
  )
  expect_equal(
    scaled$summary_stats$Pop,
    100 * original$summary_stats$Pop,
    tolerance = 2e-12
  )
})

test_that("KE phase 8 validates the survey activation flag", {
  data <- .phase8_test_data()
  expect_error(
    ke_decomp(
      data = data,
      dep_var = "health",
      indep_vars = c("x1", "x2"),
      ses_var = "ses",
      correction = "generalized",
      use_svy = NA,
      quiet = TRUE
    ),
    "use_svy must be exactly TRUE or FALSE"
  )
  expect_error(
    ke_decomp(
      data = data,
      dep_var = "health",
      indep_vars = c("x1", "x2"),
      ses_var = "ses",
      correction = "generalized",
      use_svy = c(TRUE, FALSE),
      quiet = TRUE
    ),
    "use_svy must be exactly TRUE or FALSE"
  )
})
