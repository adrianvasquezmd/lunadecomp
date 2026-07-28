.phase5_model_data <- function() {
  n <- 144L
  source_id <- seq_len(n)
  education <- rep(c(1L, 2L, 3L), length.out = n)
  region <- rep(
    rep(c(1L, 2L, 3L), each = 4L),
    length.out = n
  )
  ses <- rep(
    c(1, 2, 4, 7, 11, 16, 23, 31, 42, 55, 70, 88),
    each = 12L
  )
  income <- 420 + 7.5 * source_id +
    31 * sin(source_id / 8)
  x1 <- sin(source_id / 6) + source_id / 55 -
    0.16 * (education == 3L)
  dummy <- as.integer(
    (5 * source_id + education) %% 9 < 4
  )
  health <- 0.31 +
    0.022 * x1 -
    0.017 * dummy +
    0.026 * (education == 2L) +
    0.051 * (education == 3L) -
    0.019 * (region == 2L) +
    0.013 * (region == 3L) +
    0.018 * x1 * (education == 2L) -
    0.012 * x1 * (education == 3L) +
    0.0011 * ses +
    0.014 * cos(source_id / 5)
  input_order <- c(
    seq(3, n, by = 3),
    seq(1, n, by = 3),
    seq(2, n, by = 3)
  )
  data.frame(
    source_id = source_id,
    health = health,
    ses = ses,
    income = income,
    x1 = x1,
    x1_duplicate = x1,
    dummy = dummy,
    education = factor(education, levels = 1:3),
    education_ref3 = factor(
      education,
      levels = c(3, 1, 2)
    ),
    region = factor(region, levels = 1:3),
    stringsAsFactors = FALSE
  )[input_order, , drop = FALSE]
}

.phase5_model_fit <- function(
    terms,
    index_type = "rank") {
  data <- .phase5_model_data()
  suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = terms,
    ses_var = if (index_type == "rank") "ses" else "income",
    index_type = index_type,
    correction = "generalized",
    use_svy = FALSE,
    vce_method = "linearized",
    quiet = TRUE
  ))
}

test_that("KE phase 5 accepts factor interactions as model terms", {
  fit <- .phase5_model_fit("x1 * education")

  expect_identical(
    fit$raw$model$term_labels,
    c("x1", "education", "x1:education")
  )
  expect_identical(
    names(fit$raw$model$coefficients),
    c(
      "(Intercept)", "x1", "education2", "education3",
      "x1:education2", "x1:education3"
    )
  )
  expect_identical(
    names(fit$raw$model$term_groups),
    c("x1", "education", "x1:education")
  )
  expect_identical(
    fit$raw$model$factor_levels$education,
    c("1", "2", "3")
  )
  expect_true("Omitted" %in% names(fit$results_detailed))

  separate_terms <- .phase5_model_fit(
    c("x1", "education", "x1:education")
  )
  expect_equal(
    separate_terms$raw$model$matrix,
    fit$raw$model$matrix,
    tolerance = 0
  )
  expect_equal(
    separate_terms$raw$model$coefficients,
    fit$raw$model$coefficients,
    tolerance = 0
  )
})

test_that("KE phase 5 interaction coefficients match frozen Stata HC1", {
  fit <- .phase5_model_fit("x1 * education")
  expected <- c(
    "(Intercept)" = -0.25967960401743323,
    x1 = 0.19333779132201431,
    education2 = 0.0018627041021214389,
    education3 = 0.058655218457211472,
    "x1:education2" = 0.036031579255758471,
    "x1:education3" = 0.000066777092811667492
  )
  expected_se <- c(
    "(Intercept)" = 0.031705184988673231,
    x1 = 0.017102445268154135,
    education2 = 0.04774067485302181,
    education3 = 0.046031626530587053,
    "x1:education2" = 0.025907116425704913,
    "x1:education3" = 0.025010376908325863
  )

  expect_equal(
    fit$raw$model$coefficients,
    expected,
    tolerance = 5e-14
  )
  expect_equal(
    fit$results_detailed$Std_Error,
    unname(expected_se),
    tolerance = 5e-14
  )
  expect_equal(fit$summary_stats$DF, 138)
  expect_match(
    fit$raw$model$vcov_type,
    "Stata-style HC1"
  )
})

test_that("KE phase 5 ordinary linearized covariance is HC1", {
  fit <- .phase5_model_fit(
    "x1 + dummy + education + region + x1:education"
  )
  matrix <- fit$raw$model$matrix
  residual <- fit$raw$model$residuals
  model_rank <- qr(matrix)$rank
  bread <- solve(crossprod(matrix))
  hc1 <- nrow(matrix) / (nrow(matrix) - model_rank) *
    bread %*%
    crossprod(matrix * residual) %*%
    bread

  expect_equal(
    fit$raw$model$vcov,
    hc1,
    tolerance = 5e-14
  )
  expect_equal(
    fit$summary_stats$DF,
    nrow(matrix) - model_rank
  )
  expect_equal(
    fit$results_detailed$P_Value,
    2 * stats::pt(
      abs(fit$results_detailed$Statistic),
      df = fit$summary_stats$DF,
      lower.tail = FALSE
    ),
    tolerance = 0
  )
})

test_that("KE phase 5 block F tests and logworth match Stata", {
  fit <- .phase5_model_fit("x1 * education")
  anova <- as.data.frame(fit$results_anova)
  anova <- anova[match(
    c("x1", "education", "x1:education"),
    anova$Variable
  ), ]

  expect_equal(
    anova$DF,
    c(1, 2, 2)
  )
  expect_equal(
    anova$F_Statistic,
    c(
      127.79594900201764,
      0.99885373847776215,
      1.2125836596202226
    ),
    tolerance = 6e-12
  )
  expect_equal(
    anova$Prob_F,
    c(
      2.2118779779116892e-21,
      0.37094814813338095,
      0.30057676805924055
    ),
    tolerance = 2e-13
  )
  expect_equal(
    anova$Logworth,
    -log10(anova$Prob_F),
    tolerance = 2e-14
  )
  expect_gt(
    fit$results_detailed$P_Value[
      fit$raw$model$term_groups$x1
    ],
    0
  )
})

test_that("KE phase 5 factor block is invariant to reference level", {
  base <- .phase5_model_fit("x1 + education")
  ref3 <- .phase5_model_fit("x1 + education_ref3")

  expect_equal(
    base$results_overall$Estimate,
    ref3$results_overall$Estimate,
    tolerance = 1e-15
  )
  expect_equal(
    base$model_metrics$R_Squared,
    ref3$model_metrics$R_Squared,
    tolerance = 1e-15
  )
  expect_equal(
    base$raw$model$fitted,
    ref3$raw$model$fitted,
    tolerance = 2e-15
  )
  base_factor <- base$results_anova[
    base$results_anova$Variable == "education",
  ]
  ref3_factor <- ref3$results_anova[
    ref3$results_anova$Variable == "education_ref3",
  ]
  expect_equal(
    base_factor$F_Statistic,
    ref3_factor$F_Statistic,
    tolerance = 2e-14
  )
  expect_equal(
    ref3$raw$model$coefficients[c(
      "education_ref31", "education_ref32"
    )],
    c(
      education_ref31 = -0.060380009126072651,
      education_ref32 = -0.010407172458313749
    ),
    tolerance = 5e-14
  )
})

test_that("KE phase 5 marks aliased predictors as omitted", {
  fit <- .phase5_model_fit(
    c("x1", "x1_duplicate", "dummy")
  )
  duplicate_row <- match(
    "x1_duplicate",
    names(fit$raw$model$coefficients)
  )

  expect_true(fit$raw$model$aliased["x1_duplicate"])
  expect_equal(
    unname(fit$raw$model$coefficients["x1_duplicate"]),
    0,
    tolerance = 0
  )
  expect_equal(
    unname(fit$raw$model$vcov["x1_duplicate", ]),
    rep(0, ncol(fit$raw$model$vcov)),
    tolerance = 0
  )
  expect_true(fit$results_detailed$Omitted[duplicate_row])
  expect_true(is.na(
    fit$results_detailed$P_Value[duplicate_row]
  ))
  expect_match(
    fit$diagnostics$collinearity,
    "reported as omitted"
  )
  expect_equal(fit$summary_stats$DF, 141)
})

test_that("KE phase 5 full factor model also matches level target", {
  fit <- .phase5_model_fit(
    "x1 + dummy + education + region + x1:education",
    index_type = "level"
  )

  expect_equal(
    fit$results_overall$Estimate,
    0.014028379115368242,
    tolerance = 1e-12
  )
  expect_equal(
    fit$raw$model$coefficients[c("x1", "education3")],
    c(
      x1 = 0.12460768413259406,
      education3 = 0.057392218882205914
    ),
    tolerance = 5e-14
  )
  expect_equal(
    fit$results_anova$F_Statistic[
      fit$results_anova$Variable == "x1:education"
    ],
    2.1369862065516121,
    tolerance = 6e-12
  )
})
