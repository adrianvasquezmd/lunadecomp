.phase6_scope_data <- function(n = 36L) {
  source_id <- seq_len(n)
  ses <- 300 * exp(
    seq(-0.5, 0.8, length.out = n) +
      0.04 * sin(source_id / 4)
  )
  x1 <- sin(source_id / 5) + source_id / 50
  x2 <- as.integer(source_id %% 3 == 0)
  health <- stats::plogis(
    -0.9 +
      0.4 * x1 -
      0.2 * x2 +
      0.3 * as.numeric(scale(log(ses)))
  )
  data.frame(
    health = health,
    ses = ses,
    x1 = x1,
    x2 = x2
  )
}

.phase6_scope_fit <- function(index_type = "rank") {
  suppressWarnings(ke_decomp(
    data = .phase6_scope_data(),
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = index_type,
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))
}

test_that("KE phase 6 labels conditional linearization precisely", {
  for (index_type in c("rank", "level")) {
    fit <- .phase6_scope_fit(index_type)

    expect_identical(
      fit$raw$vcov$scope,
      "conditional_on_empirically_constructed_target"
    )
    expect_false(
      fit$raw$vcov$includes_generated_target_uncertainty
    )
    expect_true(
      fit$raw$vcov$model_complete_for_constructed_target
    )
    expect_false(
      fit$raw$vcov$model_complete_for_full_estimator
    )
    expect_false(fit$raw$vcov$overall_complete)
    expect_null(
      fit$raw$vcov$joint_index_coefficient_covariance
    )
    expect_identical(
      fit$raw$vcov$recommended_complete_methods,
      c("jackknife", "bootstrap")
    )
    expect_true(is.na(fit$results_overall$Std_Error))
    expect_true(all(is.na(fit$raw$vcov$overall)))
    expect_match(
      fit$diagnostics$vce,
      "package extension"
    )
    expect_match(
      fit$raw$known_limitations$linearized,
      "not specified"
    )
    expect_match(
      fit$raw$known_limitations$overall_inference,
      "point estimate without standard error"
    )
    expect_match(
      fit$diagnostics$overall_inference,
      "point estimate only"
    )
  }
})

test_that("KE phase 6 preserves Stata HC1 under conditional scope", {
  fit <- .phase6_scope_fit("rank")
  matrix <- fit$raw$model$matrix
  residual <- fit$raw$model$residuals
  model_rank <- qr(matrix)$rank
  bread <- solve(crossprod(matrix))
  expected <- nrow(matrix) / (
    nrow(matrix) - model_rank
  ) *
    bread %*%
    crossprod(matrix * residual) %*%
    bread

  expect_equal(
    fit$raw$model$vcov,
    expected,
    tolerance = 2e-14
  )
  expect_match(
    fit$raw$model$vcov_type,
    "conditional on the constructed target"
  )
})

test_that("KE phase 6 limits complete recomputation to coefficients", {
  fit <- suppressWarnings(ke_decomp(
    data = .phase6_scope_data(24L),
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "jackknife",
    quiet = TRUE
  ))

  expect_identical(
    fit$raw$vcov$scope,
    "coefficient_covariance_with_complete_target_recomputation"
  )
  expect_true(
    fit$raw$vcov$includes_generated_target_uncertainty
  )
  expect_true(
    fit$raw$vcov$model_complete_for_full_estimator
  )
  expect_true(
    fit$raw$vcov[[
      "coefficient_covariance_includes_target_reconstruction"
    ]]
  )
  expect_false(fit$raw$vcov$overall_complete)
  expect_null(
    fit$raw$vcov$joint_index_coefficient_covariance
  )
  expect_true(is.na(fit$results_overall$Std_Error))
  expect_null(fit$raw$replication$index_variance)
  expect_identical(
    colnames(fit$raw$replication$vcov),
    c("(Intercept)", "x1", "x2")
  )
  expect_identical(
    fit$raw$replication$vcov,
    fit$raw$model$vcov
  )
  expect_identical(
    fit$raw$vcov$recommended_complete_methods,
    character()
  )
})

test_that("KE phase 6 keeps the index point-only for every ordinary VCE", {
  for (method in c("linearized", "jackknife", "bootstrap")) {
    fit <- suppressWarnings(ke_decomp(
      data = .phase6_scope_data(24L),
      dep_var = "health",
      indep_vars = c("x1", "x2"),
      ses_var = "ses",
      index_type = "rank",
      correction = "generalized",
      vce_method = method,
      boot_reps = 10L,
      seed = 62026L,
      quiet = TRUE
    ))

    expect_true(is.finite(fit$results_overall$Estimate))
    expect_true(all(is.na(
      fit$results_overall[
        c(
          "Std_Error", "Statistic", "P_Value",
          "Conf_Low", "Conf_High"
        )
      ]
    )))
    expect_false(fit$raw$vcov$overall_complete)
  }
})

test_that("KE phase 6 keeps the complex-survey index point-only", {
  data <- .phase6_scope_data(40L)
  data$weight <- 1 + (seq_len(nrow(data)) %% 5L) / 5
  data$strata <- rep(seq_len(4L), each = 10L)
  data$psu <- rep(rep(seq_len(5L), each = 2L), 4L)

  for (method in c("linearized", "jackknife", "bootstrap")) {
    fit <- suppressWarnings(ke_decomp(
      data = data,
      dep_var = "health",
      indep_vars = c("x1", "x2"),
      ses_var = "ses",
      index_type = "rank",
      correction = "generalized",
      use_svy = TRUE,
      weight_var = "weight",
      strata_var = "strata",
      psu_var = "psu",
      vce_method = method,
      boot_reps = 10L,
      seed = 62026L,
      quiet = TRUE
    ))

    expect_true(is.finite(fit$results_overall$Estimate))
    expect_true(all(is.na(
      fit$results_overall[
        c(
          "Std_Error", "Statistic", "P_Value",
          "Conf_Low", "Conf_High"
        )
      ]
    )))
    expect_identical(
      colnames(fit$raw$model$vcov),
      colnames(fit$raw$model$matrix)
    )
  }
})
