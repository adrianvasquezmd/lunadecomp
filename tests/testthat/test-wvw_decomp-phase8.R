phase8_ols_data <- function() {
  n <- 72
  id <- seq_len(n)
  x1 <- sin(id / 7) + id / 40
  x2 <- as.numeric(id %% 3 == 0)
  data.frame(
    y = 3.4 + 0.65 * x1 - 0.42 * x2 +
      0.18 * cos(id / 5),
    x1 = x1,
    x2 = x2,
    ses = id + (id %% 5) / 10,
    weight = 1 + (id %% 7) / 3
  )
}

phase8_binary_data <- function() {
  n <- 180
  id <- seq_len(n)
  x1 <- sin(id / 11) + (id %% 13) / 10
  x2 <- as.numeric(id %% 4 < 2)
  region <- factor(
    rep(c("north", "central", "south"), length.out = n),
    levels = c("north", "central", "south")
  )
  probability <- stats::plogis(
    -0.55 + 0.62 * x1 - 0.48 * x2 +
      0.28 * (region == "central") -
      0.19 * (region == "south")
  )
  uniform_grid <- ((id * 47) %% 181 + 0.5) / 181
  data.frame(
    y = as.numeric(uniform_grid < probability),
    x1 = x1,
    x2 = x2,
    region = region,
    ses = rep(seq_len(n / 3), each = 3),
    weight = 0.8 + (id %% 9) / 5
  )
}

phase8_estimate_vector <- function(fit) {
  c(
    fit$raw$estimates$detailed,
    `Total Index` = fit$raw$index$total,
    Residual = fit$raw$index$residual
  )
}

test_that("phase 8 OLS linearization is complete and preserves closure", {
  fit <- wvw_decomp(
    data = phase8_ols_data(),
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    weight_var = "weight",
    model_type = "ols",
    correction = "standard",
    vce_method = "linearized",
    quiet = TRUE
  )

  influence <- fit$raw$linearization$influence_function
  covariance <- fit$raw$linearization$vcov
  normalized_weight <- fit$raw$design$normalized_weights

  expect_true(fit$raw$vcov$detailed_complete)
  expect_true(fit$raw$vcov$grouped_complete)
  expect_true(fit$raw$vcov$overall_complete)
  expect_equal(
    unname(colSums(normalized_weight * influence)),
    rep(0, ncol(influence)),
    tolerance = 1e-13
  )
  expect_equal(covariance, t(covariance), tolerance = 1e-14)
  expect_gte(
    min(eigen(
      covariance,
      symmetric = TRUE,
      only.values = TRUE
    )$values),
    -1e-12
  )
  expect_equal(
    fit$raw$index$nonlinear_approximation_error,
    0,
    tolerance = 1e-12
  )
  expect_equal(
    unname(
      fit$raw$linearization$approximation_error_influence
    ),
    rep(0, nrow(influence)),
    tolerance = 1e-11
  )
  expect_equal(
    influence[, "Total Index"],
    rowSums(influence[, c("(Intercept)", "x1", "x2")]) +
      influence[, "Residual"],
    tolerance = 1e-11
  )
})

test_that("phase 8 reproduces the Stata conindex robust benchmark", {
  fit <- wvw_decomp(
    data = phase8_ols_data(),
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    weight_var = "weight",
    model_type = "ols",
    correction = "standard",
    vce_method = "linearized",
    quiet = TRUE
  )
  raw <- fit$raw
  rank <- raw$rank$fractional_rank
  weight <- as.numeric(raw$design$weights)
  outcome <- phase8_ols_data()$y[
    raw$sample$analytic_source_rows
  ]
  square_root_weight <- sqrt(weight)
  model_matrix <- cbind(
    intercept = square_root_weight,
    rank = rank * square_root_weight
  )
  response <- 2 * raw$index$rank_variance *
    raw$index$correction_scale * outcome *
    square_root_weight
  bread <- solve(crossprod(model_matrix))
  coefficients <- as.vector(
    bread %*% crossprod(model_matrix, response)
  )
  residual <- response -
    as.vector(model_matrix %*% coefficients)
  n <- length(outcome)
  robust_vcov <- n / (n - 2) * bread %*%
    crossprod(model_matrix * residual) %*% bread
  benchmark <-
    raw$linearization$conindex_robust_benchmark$total

  expect_equal(benchmark$estimate, coefficients[2], tolerance = 1e-14)
  expect_equal(
    benchmark$std_error,
    sqrt(robust_vcov[2, 2]),
    tolerance = 1e-14
  )
  expect_match(benchmark$finite_sample_correction, "N/\\(N-2\\)")
})

test_that("phase 8 nonlinear influence matches numerical contamination", {
  data <- phase8_binary_data()
  fit <- wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2", "region"),
    ses_var = "ses",
    weight_var = "weight",
    model_type = "logit",
    correction = "wagstaff",
    vce_method = "linearized",
    quiet = TRUE
  )

  sorted_position <- 109L
  source_position <-
    fit$raw$sample$analytic_source_rows[sorted_position]
  base_probability <- data$weight / sum(data$weight)
  epsilon <- 1e-5
  plus_probability <- (1 - epsilon) * base_probability
  minus_probability <- (1 + epsilon) * base_probability
  plus_probability[source_position] <-
    plus_probability[source_position] + epsilon
  minus_probability[source_position] <-
    minus_probability[source_position] - epsilon

  plus_data <- data
  minus_data <- data
  plus_data$weight <- plus_probability
  minus_data$weight <- minus_probability
  plus_fit <- wvw_decomp(
    plus_data,
    "y",
    c("x1", "x2", "region"),
    "ses",
    correction = "wagstaff",
    model_type = "logit",
    weight_var = "weight",
    quiet = TRUE
  )
  minus_fit <- wvw_decomp(
    minus_data,
    "y",
    c("x1", "x2", "region"),
    "ses",
    correction = "wagstaff",
    model_type = "logit",
    weight_var = "weight",
    quiet = TRUE
  )

  numerical_derivative <- (
    phase8_estimate_vector(plus_fit) -
      phase8_estimate_vector(minus_fit)
  ) / (2 * epsilon)
  analytic_influence <-
    fit$raw$linearization$influence_function[
      sorted_position,
      names(numerical_derivative)
    ]

  expect_equal(
    analytic_influence,
    numerical_derivative,
    tolerance = 2e-7
  )
  expect_true(all(is.finite(
    fit$results_overall$Std_Error
  )))
  expect_gt(
    abs(
      fit$raw$linearization$marginal_effect_gradient_beta[
        "x1",
        "x2"
      ]
    ),
    1e-8
  )
})

test_that("phase 8 conditions on a supplied precalculated rank", {
  data <- phase8_ols_data()
  data$rank_precalculated <- rank(data$ses) /
    (nrow(data) + 1)
  fit <- wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    precalc_rank_var = "rank_precalculated",
    weight_var = "weight",
    model_type = "ols",
    correction = "standard",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_identical(
    fit$raw$linearization$rank_treated_as,
    "fixed precalculated rank"
  )
  expect_match(
    fit$raw$known_limitations$precalculated_rank,
    "conditioned on"
  )
})
