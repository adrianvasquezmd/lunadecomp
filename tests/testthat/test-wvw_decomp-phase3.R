.phase3_wvw_data <- function() {
  id <- seq_len(60L)
  data <- data.frame(
    id = id,
    ses = ((17L * id) %% 61L) + 1L,
    x1 = 1 + ((7L * id) %% 13L) / 4,
    x2 = 2 + ((11L * id) %% 17L) / 5,
    x3 = as.integer((id %% 4L) %in% c(0L, 1L)),
    weight = 1L + ((5L * id) %% 7L)
  )
  data$noise <- (((19L * id) %% 23L) - 11) / 20
  data$y <- with(
    data,
    5 + 0.7 * x1 - 0.45 * x2 + 0.85 * x3 + noise
  )
  data
}

.phase3_wvw_fit <- function(data, determinants = c("x1", "x2", "x3")) {
  wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = determinants,
    ses_var = "ses",
    correction = "standard",
    model_type = "ols",
    use_svy = FALSE,
    weight_var = "weight",
    vce_method = "linearized",
    quiet = TRUE
  )
}

test_that("phase 3 weighted OLS core matches the frozen Stata reference", {
  data <- .phase3_wvw_data()
  fit <- .phase3_wvw_fit(data)
  raw <- fit$raw
  determinant <- raw$determinants[
    match(c("x1", "x2", "x3"), raw$determinants$column),
  ]

  expect_equal(raw$sample$analytic_n, 60)
  expect_equal(raw$design$population_weight, 243)
  expect_equal(
    raw$index$total,
    -0.0070343621868394924,
    tolerance = 1e-12
  )
  expect_equal(
    raw$index$explained,
    -0.00081637924749773877,
    tolerance = 1e-12
  )
  expect_equal(
    raw$index$residual,
    -0.0062179829393417541,
    tolerance = 1e-12
  )
  expect_equal(
    unname(raw$model$coefficients),
    c(
      5.0224672072861809,
      0.62535675323223106,
      -0.4169887735900476,
      0.94206835413112011
    ),
    tolerance = 1e-12
  )
  expect_equal(
    determinant$contribution,
    c(
      -0.0030987988037449272,
      -0.00097559884464647108,
      0.0032580184008936595
    ),
    tolerance = 1e-12
  )
})

test_that("phase 3 probability-weighted points are scale invariant", {
  data <- .phase3_wvw_data()
  fit_1 <- .phase3_wvw_fit(data)$raw
  data$weight <- 10 * data$weight
  fit_10 <- .phase3_wvw_fit(data)$raw

  expect_equal(
    fit_10$rank$fractional_rank,
    fit_1$rank$fractional_rank,
    tolerance = 1e-15
  )
  expect_equal(
    fit_10$model$coefficients,
    fit_1$model$coefficients,
    tolerance = 1e-12
  )
  expect_equal(
    c(
      fit_10$index$total,
      fit_10$index$explained,
      fit_10$index$residual
    ),
    c(
      fit_1$index$total,
      fit_1$index$explained,
      fit_1$index$residual
    ),
    tolerance = 1e-12
  )
  expect_equal(
    fit_10$estimates$detailed,
    fit_1$estimates$detailed,
    tolerance = 1e-12
  )
})

test_that("phase 3 distinguishes population and fweight denominators", {
  skip_if_not_installed("rineq")
  data <- .phase3_wvw_data()
  fit <- .phase3_wvw_fit(data)$raw
  weight_normalized <- data$weight / sum(data$weight)
  fractional_rank <- rineq::rank_gwt(data$ses, data$weight)
  outcome_mean <- stats::weighted.mean(data$y, data$weight)
  population_covariance <- sum(
    weight_normalized *
      (data$y - outcome_mean) *
      (
        fractional_rank -
          stats::weighted.mean(fractional_rank, data$weight)
      )
  )
  fweight_factor <- sum(data$weight) / (sum(data$weight) - 1)
  fweight_total <- 2 * population_covariance *
    fweight_factor / outcome_mean

  expect_equal(
    fit$index$total,
    2 * population_covariance / outcome_mean,
    tolerance = 1e-12
  )
  expect_equal(
    fweight_total / fit$index$total,
    fweight_factor,
    tolerance = 1e-12
  )

  data_10 <- transform(data, weight = 10 * weight)
  factor_10 <- sum(data_10$weight) /
    (sum(data_10$weight) - 1)
  fweight_total_10 <- fit$index$total * factor_10
  expect_gt(abs(fweight_total_10 - fweight_total), 1e-6)
})

test_that("phase 3 weighted factor blocks are scale invariant", {
  data <- .phase3_wvw_data()
  data$region <- factor(1L + ((3L * data$id) %% 5L))
  determinants <- c("x1", "x2", "x3", "region")
  fit_1 <- .phase3_wvw_fit(data, determinants)$raw
  data$weight <- 10 * data$weight
  fit_10 <- .phase3_wvw_fit(data, determinants)$raw

  region_1 <- fit_1$determinants$source_term == "region"
  region_10 <- fit_10$determinants$source_term == "region"

  expect_equal(sum(region_1), 4)
  expect_equal(sum(region_10), 4)
  expect_equal(
    fit_10$determinants$contribution[region_10],
    fit_1$determinants$contribution[region_1],
    tolerance = 1e-12
  )
  expect_equal(
    sum(fit_10$determinants$contribution[region_10]),
    sum(fit_1$determinants$contribution[region_1]),
    tolerance = 1e-12
  )
  expect_equal(
    fit_10$index$total,
    fit_1$index$total,
    tolerance = 1e-12
  )
})
