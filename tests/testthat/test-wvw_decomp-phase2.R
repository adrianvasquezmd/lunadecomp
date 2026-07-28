.phase2_wvw_data <- function() {
  id <- seq_len(60L)
  data <- data.frame(
    id = id,
    ses = ((17L * id) %% 61L) + 1L,
    x1 = 1 + ((7L * id) %% 13L) / 4,
    x2 = 2 + ((11L * id) %% 17L) / 5,
    x3 = as.integer((id %% 4L) %in% c(0L, 1L))
  )
  data$noise <- (((19L * id) %% 23L) - 11) / 20
  data$y <- with(
    data,
    5 + 0.7 * x1 - 0.45 * x2 + 0.85 * x3 + noise
  )
  data
}

.phase2_wvw_fit <- function(data) {
  wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2", "x3"),
    ses_var = "ses",
    correction = "standard",
    model_type = "ols",
    use_svy = FALSE,
    vce_method = "linearized",
    quiet = TRUE
  )
}

test_that("phase 2 OLS core matches the frozen Stata reference", {
  data <- .phase2_wvw_data()
  fit <- .phase2_wvw_fit(data)
  raw <- fit$raw
  determinant <- raw$determinants[
    match(c("x1", "x2", "x3"), raw$determinants$column),
  ]

  expect_equal(raw$sample$analytic_n, 60)
  expect_equal(
    raw$index$total,
    -0.0025382262996941911,
    tolerance = 1e-12
  )
  expect_equal(
    raw$index$explained,
    0.0007576640362588398,
    tolerance = 1e-12
  )
  expect_equal(
    raw$index$residual,
    -0.0032958903359530272,
    tolerance = 1e-12
  )
  expect_equal(
    unname(raw$model$coefficients),
    c(
      4.8645570107205414,
      0.67596370171523912,
      -0.39012992671058955,
      0.84321472502720007
    ),
    tolerance = 1e-12
  )
  expect_equal(
    determinant$contribution,
    c(
      -0.00027110390603709401,
      0.0010287679422959347,
      -8.4477942766867291e-19
    ),
    tolerance = 1e-12
  )

  expected_rank <- (rank(data$ses, ties.method = "average") - 0.5) /
    nrow(data)
  rank_by_source_row <- raw$rank$fractional_rank[
    match(seq_len(nrow(data)), raw$rank$source_row)
  ]
  expect_equal(rank_by_source_row, expected_rank, tolerance = 2e-15)
})

test_that("phase 2 OLS components reproduce the algebraic decomposition", {
  data <- .phase2_wvw_data()
  fit <- .phase2_wvw_fit(data)
  raw <- fit$raw
  matrix_x <- stats::model.matrix(~ x1 + x2 + x3, data)
  model <- stats::lm(y ~ x1 + x2 + x3, data = data)
  fractional_rank <- (rank(data$ses) - 0.5) / nrow(data)
  outcome_mean <- mean(data$y)
  rank_mean <- mean(fractional_rank)
  gc_x <- 2 * (
    colMeans(matrix_x * fractional_rank) -
      colMeans(matrix_x) * rank_mean
  )
  manual_contribution <-
    stats::coef(model) * gc_x / outcome_mean
  manual_contribution["(Intercept)"] <- 0
  manual_total <- 2 * mean(
    (data$y - outcome_mean) *
      (fractional_rank - rank_mean)
  ) / outcome_mean
  model_residual <- stats::residuals(model)
  manual_residual <- 2 * mean(
    (model_residual - mean(model_residual)) *
      (fractional_rank - rank_mean)
  ) / outcome_mean

  expect_equal(
    raw$estimates$detailed[names(manual_contribution)],
    manual_contribution,
    tolerance = 1e-12
  )
  expect_equal(raw$index$total, manual_total, tolerance = 1e-12)
  expect_equal(raw$index$residual, manual_residual, tolerance = 1e-12)
  expect_equal(
    raw$index$total,
    raw$index$explained + raw$index$residual,
    tolerance = 1e-12
  )
})

test_that("phase 2 OLS decomposition is invariant to order and scale", {
  data <- .phase2_wvw_data()
  baseline <- .phase2_wvw_fit(data)$raw

  variants <- list(
    reverse_rows = data[rev(seq_len(nrow(data))), , drop = FALSE],
    x1_times_100 = transform(data, x1 = 100 * x1),
    y_times_10 = transform(data, y = 10 * y),
    ses_affine = transform(data, ses = 10 + 7 * ses)
  )
  variant_fits <- lapply(variants, function(value) {
    rownames(value) <- NULL
    .phase2_wvw_fit(value)$raw
  })

  baseline_overall <- c(
    baseline$index$total,
    baseline$index$explained,
    baseline$index$residual
  )
  baseline_detailed <- baseline$estimates$detailed[
    c("x1", "x2", "x3")
  ]

  for (name in names(variant_fits)) {
    current <- variant_fits[[name]]
    expect_equal(
      c(
        current$index$total,
        current$index$explained,
        current$index$residual
      ),
      baseline_overall,
      tolerance = 1e-12,
      info = name
    )
    expect_equal(
      current$estimates$detailed[c("x1", "x2", "x3")],
      baseline_detailed,
      tolerance = 1e-12,
      info = name
    )
  }

  expect_equal(
    variant_fits$x1_times_100$model$coefficients["x1"] * 100,
    baseline$model$coefficients["x1"],
    tolerance = 1e-12
  )
  expect_equal(
    variant_fits$y_times_10$model$coefficients / 10,
    baseline$model$coefficients,
    tolerance = 1e-12
  )
})

test_that("phase 2 OLS point decomposition is triangulated with rineq", {
  skip_if_not_installed("rineq")
  data <- .phase2_wvw_data()
  fit <- .phase2_wvw_fit(data)
  model <- stats::lm(y ~ x1 + x2 + x3, data = data)
  rineq_decomposition <- rineq::contribution(
    model,
    ranker = data$ses,
    correction = FALSE,
    type = "CI",
    intercept = "exclude"
  )
  rineq_total <- rineq::ci(
    ineqvar = data$ses,
    outcome = data$y,
    method = "direct",
    rank_function = rineq::rank_gwt
  )

  expect_equal(
    fit$raw$index$total,
    rineq_total$concentration_index,
    tolerance = 1e-12
  )
  expect_equal(
    fit$raw$index$total,
    rineq_decomposition$overall_ci$concentration_index,
    tolerance = 1e-12
  )
  expect_equal(
    fit$raw$estimates$detailed[c("x1", "x2", "x3")],
    rineq_decomposition$ci_contribution[c("x1", "x2", "x3")],
    tolerance = 1e-12
  )
  expect_equal(
    fit$raw$index$residual,
    unname(rineq_decomposition$ci_contribution["residual"]),
    tolerance = 1e-12
  )
})
