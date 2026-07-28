make_fairlie_phase5_data <- function() {
  set.seed(20260729)
  n0 <- 40L
  n1 <- 52L
  n <- n0 + n1
  id <- seq_len(n)
  group <- c(rep(0L, n0), rep(1L, n1))
  x1 <- as.numeric(scale(
    seq(-2.7, 2.7, length.out = n) +
      0.38 * sin(id / 5) - 0.50 * group
  ))
  x2 <- as.numeric(scale(
    cos(id / 7) + ((id * 13L) %% 31L) / 9 +
      0.33 * group
  ))
  x3 <- as.numeric(scale(
    sin(id / 10) - ((id * 11L) %% 23L) / 13 -
      0.24 * group
  ))
  eta <- -0.28 + 0.66 * x1 - 0.41 * x2 + 0.34 * x3 +
    0.44 * group
  weight <- 0.35 + ((id * 17L) %% 29L) / 7 +
    0.18 * (x1 - min(x1))
  data.frame(
    id = id,
    group = group,
    y_binary = rbinom(n, size = 1L, prob = plogis(eta)),
    x1 = x1,
    x2 = x2,
    x3 = x3,
    weight = weight,
    constant_weight = 2
  )
}

fit_fairlie_phase5 <- function(
    data,
    weight_var = "weight",
    reps = 83L,
    seed = 510103L,
    reference = "pooled"
) {
  f_decomp(
    data = data,
    dep_var = "y_binary",
    group_var = "group",
    group_levels = c(0, 1),
    indep_vars = c("x1", "x2", "x3"),
    ref_method = reference,
    pooled_anchor = "favored",
    model_type = "logit",
    randomize_order = FALSE,
    weight_var = weight_var,
    vce_method = "linearized",
    reps = reps,
    seed = seed,
    quiet = TRUE
  )
}

manual_fairlie_phase5 <- function(data, fit, reps, seed) {
  source_rows <- fit$raw$analytic_sample
  working <- data[source_rows, , drop = FALSE]
  x <- stats::model.matrix(~x1 + x2 + x3, data = working)
  beta <- fit$raw$reference_coefficients
  x0 <- x[working$group == 0, , drop = FALSE]
  x1 <- x[working$group == 1, , drop = FALSE]
  w0 <- working$weight[working$group == 0]
  w1 <- working$weight[working$group == 1]
  order0 <- order(stats::plogis(as.vector(x0 %*% beta)))
  order1 <- order(stats::plogis(as.vector(x1 %*% beta)))
  x0 <- x0[order0, , drop = FALSE]
  x1 <- x1[order1, , drop = FALSE]
  w0 <- w0[order0]
  w1 <- w1[order1]
  n_match <- trunc((nrow(x0) + nrow(x1)) / 2)
  contributions <- matrix(NA_real_, nrow = reps, ncol = 3)

  sample_pps <- function(weights) {
    upper_bound <- cumsum(weights) / sum(weights)
    uniforms <- sort(stats::runif(n_match))
    findInterval(
      uniforms,
      upper_bound,
      left.open = TRUE
    ) + 1L
  }

  set.seed(seed)
  for (replication in seq_len(reps)) {
    selected0 <- sample_pps(w0)
    selected1 <- sample_pps(w1)
    left <- x0[selected0, , drop = FALSE]
    matched1 <- x1[selected1, , drop = FALSE]
    for (block in seq_len(3)) {
      right <- left
      right[, block + 1] <- matched1[, block + 1]
      contributions[replication, block] <- mean(
        stats::plogis(as.vector(left %*% beta)) -
          stats::plogis(as.vector(right %*% beta))
      )
      left <- right
    }
  }
  colnames(contributions) <- c("Exp_x1", "Exp_x2", "Exp_x3")
  list(
    estimate = colMeans(contributions),
    monte_carlo_se = apply(contributions, 2, stats::sd) / sqrt(reps)
  )
}

test_that("weighted Fairlie matching reproduces ranked PPS draws", {
  data <- make_fairlie_phase5_data()
  reps <- 83L
  seed <- 510103L
  fit <- fit_fairlie_phase5(data, reps = reps, seed = seed)
  manual <- manual_fairlie_phase5(data, fit, reps, seed)

  expect_equal(
    fit$raw$detailed_estimates,
    manual$estimate,
    tolerance = 1e-14
  )
  expect_equal(
    fit$raw$matching$detailed_monte_carlo_se,
    manual$monte_carlo_se,
    tolerance = 1e-14
  )
  expect_equal(fit$raw$matching$sample_size, 46)
  expect_identical(
    fit$raw$matching$sampling_method,
    "ranked_pps_with_replacement"
  )
  expect_identical(fit$raw$matching$sampled_group, "both")
  expect_true(fit$raw$matching$stochastic_subsampling)
  expect_identical(fit$raw$weighting$interpretation, "pweight")
})

test_that("constant pweights preserve unweighted model point estimands", {
  data <- make_fairlie_phase5_data()
  weighted <- fit_fairlie_phase5(
    data,
    weight_var = "constant_weight",
    reps = 17L,
    reference = "group0"
  )
  unweighted <- fit_fairlie_phase5(
    data,
    weight_var = NULL,
    reps = 1L,
    reference = "group0"
  )

  expect_equal(
    weighted$raw$overall_estimates,
    unweighted$raw$overall_estimates,
    tolerance = 0
  )
  expect_equal(
    weighted$raw$reference_coefficients,
    unweighted$raw$reference_coefficients,
    tolerance = 0
  )
  expect_equal(
    weighted$raw$reference_predictions,
    unweighted$raw$reference_predictions,
    tolerance = 0
  )
})

test_that("global pweight rescaling leaves every point estimate unchanged", {
  data <- make_fairlie_phase5_data()
  data$rescaled_weight <- 1000 * data$weight
  original <- fit_fairlie_phase5(data, weight_var = "weight")
  rescaled <- fit_fairlie_phase5(
    data,
    weight_var = "rescaled_weight"
  )

  expect_equal(
    original$raw$overall_estimates,
    rescaled$raw$overall_estimates,
    tolerance = 2e-13
  )
  expect_equal(
    original$raw$detailed_estimates,
    rescaled$raw$detailed_estimates,
    tolerance = 2e-13
  )
  expect_equal(
    original$raw$reference_model_coefficients,
    rescaled$raw$reference_model_coefficients,
    tolerance = 2e-13
  )
  expect_equal(
    original$raw$reference_predictions,
    rescaled$raw$reference_predictions,
    tolerance = 2e-13
  )
  expect_true(
    original$raw$weighting$scale_invariant_point_estimator
  )
})

test_that("zero weights are excluded and recorded like Stata pweights", {
  data <- make_fairlie_phase5_data()
  zero_rows <- c(3L, 47L)
  data$weight[zero_rows] <- 0
  fit <- fit_fairlie_phase5(data)

  expect_equal(fit$summary_stats$N, nrow(data) - 2L)
  expect_equal(fit$summary_stats$N0, 39)
  expect_equal(fit$summary_stats$N1, 51)
  expect_equal(fit$summary_stats$n_zero_weight_dropped, 2)
  expect_equal(fit$raw$weighting$zero_weight_rows_excluded, 2)
  expect_false(any(zero_rows %in% fit$raw$analytic_sample))
  expect_equal(
    as.numeric(fit$raw$weighting$group_total_weight),
    as.numeric(c(
      sum(data$weight[data$group == 0]),
      sum(data$weight[data$group == 1])
    )),
    tolerance = 1e-12
  )
})

test_that("invalid pweights and a zero-weight group fail explicitly", {
  data <- make_fairlie_phase5_data()

  negative <- data
  negative$weight[[1]] <- -1
  expect_error(
    fit_fairlie_phase5(negative),
    "cannot contain negative"
  )

  nonfinite <- data
  nonfinite$weight[[1]] <- Inf
  expect_error(
    fit_fairlie_phase5(nonfinite),
    "only finite"
  )

  nonnumeric <- data
  nonnumeric$weight <- as.character(nonnumeric$weight)
  nonnumeric$weight[[1]] <- "not-a-weight"
  expect_error(
    fit_fairlie_phase5(nonnumeric),
    "must contain numeric"
  )

  zero_group <- data
  zero_group$weight[zero_group$group == 0] <- 0
  expect_error(
    fit_fairlie_phase5(zero_group),
    "both selected groups.*positive weight"
  )
})

test_that("missing pweights follow complete-case exclusion", {
  data <- make_fairlie_phase5_data()
  data$weight[[1]] <- NA_real_
  fit <- fit_fairlie_phase5(data)

  expect_equal(fit$summary_stats$N, nrow(data) - 1L)
  expect_equal(fit$summary_stats$n_missing_dropped, 1)
  expect_equal(fit$summary_stats$n_zero_weight_dropped, 0)
  expect_false(1L %in% fit$raw$analytic_sample)
})
