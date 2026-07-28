make_fairlie_phase3_data <- function() {
  set.seed(20260727)
  n_per_group <- 90L
  n <- 2L * n_per_group
  id <- seq_len(n)
  group <- rep(0:1, each = n_per_group)
  within_group <- ave(id, group, FUN = seq_along)
  x1 <- as.numeric(scale(
    seq(-2.8, 2.8, length.out = n) +
      0.42 * sin(id / 5) - 0.55 * group
  ))
  x2 <- as.numeric(scale(
    cos(id / 8) + ((id * 17L) %% 37L) / 11 +
      0.35 * group
  ))
  x3 <- as.numeric(scale(
    sin(id / 11) - ((id * 11L) %% 29L) / 15 -
      0.28 * group
  ))
  eta <- -0.20 + 0.62 * x1 - 0.38 * x2 + 0.31 * x3 +
    0.48 * group
  generated <- data.frame(
    id = id,
    group = group,
    y_binary = rbinom(n, size = 1L, prob = plogis(eta)),
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
  generated[
    (group == 0 & within_group <= 48) |
      (group == 1 & within_group <= 84),
    ,
    drop = FALSE
  ]
}

fit_fairlie_phase3 <- function(
    data,
    reps = 257,
    seed = 310097
) {
  f_decomp(
    data = data,
    dep_var = "y_binary",
    group_var = "group",
    group_levels = c(0, 1),
    indep_vars = c("x1", "x2", "x3"),
    ref_method = "group0",
    model_type = "logit",
    randomize_order = FALSE,
    vce_method = "linearized",
    reps = reps,
    seed = seed,
    quiet = TRUE
  )
}

manual_fairlie_phase3 <- function(data, fit, reps, seed) {
  x <- stats::model.matrix(~x1 + x2 + x3, data = data)
  x0 <- x[data$group == 0, , drop = FALSE]
  x1 <- x[data$group == 1, , drop = FALSE]
  beta <- fit$raw$reference_coefficients
  order0 <- order(stats::plogis(as.vector(x0 %*% beta)))
  order1 <- order(stats::plogis(as.vector(x1 %*% beta)))
  n_match <- nrow(x0)
  contributions <- matrix(NA_real_, nrow = reps, ncol = 3)

  set.seed(seed)
  for (replication in seq_len(reps)) {
    selected1 <- order1[sort(sample.int(
      nrow(x1),
      n_match,
      replace = FALSE
    ))]
    left <- x0[order0, , drop = FALSE]
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

test_that("unequal matching samples only the larger ranked group", {
  data <- make_fairlie_phase3_data()
  fit <- fit_fairlie_phase3(data)

  expect_equal(fit$summary_stats$N0, 48)
  expect_equal(fit$summary_stats$N1, 84)
  expect_equal(fit$raw$matching$sample_size, 48)
  expect_identical(
    fit$raw$matching$sampling_method,
    "ranked_srs_without_replacement"
  )
  expect_identical(fit$raw$matching$sampled_group, "group_1")
  expect_true(fit$raw$matching$stochastic_subsampling)
  expect_identical(
    unname(fit$raw$matching$group_sizes),
    c(48L, 84L)
  )
})

test_that("Phase 3 matching and Monte Carlo SE reproduce a manual harness", {
  data <- make_fairlie_phase3_data()
  reps <- 257L
  seed <- 310097L
  fit <- fit_fairlie_phase3(data, reps = reps, seed = seed)
  manual <- manual_fairlie_phase3(data, fit, reps, seed)

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
  expect_true(all(
    fit$raw$matching$detailed_monte_carlo_se > 0
  ))
})

test_that("the complete smaller group no longer consumes random draws", {
  data <- make_fairlie_phase3_data()
  reps <- 11L

  set.seed(81231)
  fit_fairlie_phase3(data, reps = reps, seed = NULL)
  observed_next_random <- stats::runif(1)

  set.seed(81231)
  for (replication in seq_len(reps)) {
    sort(sample.int(84, 48, replace = FALSE))
  }
  expected_next_random <- stats::runif(1)

  expect_identical(observed_next_random, expected_next_random)
})

test_that("explicit Phase 3 seeds are reproducible and restore caller RNG", {
  data <- make_fairlie_phase3_data()

  set.seed(99173)
  caller_state <- .Random.seed
  fit1 <- fit_fairlie_phase3(data, reps = 37, seed = 7712)
  expect_identical(.Random.seed, caller_state)
  fit2 <- fit_fairlie_phase3(data, reps = 37, seed = 7712)

  expect_equal(
    fit1$raw$detailed_estimates,
    fit2$raw$detailed_estimates,
    tolerance = 0
  )
  expect_equal(
    fit1$raw$matching$detailed_monte_carlo_se,
    fit2$raw$matching$detailed_monte_carlo_se,
    tolerance = 0
  )
})

test_that("one matching replication has undefined Monte Carlo error", {
  data <- make_fairlie_phase3_data()
  fit <- fit_fairlie_phase3(data, reps = 1)

  expect_true(all(is.na(
    fit$raw$matching$detailed_monte_carlo_se
  )))
})
