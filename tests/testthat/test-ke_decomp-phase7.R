.phase7_test_data <- function(n = 30L) {
  source_id <- seq_len(n)
  data.frame(
    health = 0.31 +
      0.04 * sin(source_id / 3) +
      0.006 * source_id +
      0.025 * as.numeric(source_id %% 4 == 0),
    ses = rep(seq_len(n / 3), each = 3),
    income = 180 + 7 * source_id + 8 * cos(source_id / 4),
    x1 = sin(source_id / 5) + source_id / 45,
    x2 = as.numeric(source_id %% 4 == 0),
    category = factor(
      rep(c("base", "middle", "high"), length.out = n),
      levels = c("base", "middle", "high")
    )
  )
}

.phase7_grouped_rank <- function(ses, weights) {
  order_index <- order(ses, seq_along(ses))
  ordered_ses <- ses[order_index]
  normalized <- weights[order_index] / sum(weights)
  group_mass <- ave(normalized, ordered_ses, FUN = sum)
  first <- !duplicated(ordered_ses)
  cumulative <- cumsum(ifelse(first, group_mass, 0))
  ordered_rank <- cumulative - 0.5 * group_mass
  result <- numeric(length(ses))
  result[order_index] <- ordered_rank
  result
}

.phase7_manual_replicates <- function(
    fit, data, index_type, correction,
    dep_min = NA_real_, dep_max = NA_real_) {
  replication <- fit$raw$replication
  weights <- replication$replicate_weights
  ordered <- data[
    fit$raw$sample$analytic_source_rows,
    ,
    drop = FALSE
  ]
  matrix <- fit$raw$model$matrix

  result <- t(vapply(
    seq_len(ncol(weights)),
    function(replicate_number) {
      weight <- weights[, replicate_number]
      positive <- weight > 0
      normalized <- weight / sum(weight)
      health_mean <- sum(normalized * ordered$health)
      if (index_type == "rank") {
        position <- .phase7_grouped_rank(
          ordered$ses,
          weight
        )
        target_base <-
          2 * position * ordered$health - health_mean
      } else {
        income_mean <- sum(normalized * ordered$income)
        target_base <-
          ordered$income / income_mean * ordered$health -
          health_mean
      }
      scale <- switch(
        correction,
        generalized = 1,
        standard = 1 / health_mean,
        erreygers = if (index_type == "rank") {
          4 / (dep_max - dep_min)
        } else {
          1 / (dep_max - dep_min)
        },
        wagstaff = (dep_max - dep_min) / (
          (dep_max - health_mean) *
            (health_mean - dep_min)
        )
      )
      model <- stats::lm.wfit(
        x = matrix[positive, , drop = FALSE],
        y = target_base[positive] * scale,
        w = weight[positive]
      )
      coefficients <- model$coefficients
      coefficients[is.na(coefficients)] <- 0
      coefficients
    },
    numeric(ncol(matrix))
  ))
  colnames(result) <- colnames(matrix)
  result
}

.phase7_covariance <- function(estimates, method) {
  estimates <- estimates[
    stats::complete.cases(estimates),
    ,
    drop = FALSE
  ]
  scale <- if (method == "bootstrap") {
    1 / (nrow(estimates) - 1)
  } else {
    (nrow(estimates) - 1) / nrow(estimates)
  }
  scale * crossprod(
    sweep(estimates, 2, colMeans(estimates), "-")
  )
}

test_that("KE phase 7 ordinary bootstrap rebuilds the complete target", {
  data <- .phase7_test_data()
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = "x1 * category",
    ses_var = "ses",
    index_type = "rank",
    correction = "standard",
    vce_method = "bootstrap",
    boot_reps = 25L,
    seed = 71051L,
    quiet = TRUE
  ))
  replication <- fit$raw$replication
  manual <- .phase7_manual_replicates(
    fit, data, "rank", "standard"
  )

  expect_identical(replication$engine, "boot::boot")
  expect_true(replication$ordinary_bootstrap)
  expect_false(replication$ordinary_jackknife)
  expect_identical(replication$simulation, "ordinary")
  expect_identical(replication$statistic_type, "indices")
  expect_identical(replication$fit_weight_type, "frequency")
  expect_equal(
    colSums(replication$replicate_weights),
    rep(nrow(data), 25),
    tolerance = 0
  )
  expect_true(all(
    replication$replicate_weights >= 0 &
      replication$replicate_weights ==
        round(replication$replicate_weights)
  ))
  expect_equal(
    replication$estimates,
    manual,
    tolerance = 2e-13
  )
  expect_equal(
    replication$vcov,
    stats::cov(manual),
    tolerance = 2e-14
  )
  expect_identical(
    replication$center,
    colMeans(replication$valid_estimates)
  )
  expect_equal(replication$effective_scale, 1 / 24)
  expect_true(is.infinite(fit$summary_stats$DF))
  expect_true(is.na(fit$results_overall$Std_Error))
  expect_null(replication$index_variance)
  expect_identical(
    colnames(replication$estimates),
    colnames(fit$raw$model$matrix)
  )
})

test_that("KE phase 7 ordinary jackknife is deterministic delete-one", {
  data <- .phase7_test_data()
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "income",
    index_type = "level",
    correction = "erreygers",
    dep_min = 0,
    dep_max = 1,
    vce_method = "jackknife",
    quiet = TRUE
  ))
  replication <- fit$raw$replication
  manual <- .phase7_manual_replicates(
    fit, data, "level", "erreygers", 0, 1
  )
  expected_weights <- matrix(
    nrow(data) / (nrow(data) - 1),
    nrow = nrow(data),
    ncol = nrow(data)
  )
  diag(expected_weights) <- 0

  expect_identical(
    replication$engine,
    "resample::jackknife"
  )
  expect_true(replication$ordinary_jackknife)
  expect_false(replication$ordinary_bootstrap)
  expect_identical(replication$simulation, "delete-one")
  expect_identical(
    replication$statistic_type,
    "retained row indices"
  )
  expect_identical(
    replication$fit_weight_type,
    "JK1 replicate factor"
  )
  expect_equal(
    replication$replicate_weights,
    expected_weights,
    tolerance = 0
  )
  expect_equal(
    replication$estimates,
    manual,
    tolerance = 2e-13
  )
  expect_equal(
    replication$vcov,
    .phase7_covariance(manual, "jackknife"),
    tolerance = 2e-14
  )
  expect_equal(
    replication$effective_scale,
    (nrow(data) - 1) / nrow(data)
  )
  expect_equal(fit$summary_stats$DF, nrow(data) - 1)
  expect_true(is.na(fit$results_overall$Std_Error))
  expect_null(replication$index_variance)
})

test_that("KE phase 7 uses Stata valid-replicate scaling", {
  n <- 30L
  source_id <- seq_len(n)
  data <- data.frame(
    health_bootstrap = c(rep(1, 3), rep(0, n - 3)),
    health_jackknife = c(1, rep(0, n - 1)),
    ses = rep(seq_len(10), each = 3),
    x1 = sin(source_id / 3) + source_id / 40,
    x2 = as.numeric(source_id %% 4 == 0)
  )

  bootstrap <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health_bootstrap",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "wagstaff",
    dep_min = 0,
    dep_max = 1,
    vce_method = "bootstrap",
    boot_reps = 100L,
    seed = 1L,
    relax = TRUE,
    quiet = TRUE
  ))
  bootstrap_replication <- bootstrap$raw$replication
  expect_gt(bootstrap_replication$failed_replicates, 0)
  expect_equal(
    bootstrap_replication$effective_scale,
    1 / (bootstrap_replication$valid_replicates - 1)
  )
  expect_identical(
    bootstrap_replication$center,
    colMeans(bootstrap_replication$valid_estimates)
  )
  expect_equal(
    bootstrap_replication$vcov,
    stats::cov(bootstrap_replication$valid_estimates),
    tolerance = 2e-14
  )

  jackknife <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health_jackknife",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "standard",
    vce_method = "jackknife",
    relax = TRUE,
    quiet = TRUE
  ))
  jackknife_replication <- jackknife$raw$replication
  expect_equal(jackknife_replication$failed_replicates, 1)
  expect_equal(jackknife_replication$valid_replicates, n - 1)
  expect_equal(
    jackknife_replication$effective_scale,
    (jackknife_replication$valid_replicates - 1) /
      jackknife_replication$valid_replicates
  )
  expect_equal(
    jackknife_replication$vcov,
    .phase7_covariance(
      jackknife_replication$valid_estimates,
      "jackknife"
    ),
    tolerance = 2e-14
  )
  expect_match(
    jackknife$diagnostics$rep_adjust,
    "1 failed replicate"
  )
})

test_that("KE phase 7 obeys its caller RNG contract", {
  data <- .phase7_test_data()
  arguments <- list(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "bootstrap",
    boot_reps = 12L,
    seed = 71051L,
    quiet = TRUE
  )

  set.seed(77001)
  state_before <- .Random.seed
  first <- suppressWarnings(do.call(ke_decomp, arguments))
  expect_identical(.Random.seed, state_before)
  second <- suppressWarnings(do.call(ke_decomp, arguments))
  expect_identical(
    first$raw$replication$estimates,
    second$raw$replication$estimates
  )

  set.seed(77002)
  state_before <- .Random.seed
  arguments$seed <- NULL
  invisible(suppressWarnings(do.call(ke_decomp, arguments)))
  expect_false(identical(.Random.seed, state_before))

  set.seed(77003)
  state_before <- .Random.seed
  arguments$vce_method <- "jackknife"
  invisible(suppressWarnings(do.call(ke_decomp, arguments)))
  expect_identical(.Random.seed, state_before)
})

test_that("KE phase 7 keeps aliased replicate coefficients explicit", {
  data <- .phase7_test_data()
  data$x1_duplicate <- data$x1
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x1_duplicate", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "bootstrap",
    boot_reps = 12L,
    seed = 71051L,
    quiet = TRUE
  ))

  expect_true(fit$raw$model$aliased[["x1_duplicate"]])
  expect_equal(
    fit$raw$replication$estimates[, "x1_duplicate"],
    rep(0, 12),
    tolerance = 0
  )
  expect_equal(
    unname(fit$raw$model$vcov["x1_duplicate", ]),
    rep(0, ncol(fit$raw$model$vcov)),
    tolerance = 0
  )
})
