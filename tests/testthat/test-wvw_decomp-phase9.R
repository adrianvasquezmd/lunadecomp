phase9_ordinary_data <- function(n = 36L) {
  id <- seq_len(n)
  x1 <- sin(id / 5) + id / 30
  x2 <- as.numeric(id %% 3 == 0)
  data.frame(
    y = 2.8 + 0.54 * x1 - 0.31 * x2 +
      0.15 * cos(id / 4),
    x1 = x1,
    x2 = x2,
    ses = rep(seq_len(n / 2), each = 2)
  )
}

phase9_binary_data <- function(n = 120L) {
  id <- seq_len(n)
  x1 <- sin(id / 8) + (id %% 11) / 12
  x2 <- as.numeric(id %% 4 < 2)
  probability <- stats::plogis(
    -0.45 + 0.58 * x1 - 0.41 * x2
  )
  uniform_grid <- ((id * 43) %% 127 + 0.5) / 127
  data.frame(
    y = as.numeric(uniform_grid < probability),
    x1 = x1,
    x2 = x2,
    ses = rep(seq_len(n / 3), each = 3)
  )
}

phase9_reconstruct_vcov <- function(replication) {
  valid_estimates <- replication$estimates[
    replication$valid,
    ,
    drop = FALSE
  ]
  centered <- sweep(
    valid_estimates,
    2,
    colMeans(valid_estimates),
    "-"
  )
  replication$effective_scale * crossprod(
    centered * sqrt(replication$valid_rscales)
  )
}

test_that("phase 9 ordinary jackknife matches Stata non-MSE conventions", {
  data <- phase9_ordinary_data()
  fit <- wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    vce_method = "jackknife",
    quiet = TRUE
  )
  replication <- fit$raw$replication
  n <- nrow(data)
  weights <- replication$replicate_weights

  expect_true(replication$ordinary_resampling)
  expect_identical(replication$generator_type,
                   "ordinary delete-one jackknife")
  expect_true(replication$ordinary_jackknife)
  expect_false(replication$ordinary_bootstrap)
  expect_identical(replication$engine, "resample::jackknife")
  expect_identical(replication$simulation, "delete-one")
  expect_identical(
    replication$statistic_type,
    "retained row indices"
  )
  expect_identical(
    replication$fit_weight_type,
    "JK1 replicate factor"
  )
  expect_identical(replication$replicates, replication$estimates)
  expect_false(replication$mse)
  expect_identical(replication$center_method, "replicate mean")
  expect_equal(replication$center,
               colMeans(replication$estimates),
               tolerance = 0)
  expect_equal(replication$requested_replicates, n)
  expect_equal(replication$valid_replicates, n)
  expect_equal(replication$failed_replicates, 0)
  expect_equal(replication$generator_scale, (n - 1) / n,
               tolerance = 0)
  expect_equal(replication$effective_scale, (n - 1) / n,
               tolerance = 0)
  expect_equal(dim(weights), c(n, n))
  expect_equal(diag(weights), rep(0, n), tolerance = 0)
  expect_equal(
    weights[row(weights) != col(weights)],
    rep(n / (n - 1), n * (n - 1)),
    tolerance = 0
  )
  expect_equal(
    unname(colSums(weights)),
    rep(n, n),
    tolerance = 1e-14
  )
  expect_equal(
    replication$vcov,
    phase9_reconstruct_vcov(replication),
    tolerance = 1e-15
  )
  expect_equal(fit$summary_stats$DF, n - 1)
  expect_false(replication$pending_methodological_review)
  expect_null(fit$raw$known_limitations$replication)
})

test_that("phase 9 bootstrap uses ordinary multinomial samples", {
  data <- phase9_ordinary_data()
  fit <- wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    vce_method = "bootstrap",
    boot_reps = 30,
    seed = 9021,
    quiet = TRUE
  )
  repeated_fit <- wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    vce_method = "bootstrap",
    boot_reps = 30,
    seed = 9021,
    quiet = TRUE
  )
  replication <- fit$raw$replication
  weights <- replication$replicate_weights

  expect_true(replication$ordinary_resampling)
  expect_match(replication$generator_type, "multinomial bootstrap")
  expect_true(replication$ordinary_bootstrap)
  expect_false(replication$ordinary_jackknife)
  expect_identical(replication$engine, "boot::boot")
  expect_identical(replication$simulation, "ordinary")
  expect_identical(replication$statistic_type, "indices")
  expect_identical(replication$fit_weight_type, "frequency")
  expect_identical(replication$replicates, replication$estimates)
  expect_false(replication$mse)
  expect_identical(replication$center_method, "replicate mean")
  expect_equal(replication$center,
               colMeans(replication$estimates),
               tolerance = 0)
  expect_equal(replication$requested_replicates, 30)
  expect_equal(replication$valid_replicates, 30)
  expect_equal(replication$failed_replicates, 0)
  expect_equal(replication$generator_scale, 1 / 29,
               tolerance = 0)
  expect_equal(replication$effective_scale, 1 / 29,
               tolerance = 0)
  expect_equal(dim(weights), c(nrow(data), 30))
  expect_equal(weights, round(weights), tolerance = 0)
  expect_true(all(weights >= 0))
  expect_equal(
    unname(colSums(weights)),
    rep(nrow(data), 30),
    tolerance = 0
  )
  expect_equal(
    replication$vcov,
    phase9_reconstruct_vcov(replication),
    tolerance = 1e-15
  )
  expect_equal(
    weights,
    repeated_fit$raw$replication$replicate_weights,
    tolerance = 0
  )
  expect_equal(fit$summary_stats$DF, Inf)
  expect_false(replication$pending_methodological_review)
})

test_that("phase 9 restores the caller RNG state", {
  data <- phase9_ordinary_data()
  set.seed(193)
  before <- .Random.seed

  invisible(wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    vce_method = "bootstrap",
    boot_reps = 8,
    seed = 771,
    quiet = TRUE
  ))

  expect_identical(.Random.seed, before)
})

test_that("phase 9 jackknife does not initialize an absent RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv)
  previous_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  invisible(wvw_decomp(
    data = phase9_ordinary_data(),
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    vce_method = "jackknife",
    quiet = TRUE
  ))

  expect_false(exists(".Random.seed", envir = .GlobalEnv))
})

test_that("phase 9 recomputes nonlinear AMEs, ranks, and correction", {
  fit <- wvw_decomp(
    data = phase9_binary_data(),
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    correction = "wagstaff",
    model_type = "logit",
    vce_method = "bootstrap",
    boot_reps = 20,
    seed = 482,
    quiet = TRUE
  )
  replication <- fit$raw$replication

  expect_equal(replication$valid_replicates, 20)
  expect_equal(replication$failed_replicates, 0)
  expect_true(all(is.finite(replication$estimates)))
  expect_true(all(is.finite(replication$vcov)))
  expect_equal(
    replication$vcov,
    phase9_reconstruct_vcov(replication),
    tolerance = 1e-14
  )
  expect_true(fit$raw$vcov$detailed_complete)
  expect_true(fit$raw$vcov$grouped_complete)
  expect_true(fit$raw$vcov$overall_complete)
})

test_that("phase 9 uses Stata's valid-replicate scale after failures", {
  x <- seq(-2, 2, length.out = 25)
  data <- data.frame(
    y = as.numeric(x > 0),
    x = x,
    ses = seq_along(x)
  )
  data$y[3] <- 1

  fit <- wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = "x",
    ses_var = "ses",
    correction = "generalized",
    model_type = "logit",
    vce_method = "jackknife",
    quiet = TRUE
  )
  replication <- fit$raw$replication
  valid_count <- replication$valid_replicates

  expect_gt(replication$failed_replicates, 0)
  expect_equal(
    replication$effective_scale,
    (valid_count - 1) / valid_count,
    tolerance = 0
  )
  expect_equal(fit$summary_stats$DF, valid_count - 1)
  expect_equal(
    replication$vcov,
    phase9_reconstruct_vcov(replication),
    tolerance = 1e-14
  )
})

test_that("phase 9 rejects invalid bootstrap replication counts", {
  data <- phase9_ordinary_data()
  for (bad_reps in list(1, 2.5, NA_real_, Inf, "20")) {
    expect_error(
      wvw_decomp(
        data = data,
        dep_var = "y",
        indep_vars = c("x1", "x2"),
        ses_var = "ses",
        vce_method = "bootstrap",
        boot_reps = bad_reps,
        quiet = TRUE
      ),
      "boot_reps"
    )
  }
})
