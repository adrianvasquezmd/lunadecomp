.phase10_wvw_data <- function(cluster_size = 8L) {
  n <- 4L * 4L * cluster_size
  id <- seq_len(n)
  strata <- rep(seq_len(4L), each = 4L * cluster_size)
  psu_within_stratum <- rep(
    rep(seq_len(4L), each = cluster_size),
    times = 4L
  )
  psu <- (strata - 1L) * 4L + psu_within_stratum
  x1 <- sin(id / 9) + (id %% 13L) / 10 -
    0.08 * strata
  x2 <- as.numeric((id * 5L + strata) %% 9L < 4L)

  data.frame(
    id = id,
    strata = strata,
    psu = psu,
    psu_reused = psu_within_stratum,
    ses = rep(seq_len(n / 2L), each = 2L),
    weight = 0.65 + (id %% 11L) / 5 +
      0.07 * strata + 0.03 * psu_within_stratum,
    x1 = x1,
    x2 = x2,
    y = 4.5 + 0.68 * x1 - 0.47 * x2 +
      0.16 * cos(id / 6) + 0.05 * strata
  )
}

.phase10_wvw_fit <- function(
    data = .phase10_wvw_data(4L),
    method = "linearized",
    psu = "psu",
    boot_reps = 12L,
    seed = 104732L,
    lonely_psu = "adjust",
    bootstrap_singleton = "fail"
) {
  wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    weight_var = "weight",
    strata_var = "strata",
    psu_var = psu,
    vce_method = method,
    boot_reps = boot_reps,
    seed = seed,
    lonely_psu = lonely_psu,
    bootstrap_singleton = bootstrap_singleton,
    quiet = TRUE
  )
}

.phase10_reconstruct_survey_vcov <- function(replication) {
  valid <- replication$valid
  estimates <- replication$estimates[valid, , drop = FALSE]
  centers <- replication$center
  if (is.vector(centers)) {
    centers <- matrix(
      centers,
      nrow = nrow(estimates),
      ncol = ncol(estimates),
      byrow = TRUE
    )
  }
  replication$effective_scale * crossprod(
    (estimates - centers) *
      sqrt(replication$valid_rscales)
  )
}

test_that("phase 10 design variables activate survey mode", {
  data <- .phase10_wvw_data(4L)
  fit <- wvw_decomp(
    data,
    "y",
    c("x1", "x2"),
    "ses",
    strata_var = "strata",
    psu_var = "psu",
    quiet = TRUE
  )

  expect_true(fit$raw$settings$use_svy)
  expect_equal(fit$summary_stats$Strata, 4)
  expect_equal(fit$summary_stats$PSUs, 16)
  expect_equal(fit$summary_stats$DF, 12)
  expect_equal(
    unname(fit$raw$design$weights),
    rep(1, nrow(data)),
    tolerance = 0
  )
  expect_match(
    fit$raw$known_limitations$survey_design,
    "ultimate-cluster"
  )

  expect_error(
    wvw_decomp(data, "y", c("x1", "x2"), "ses",
               use_svy = NA, quiet = TRUE),
    "use_svy"
  )
  expect_error(
    wvw_decomp(data, "y", c("x1", "x2"), "ses",
               lonely_psu = "invented", quiet = TRUE),
    "lonely_psu"
  )
})

test_that("phase 10 Taylor covariance matches the frozen Stata design", {
  fit <- .phase10_wvw_fit(.phase10_wvw_data(8L))
  observed <- fit$raw$linearization$vcov
  expected <- matrix(
    c(
      0, 0, 0, 0, 0,
      0, 2.1540629296100001e-04, 4.2310938787046146e-06,
      2.1387285564016380e-04, -5.7645311995410015e-06,
      0, 4.2310938787046146e-06, 9.5803531651422774e-07,
      4.8601238641775072e-06, -3.2900533104133024e-07,
      0, 2.1387285564016380e-04, 4.8601238641775072e-06,
      2.2269962699656708e-04, 3.9666474922255693e-06,
      0, -5.7645311995410015e-06, -3.2900533104133024e-07,
      3.9666474922255693e-06, 1.0060184022807880e-05
    ),
    nrow = 5L,
    byrow = TRUE,
    dimnames = dimnames(observed)
  )
  benchmark <-
    fit$raw$linearization$conindex_survey_benchmark$total

  expect_equal(observed, expected, tolerance = 1e-14)
  expect_equal(
    benchmark$estimate,
    -0.0050375773304648375,
    tolerance = 1e-8
  )
  expect_equal(
    benchmark$std_error,
    0.01541175487672073,
    tolerance = 1e-8
  )
  expect_equal(benchmark$degrees_of_freedom, 12)
  expect_match(
    fit$raw$linearization$finite_sample_correction,
    "survey Taylor"
  )
})

test_that("phase 10 JKn uses stratum centers and design scale", {
  fit <- .phase10_wvw_fit(method = "jackknife")
  replication <- fit$raw$replication
  factors <- replication$replicate_factors
  strata <- unname(fit$raw$design$strata)

  center_groups <- vapply(
    seq_len(ncol(factors)),
    function(column) {
      unique(strata[abs(factors[, column] - 1) > 1e-12])
    },
    character(1)
  )
  expected_centers <- do.call(
    rbind,
    lapply(center_groups, function(group) {
      colMeans(
        replication$estimates[
          center_groups == group,
          ,
          drop = FALSE
        ]
      )
    })
  )

  expect_false(replication$ordinary_resampling)
  expect_false(replication$ordinary_bootstrap)
  expect_false(replication$ordinary_jackknife)
  expect_identical(
    replication$engine,
    "survey::as.svrepdesign"
  )
  expect_identical(
    replication$fit_weight_type,
    "survey replicate weight"
  )
  expect_identical(replication$replicates, replication$estimates)
  expect_identical(replication$survey_replicate_type, "JKn")
  expect_identical(
    replication$center_method,
    "stratum-specific replicate mean"
  )
  expect_equal(replication$center, expected_centers, tolerance = 0)
  expect_equal(replication$generator_scale, 1, tolerance = 0)
  expect_equal(replication$effective_scale, 1, tolerance = 0)
  expect_equal(replication$requested_replicates, 16)
  expect_equal(replication$failed_replicates, 0)
  expect_equal(fit$summary_stats$DF, 12)
  expect_equal(
    replication$vcov,
    .phase10_reconstruct_survey_vcov(replication),
    tolerance = 1e-14
  )

  for (column in seq_len(ncol(factors))) {
    changed <- abs(factors[, column] - 1) > 1e-12
    expect_length(unique(strata[changed]), 1)
    expect_equal(sum(factors[, column] == 0), 4)
    expect_equal(
      unique(factors[changed & factors[, column] > 0, column]),
      4 / 3,
      tolerance = 1e-14
    )
  }
})

test_that("phase 10 bootstrap is Rao-Wu rescaled", {
  data <- .phase10_wvw_data(4L)
  fit <- .phase10_wvw_fit(
    data = data,
    method = "bootstrap",
    boot_reps = 20L
  )
  replication <- fit$raw$replication
  factors <- replication$replicate_factors
  strata <- unname(fit$raw$design$strata)
  psu <- unname(fit$raw$design$psu)

  expect_false(replication$ordinary_resampling)
  expect_false(replication$ordinary_bootstrap)
  expect_false(replication$ordinary_jackknife)
  expect_identical(
    replication$engine,
    "survey::as.svrepdesign"
  )
  expect_identical(
    replication$fit_weight_type,
    "survey replicate weight"
  )
  expect_identical(replication$replicates, replication$estimates)
  expect_identical(
    replication$survey_replicate_type,
    "Rao-Wu rescaled bootstrap"
  )
  expect_identical(replication$center_method, "replicate mean")
  expect_equal(replication$generator_scale, 1 / 19,
               tolerance = 1e-14)
  expect_equal(replication$effective_scale, 1 / 20,
               tolerance = 1e-14)
  expect_equal(replication$requested_replicates, 20)
  expect_equal(replication$failed_replicates, 0)
  expect_equal(fit$summary_stats$DF, 12)
  expect_equal(
    replication$vcov,
    .phase10_reconstruct_survey_vcov(replication),
    tolerance = 1e-14
  )

  for (column in seq_len(ncol(factors))) {
    for (stratum in unique(strata)) {
      rows <- strata == stratum
      psu_factors <- vapply(
        split(factors[rows, column], psu[rows]),
        function(values) unique(values)[1],
        numeric(1)
      )
      expect_equal(sum(psu_factors), 4, tolerance = 1e-14)
      expect_equal(
        psu_factors / (4 / 3),
        round(psu_factors / (4 / 3)),
        tolerance = 1e-14
      )
    }
  }
})

test_that("phase 10 nests reused PSU labels within strata", {
  data <- .phase10_wvw_data(4L)
  for (method in c("linearized", "jackknife", "bootstrap")) {
    unique_fit <- .phase10_wvw_fit(
      data,
      method,
      psu = "psu",
      boot_reps = 10L
    )
    reused_fit <- .phase10_wvw_fit(
      data,
      method,
      psu = "psu_reused",
      boot_reps = 10L
    )

    expect_equal(
      unique_fit$raw$estimates$overall,
      reused_fit$raw$estimates$overall,
      tolerance = 0,
      info = method
    )
    expect_equal(
      unique_fit$raw$vcov$overall,
      reused_fit$raw$vcov$overall,
      tolerance = 0,
      info = method
    )
    expect_equal(reused_fit$summary_stats$PSUs, 16)
    expect_equal(reused_fit$summary_stats$DF, 12)
  }
})

test_that("phase 10 lonely-PSU rules match Stata mappings", {
  data <- .phase10_wvw_data(4L)
  data$psu[data$strata == 1] <- 1

  fits <- lapply(
    c("adjust", "certainty", "remove", "average"),
    function(rule) {
      .phase10_wvw_fit(
        data = data,
        lonely_psu = rule
      )
    }
  )
  names(fits) <- c("adjust", "certainty", "remove", "average")
  covariances <- lapply(
    fits,
    function(fit) fit$raw$linearization$vcov
  )

  expect_equal(
    covariances$certainty,
    covariances$remove,
    tolerance = 0
  )
  expect_false(isTRUE(all.equal(
    covariances$adjust,
    covariances$certainty,
    tolerance = 1e-15
  )))
  expect_false(isTRUE(all.equal(
    covariances$average,
    covariances$certainty,
    tolerance = 1e-15
  )))
  expect_true(all(vapply(
    covariances,
    function(covariance) all(is.finite(covariance)),
    logical(1)
  )))
  expect_true(all(vapply(
    fits,
    function(fit) fit$summary_stats$DF == 9,
    logical(1)
  )))
  expect_error(
    .phase10_wvw_fit(data = data, lonely_psu = "fail"),
    "one PSU"
  )
})

test_that("phase 10 rejects singleton strata for survey replication", {
  data <- .phase10_wvw_data(4L)
  data$psu[data$strata == 1] <- 1

  for (method in c("jackknife", "bootstrap")) {
    expect_error(
      .phase10_wvw_fit(
        data = data,
        method = method,
        boot_reps = 8L
      ),
      "at least two PSUs",
      info = method
    )
  }
})

test_that("phase 10 bootstrap supports singleton certainty through svrep", {
  data <- .phase10_wvw_data(4L)
  data$psu[data$strata == 1] <- 1
  bootstrap_reps <- 12L

  fit <- .phase10_wvw_fit(
    data = data,
    method = "bootstrap",
    boot_reps = bootstrap_reps,
    bootstrap_singleton = "certainty"
  )
  replication <- fit$raw$replication

  expect_identical(
    replication$engine,
    "svrep::as_bootstrap_design"
  )
  expect_identical(
    replication$survey_replicate_type,
    "Rao-Wu-Yue-Beaumont bootstrap with singleton certainty"
  )
  expect_true(replication$bootstrap_singleton$applied)
  expect_equal(replication$bootstrap_singleton$n_singleton_strata, 1)
  expect_equal(replication$generator_scale, 1 / bootstrap_reps)
  expect_equal(replication$effective_scale, 1 / bootstrap_reps)

  factors <- replication$replicate_factors
  strata <- unname(fit$raw$design$strata)
  psu <- unname(fit$raw$design$psu)
  expect_equal(
    unname(factors[strata == 1, , drop = FALSE]),
    matrix(
      1,
      nrow = sum(strata == 1),
      ncol = bootstrap_reps
    ),
    tolerance = 1e-12
  )
  for (stratum in setdiff(unique(strata), 1)) {
    rows <- strata == stratum
    psu_factors <- vapply(
      split(factors[rows, 1], psu[rows]),
      function(values) unique(values)[1],
      numeric(1)
    )
    expect_equal(sum(psu_factors), 4, tolerance = 1e-12)
    expect_equal(
      psu_factors / (4 / 3),
      round(psu_factors / (4 / 3)),
      tolerance = 1e-12
    )
  }
  expect_match(
    fit$diagnostics$bootstrap_singleton_control,
    "lonely_psu does not control Rao-Wu bootstrap"
  )
  expect_match(
    fit$diagnostics$bootstrap_singleton_assumption,
    "zero first-stage variance"
  )
})

test_that("phase 10 validates bootstrap_singleton independently", {
  expect_error(
    .phase10_wvw_fit(bootstrap_singleton = "average"),
    "bootstrap_singleton must be"
  )
})

test_that("phase 10 distinguishes simple-weight inference by VCE", {
  data <- .phase10_wvw_data(2L)
  arguments <- list(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    weight_var = "weight",
    quiet = TRUE
  )
  linearized <- do.call(
    wvw_decomp,
    c(arguments, list(vce_method = "linearized"))
  )
  jackknife <- do.call(
    wvw_decomp,
    c(arguments, list(vce_method = "jackknife"))
  )
  bootstrap <- do.call(
    wvw_decomp,
    c(
      arguments,
      list(
        vce_method = "bootstrap",
        boot_reps = 8L,
        seed = 73L
      )
    )
  )

  expect_equal(linearized$summary_stats$DF, Inf)
  expect_equal(jackknife$summary_stats$DF, nrow(data) - 1)
  expect_equal(bootstrap$summary_stats$DF, nrow(data) - 1)
  expect_identical(
    jackknife$raw$replication$survey_replicate_type,
    "JK1"
  )
  expect_identical(
    bootstrap$raw$replication$survey_replicate_type,
    "Rao-Wu rescaled bootstrap"
  )
})
