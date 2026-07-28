.phase9_ke_data <- function() {
  n_strata <- 6L
  psu_per_stratum <- 4L
  observations_per_psu <- 5L
  n <- n_strata * psu_per_stratum * observations_per_psu
  source_id <- seq_len(n)
  strata <- rep(
    seq_len(n_strata),
    each = psu_per_stratum * observations_per_psu
  )
  psu <- rep(
    rep(
      seq_len(psu_per_stratum),
      each = observations_per_psu
    ),
    n_strata
  )
  ses <- rep(seq_len(30L), each = 4L)
  x1 <- sin(source_id / 7) + source_id / 90
  x2 <- as.numeric((source_id * 5L) %% 11L < 4L)
  weight <- (
    0.8 + 0.12 * strata + 0.09 * psu +
      0.04 * (source_id %% observations_per_psu)
  )
  data <- data.frame(
    source_id = source_id,
    health = 0.26 +
      0.052 * x1 -
      0.033 * x2 +
      0.0028 * ses +
      0.012 * cos(source_id / 9) +
      0.006 * strata,
    ses = ses,
    income = 220 + 6.5 * source_id +
      18 * sin(source_id / 8) + 4 * strata,
    x1 = x1,
    x2 = x2,
    weight = weight,
    weight_scaled = 100 * weight,
    strata = strata,
    psu = psu,
    psu_global = (strata - 1L) * psu_per_stratum + psu
  )
  input_order <- c(
    seq(3, n, by = 3),
    seq(1, n, by = 3),
    seq(2, n, by = 3)
  )
  data[input_order, , drop = FALSE]
}

.phase9_ke_fit <- function(
    data, method, strata = "strata", psu = "psu",
    weight = "weight", index_type = "rank",
    correction = "generalized", boot_reps = 20L,
    seed = 990901L, ...) {
  suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = if (index_type == "rank") "ses" else "income",
    index_type = index_type,
    correction = correction,
    dep_min = if (correction %in% c("erreygers", "wagstaff")) {
      0
    } else {
      NULL
    },
    dep_max = if (correction %in% c("erreygers", "wagstaff")) {
      1
    } else {
      NULL
    },
    weight_var = weight,
    strata_var = strata,
    psu_var = psu,
    vce_method = method,
    boot_reps = boot_reps,
    seed = seed,
    quiet = TRUE,
    ...
  ))
}

.phase9_ke_replicate_vcov <- function(replication) {
  replicates <- replication$valid_estimates
  groups <- replication$center_groups
  if (is.null(groups)) {
    differences <- sweep(
      replicates,
      2,
      colMeans(replicates),
      "-"
    )
  } else {
    centers <- lapply(
      unique(groups),
      function(group) {
        colMeans(
          replicates[groups == group, , drop = FALSE]
        )
      }
    )
    names(centers) <- unique(groups)
    center_matrix <- do.call(
      rbind,
      lapply(groups, function(group) centers[[group]])
    )
    differences <- replicates - center_matrix
  }
  replication$effective_scale *
    crossprod(
      differences * sqrt(replication$used_rscales)
    )
}

test_that("KE phase 9 selects JK1 for unstratified survey designs", {
  data <- .phase9_ke_data()
  observation_jk1 <- .phase9_ke_fit(
    data,
    method = "jackknife",
    strata = NULL,
    psu = NULL
  )
  replication <- observation_jk1$raw$replication

  expect_identical(replication$survey_replicate_type, "JK1")
  expect_identical(replication$center_method, "replicate mean")
  expect_equal(replication$requested_replicates, nrow(data))
  expect_equal(replication$failed_replicates, 0)
  expect_equal(
    replication$generator_scale,
    (nrow(data) - 1) / nrow(data)
  )
  expect_equal(observation_jk1$summary_stats$DF, nrow(data) - 1)
  expect_equal(
    replication$vcov,
    .phase9_ke_replicate_vcov(replication),
    tolerance = 2e-14
  )

  cluster_jk1 <- .phase9_ke_fit(
    data,
    method = "jackknife",
    strata = NULL,
    psu = "psu_global",
    index_type = "level",
    correction = "standard"
  )
  expect_identical(
    cluster_jk1$raw$replication$survey_replicate_type,
    "JK1"
  )
  expect_equal(
    cluster_jk1$raw$replication$requested_replicates,
    24
  )
  expect_equal(
    cluster_jk1$raw$replication$effective_scale,
    23 / 24
  )
  expect_equal(cluster_jk1$summary_stats$DF, 23)
})

test_that("KE phase 9 uses stratum-centered JKn", {
  data <- .phase9_ke_data()
  fit <- .phase9_ke_fit(data, method = "jackknife")
  replication <- fit$raw$replication

  expect_identical(replication$survey_replicate_type, "JKn")
  expect_identical(
    replication$center_method,
    "stratum-specific replicate mean"
  )
  expect_equal(replication$requested_replicates, 24)
  expect_equal(replication$failed_replicates, 0)
  expect_equal(length(unique(replication$center_groups)), 6)
  expect_equal(
    as.integer(table(replication$center_groups)),
    rep(4L, 6)
  )
  expect_equal(replication$generator_scale, 1)
  expect_equal(replication$effective_scale, 1)
  expect_equal(replication$used_rscales, rep(3 / 4, 24))
  expect_equal(fit$summary_stats$DF, 18)
  expect_equal(
    replication$vcov,
    .phase9_ke_replicate_vcov(replication),
    tolerance = 2e-14
  )
  expect_false(replication$pending_methodological_review)
  expect_true(is.na(fit$results_overall$Std_Error))
})

test_that("KE phase 9 uses Rao-Wu survey bootstrap and Stata scale", {
  data <- .phase9_ke_data()
  fit <- .phase9_ke_fit(
    data,
    method = "bootstrap",
    boot_reps = 25L
  )
  replication <- fit$raw$replication

  expect_identical(
    replication$survey_replicate_type,
    "Rao-Wu rescaled bootstrap"
  )
  expect_identical(replication$center_method, "replicate mean")
  expect_equal(replication$requested_replicates, 25)
  expect_equal(replication$failed_replicates, 0)
  expect_equal(replication$generator_scale, 1 / 24)
  expect_equal(replication$effective_scale, 1 / 25)
  expect_equal(replication$used_rscales, rep(1, 25))
  expect_equal(fit$summary_stats$DF, 18)
  expect_equal(
    replication$vcov,
    .phase9_ke_replicate_vcov(replication),
    tolerance = 2e-14
  )

  factors <- replication$replicate_factors
  source_rows <- fit$raw$sample$analytic_source_rows
  analytic_data <- data[source_rows, , drop = FALSE]
  nested_psu <- interaction(
    analytic_data$strata,
    analytic_data$psu,
    drop = TRUE
  )
  psu_rows <- !duplicated(nested_psu)
  psu_strata <- analytic_data$strata[psu_rows]
  for (replicate_index in seq_len(ncol(factors))) {
    psu_factor <- factors[psu_rows, replicate_index]
    for (stratum in unique(psu_strata)) {
      stratum_factor <- psu_factor[psu_strata == stratum]
      counts <- stratum_factor * 3 / 4
      expect_equal(counts, round(counts), tolerance = 2e-14)
      expect_equal(sum(stratum_factor), 4, tolerance = 2e-14)
    }
  }
})

test_that("KE phase 9 survey replicas ignore global weight scaling", {
  data <- .phase9_ke_data()
  arguments <- list(
    data = data,
    method = "bootstrap",
    boot_reps = 20L
  )
  original <- do.call(
    .phase9_ke_fit,
    c(arguments, list(weight = "weight"))
  )
  scaled <- do.call(
    .phase9_ke_fit,
    c(arguments, list(weight = "weight_scaled"))
  )

  expect_equal(
    original$raw$replication$replicate_factors,
    scaled$raw$replication$replicate_factors,
    tolerance = 0
  )
  expect_equal(
    original$raw$replication$estimates,
    scaled$raw$replication$estimates,
    tolerance = 2e-14
  )
  expect_equal(
    original$raw$model$vcov,
    scaled$raw$model$vcov,
    tolerance = 2e-14
  )
})

test_that("KE phase 9 rejects singleton strata for replication", {
  data <- .phase9_ke_data()
  data$psu[data$strata == 1] <- 1

  expect_error(
    .phase9_ke_fit(data, method = "jackknife"),
    "at least two PSUs in every stratum"
  )
  expect_error(
    .phase9_ke_fit(
      data,
      method = "bootstrap",
      boot_reps = 12L
    ),
    "at least two PSUs in every stratum"
  )
  expect_error(
    ke_decomp(
      data = data,
      dep_var = "health",
      indep_vars = c("x1", "x2"),
      ses_var = "ses",
      lonely_psu = "unknown",
      quiet = TRUE
    ),
    "lonely_psu must be"
  )
})

test_that("KE phase 9 does not silently discard survey replicates", {
  n <- 18L
  data <- data.frame(
    health = rep(c(1, 0, 0), each = 6L),
    ses = rep(seq_len(9L), each = 2L),
    x1 = sin(seq_len(n) / 3),
    x2 = as.numeric(seq_len(n) %% 4L == 0L),
    weight = 0.9 + seq_len(n) / 40,
    strata = 1,
    psu = rep(seq_len(3L), each = 6L)
  )
  arguments <- list(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "wagstaff",
    dep_min = 0,
    dep_max = 1,
    weight_var = "weight",
    strata_var = "strata",
    psu_var = "psu",
    quiet = TRUE
  )

  expect_error(
    suppressWarnings(do.call(
      ke_decomp,
      c(arguments, list(vce_method = "jackknife"))
    )),
    "survey replicate\\(s\\) failed"
  )
  jackknife <- suppressWarnings(do.call(
    ke_decomp,
    c(
      arguments,
      list(vce_method = "jackknife", relax = TRUE)
    )
  ))
  expect_equal(jackknife$raw$replication$failed_replicates, 1)
  expect_equal(jackknife$raw$replication$valid_replicates, 2)
  expect_equal(
    jackknife$raw$replication$failed_replicate_adjustment,
    3 / 2
  )

  expect_error(
    suppressWarnings(do.call(
      ke_decomp,
      c(
        arguments,
        list(
          vce_method = "bootstrap",
          boot_reps = 20L,
          seed = 990903L
        )
      )
    )),
    "survey replicate\\(s\\) failed"
  )
  bootstrap <- suppressWarnings(do.call(
    ke_decomp,
    c(
      arguments,
      list(
        vce_method = "bootstrap",
        boot_reps = 20L,
        seed = 990903L,
        relax = TRUE
      )
    )
  ))
  expect_gt(bootstrap$raw$replication$failed_replicates, 0)
  expect_gt(bootstrap$raw$replication$valid_replicates, 1)
  expect_equal(
    bootstrap$raw$replication$failed_replicate_adjustment,
    bootstrap$raw$replication$requested_replicates /
      bootstrap$raw$replication$valid_replicates
  )
})
