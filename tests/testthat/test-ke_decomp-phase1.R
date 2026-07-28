phase1_ke_data <- function() {
  n <- 24
  data.frame(
    health = pmin(
      pmax(
        0.18 +
          0.021 * seq_len(n) +
          0.06 * rep(c(0, 1), 12) -
          0.025 * rep(c(0, 0, 1), 8),
        0
      ),
      1
    ),
    ses = rep(seq_len(12), each = 2),
    income = 100 + 8 * seq_len(n),
    x1 = seq_len(n) + rep(c(-0.2, 0.3), 12),
    category = factor(
      rep(c("base", "middle", "high"), 8),
      levels = c("base", "middle", "high")
    ),
    weight = rep(c(1, 2, 3, 1), 6),
    strata = rep(c("s1", "s2", "s3"), each = 8),
    psu = rep(rep(c("p1", "p2", "p3", "p4"), each = 2), 3),
    stringsAsFactors = FALSE
  )
}

test_that("KE phase 1 returns unrounded public and raw audit results", {
  data <- phase1_ke_data()
  data <- rbind(
    data,
    data.frame(
      health = NA,
      ses = 13,
      income = 310,
      x1 = 25,
      category = factor(
        "base",
        levels = levels(data$category)
      ),
      weight = 2,
      strata = "s3",
      psu = "p4"
    )
  )

  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "category"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expect_true(is.list(fit$raw))
  expect_identical(fit$raw$schema_version, "1.0")
  expect_identical(fit$raw$settings$numeric_rounding, "none")
  expect_identical(fit$raw$settings$index_type, "rank")

  expect_equal(
    fit$results_overall$Estimate,
    unname(fit$raw$estimates$overall),
    tolerance = 0
  )
  expect_equal(
    fit$results_detailed$Estimate,
    unname(fit$raw$estimates$detailed),
    tolerance = 0
  )
  expect_true(any(
    abs(
      fit$results_detailed$Estimate -
        round(fit$results_detailed$Estimate, 6)
    ) > 1e-10
  ))

  expect_equal(fit$raw$sample$input_n, 25)
  expect_equal(fit$raw$sample$analytic_n, 24)
  expect_equal(
    fit$raw$sample$excluded_source_rows$missing_required_values,
    25
  )
  expect_equal(nrow(fit$raw$target), 24)
  expect_equal(
    fit$raw$target$source_row,
    fit$raw$sample$analytic_source_rows
  )
  expect_equal(
    rownames(fit$raw$model$matrix),
    as.character(fit$raw$sample$analytic_source_rows)
  )
  expect_equal(
    names(fit$raw$model$coefficients),
    colnames(fit$raw$model$matrix)
  )
  expect_equal(
    dim(fit$raw$vcov$model),
    rep(ncol(fit$raw$model$matrix), 2)
  )
  expect_false(fit$raw$vcov$overall_complete)
  expect_true(is.na(fit$raw$vcov$overall[1, 1]))
  expect_null(fit$raw$replication)

  expect_equal(
    sum(
      fit$raw$target$normalized_weight *
        fit$raw$target$target
    ),
    fit$results_overall$Estimate,
    tolerance = 1e-15
  )
  expect_equal(
    fit$raw$index$regression_at_covariate_means,
    fit$results_overall$Estimate,
    tolerance = 1e-12
  )
  expect_equal(
    fit$raw$index$weighted_residual_mean,
    0,
    tolerance = 1e-12
  )
})

test_that("KE phase 1 exposes generalized weighted tied ranks", {
  data <- phase1_ke_data()
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "category"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    use_svy = TRUE,
    weight_var = "weight",
    vce_method = "linearized",
    quiet = TRUE
  ))

  ordered_rows <- order(data$ses)
  ordered_data <- data[ordered_rows, , drop = FALSE]
  normalized_weight <- ordered_data$weight / sum(ordered_data$weight)
  group_mass <- ave(
    normalized_weight,
    ordered_data$ses,
    FUN = sum
  )
  first_in_group <- !duplicated(ordered_data$ses)
  cumulative_mass <- cumsum(
    ifelse(first_in_group, group_mass, 0)
  )
  expected_rank <- cumulative_mass - 0.5 * group_mass

  expect_equal(
    fit$raw$target$source_row,
    ordered_rows,
    tolerance = 0
  )
  expect_equal(
    fit$raw$target$fractional_rank,
    expected_rank,
    tolerance = 1e-15
  )
  expect_equal(
    fit$raw$target$fractional_rank[seq(1, 24, by = 2)],
    fit$raw$target$fractional_rank[seq(2, 24, by = 2)],
    tolerance = 0
  )
})

test_that("KE phase 1 exposes the level-dependent target", {
  data <- phase1_ke_data()
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "category"),
    ses_var = "income",
    index_type = "level",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expected_level <- data$income / mean(data$income)
  expected_target <- expected_level * data$health - mean(data$health)

  expect_true(all(is.na(fit$raw$target$fractional_rank)))
  expect_equal(
    fit$raw$target$relative_level,
    expected_level,
    tolerance = 1e-15
  )
  expect_equal(
    fit$raw$target$target_base,
    expected_target,
    tolerance = 1e-15
  )
  expect_equal(
    fit$raw$target$target,
    expected_target,
    tolerance = 1e-15
  )
  expect_equal(
    fit$raw$index$socioeconomic_mean,
    mean(data$income),
    tolerance = 1e-15
  )
})

test_that("KE phase 1 exposes current bootstrap replication internals", {
  data <- phase1_ke_data()
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "category"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "bootstrap",
    boot_reps = 12,
    seed = 20260728,
    quiet = TRUE
  ))

  replication <- fit$raw$replication
  expect_true(is.list(replication))
  expect_identical(replication$method, "bootstrap")
  expect_equal(replication$requested_replicates, 12)
  expect_equal(nrow(replication$estimates), 12)
  expect_equal(replication$valid_replicates, 12)
  expect_equal(replication$failed_replicates, 0)
  expect_equal(
    ncol(replication$estimates),
    ncol(fit$raw$model$matrix)
  )
  expect_identical(
    replication$center_method,
    "replicate mean"
  )
  expect_identical(replication$engine, "boot::boot")
  expect_true(replication$ordinary_bootstrap)
  expect_false(replication$pending_methodological_review)
  expect_null(fit$raw$known_limitations$replication)
  expect_true(is.na(fit$results_overall$Std_Error))
  expect_null(replication$index_variance)
  expect_equal(
    fit$raw$vcov$model,
    replication$covariance,
    tolerance = 0
  )
})
