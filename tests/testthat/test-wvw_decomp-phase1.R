phase1_wvw_data <- function() {
  n <- 18
  data.frame(
    y = 1.7 +
      0.31 * c(2, 5, 1, 7, 4, 8, 3, 9, 6, 11, 10, 14, 12, 16, 13, 18, 15, 17) -
      0.22 * rep(c(0, 1), 9) +
      c(0.11, -0.08, 0.04, -0.13, 0.16, -0.03, 0.07, -0.09, 0.05,
        -0.06, 0.12, -0.04, 0.08, -0.11, 0.03, -0.07, 0.09, -0.02),
    x1 = c(2, 5, 1, 7, 4, 8, 3, 9, 6, 11, 10, 14, 12, 16, 13, 18, 15, 17),
    x2 = rep(c(0, 1), 9),
    ses = rep(seq_len(9), each = 2),
    w = c(1, 2, 3, 1, 2, 4, 1, 3, 2, 1, 4, 2, 3, 1, 2, 3, 1, 4),
    stringsAsFactors = FALSE
  )
}

test_that("phase 1 returns unrounded public and raw audit results", {
  data <- phase1_wvw_data()
  data <- rbind(
    data,
    data.frame(y = NA, x1 = 20, x2 = 1, ses = 10, w = 2)
  )

  fit <- suppressWarnings(wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    weight_var = "w",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expect_true(is.list(fit$raw))
  expect_identical(fit$raw$schema_version, "1.0")
  expect_identical(fit$raw$settings$numeric_rounding, "none")
  expect_match(fit$raw$settings$rank_method, "generalized weighted")

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
      fit$results_overall$Estimate -
        round(fit$results_overall$Estimate, 6)
    ) > 1e-10
  ))

  expect_equal(fit$raw$sample$input_n, 19)
  expect_equal(fit$raw$sample$analytic_n, 18)
  expect_equal(
    fit$raw$sample$excluded_source_rows$missing_required_values,
    19
  )
  expect_equal(nrow(fit$raw$model$matrix), 18)
  expect_equal(ncol(fit$raw$model$matrix), 3)
  expect_equal(
    rownames(fit$raw$model$matrix),
    as.character(fit$raw$sample$analytic_source_rows)
  )
  expect_equal(
    names(fit$raw$model$coefficients),
    colnames(fit$raw$model$matrix)
  )
  expect_equal(
    dim(fit$raw$vcov$detailed),
    c(ncol(fit$raw$model$matrix), ncol(fit$raw$model$matrix))
  )
  expect_true(fit$raw$vcov$detailed_complete)
  expect_true(fit$raw$vcov$overall_complete)
  expect_true(is.matrix(
    fit$raw$linearization$influence_function
  ))
  expect_null(fit$raw$replication)
})

test_that("phase 1 exposes the generalized tied rank and source rows", {
  data <- phase1_wvw_data()
  fit <- suppressWarnings(wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    weight_var = "w",
    vce_method = "linearized",
    quiet = TRUE
  ))

  ordered_rows <- order(data$ses)
  ordered_data <- data[ordered_rows, , drop = FALSE]
  tied_weight <- ave(
    ordered_data$w,
    ordered_data$ses,
    FUN = sum
  )
  first_in_tie <- !duplicated(ordered_data$ses)
  cumulative_tied_weight <- cumsum(
    ifelse(first_in_tie, tied_weight, 0)
  )
  expected_rank <- (
    cumulative_tied_weight - 0.5 * tied_weight
  ) / sum(ordered_data$w)

  expect_equal(
    fit$raw$rank$source_row,
    ordered_rows,
    tolerance = 0
  )
  expect_equal(
    fit$raw$rank$fractional_rank,
    expected_rank,
    tolerance = 1e-15
  )
  expect_equal(
    fit$raw$rank$fractional_rank[seq(1, 18, by = 2)],
    fit$raw$rank$fractional_rank[seq(2, 18, by = 2)],
    tolerance = 0
  )
})

test_that("phase 1 exposes validated ordinary jackknife internals", {
  data <- phase1_wvw_data()
  fit <- suppressWarnings(wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    vce_method = "jackknife",
    quiet = TRUE
  ))

  replication <- fit$raw$replication
  expect_true(is.list(replication))
  expect_identical(replication$method, "jackknife")
  expect_identical(replication$mse, FALSE)
  expect_identical(
    replication$current_center_method,
    "replicate mean"
  )
  expect_false(replication$pending_methodological_review)
  expect_true(replication$ordinary_resampling)
  expect_equal(replication$requested_replicates, nrow(data))
  expect_equal(nrow(replication$estimates), nrow(data))
  expect_equal(nrow(replication$replicate_factors), nrow(data))
  expect_equal(ncol(replication$replicate_factors), nrow(data))
  expect_equal(
    fit$results_overall$Estimate,
    unname(fit$raw$estimates$overall),
    tolerance = 0
  )
  expect_true(fit$raw$vcov$overall_complete)
})

test_that("phase 1 exposes validated ordinary bootstrap internals", {
  data <- phase1_wvw_data()
  fit <- suppressWarnings(wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    vce_method = "bootstrap",
    boot_reps = 12,
    seed = 731,
    quiet = TRUE
  ))

  replication <- fit$raw$replication
  expect_identical(replication$method, "bootstrap")
  expect_equal(replication$requested_replicates, 12)
  expect_equal(nrow(replication$estimates), 12)
  expect_equal(ncol(replication$replicate_factors), 12)
  expect_length(replication$rscales, 12)
  expect_true(is.matrix(replication$vcov))
  expect_identical(replication$mse, FALSE)
  expect_identical(replication$center_method, "replicate mean")
  expect_false(replication$pending_methodological_review)
  expect_true(replication$ordinary_resampling)
})
