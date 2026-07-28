.phase3_ke_data <- function() {
  n <- 72L
  source_id <- seq_len(n)
  tie_sizes <- c(7L, 9L, 11L, 13L, 15L, 17L)
  tie_group <- rep(seq_along(tie_sizes), times = tie_sizes)
  ses_unique <- 35 + ((31 * source_id) %% 97) + source_id / 1000
  x1 <- ((11 * source_id) %% 27) - 13 + source_id / 19
  x2 <- cos(source_id / 5) + 0.17 * (source_id %% 4)
  ses_scaled <- (ses_unique - min(ses_unique)) /
    (max(ses_unique) - min(ses_unique))
  health <- 0.32 +
    0.006 * x1 -
    0.038 * x2 +
    0.15 * ses_scaled +
    0.021 * sin(source_id / 4) -
    0.009 * cos(source_id / 9)
  input_order <- c(
    seq(3, n, by = 3),
    seq(1, n, by = 3),
    seq(2, n, by = 3)
  )

  data.frame(
    source_id = source_id,
    health = health,
    ses_unique = ses_unique,
    ses_tied = c(1, 2, 4, 7, 11, 18)[tie_group],
    ses_two_ties = ifelse(source_id <= 29L, 1, 2),
    ses_all_tied = 1,
    income = 380 + 10 * source_id + 28 * sin(source_id / 6),
    x1 = x1,
    x2 = x2,
    weight = ((5 * source_id) %% 7) + 1,
    permutation_key = ((29 * source_id) %% 73) +
      source_id / 1000,
    stringsAsFactors = FALSE
  )[input_order, , drop = FALSE]
}

.phase3_ke_grouped_rank <- function(ses, weight) {
  ord <- order(ses, seq_along(ses))
  ses_ordered <- ses[ord]
  p <- weight[ord] / sum(weight)
  group_mass <- ave(p, ses_ordered, FUN = sum)
  first <- !duplicated(ses_ordered)
  cumulative <- cumsum(ifelse(first, group_mass, 0))
  ordered_rank <- cumulative - 0.5 * group_mass
  result <- numeric(length(ses))
  result[ord] <- ordered_rank
  result
}

.phase3_ke_fit <- function(
    data, ses_var = "ses_tied", index_type = "rank") {
  suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = ses_var,
    index_type = index_type,
    correction = "generalized",
    use_svy = TRUE,
    weight_var = "weight",
    vce_method = "linearized",
    quiet = TRUE
  ))
}

test_that("KE phase 3 weighted tied ranks reproduce Appendix A and software", {
  data <- .phase3_ke_data()
  fit <- .phase3_ke_fit(data)
  raw <- fit$raw
  expected_rank <- .phase3_ke_grouped_rank(
    data$ses_tied,
    data$weight
  )
  expected_by_source <- expected_rank[
    match(
      data$source_id[raw$target$source_row],
      data$source_id
    )
  ]
  rineq_rank <- rineq::rank_gwt(data$ses_tied, data$weight)
  rineq_by_source <- rineq_rank[
    match(
      data$source_id[raw$target$source_row],
      data$source_id
    )
  ]

  expect_equal(
    raw$target$fractional_rank,
    expected_by_source,
    tolerance = 2e-14
  )
  expect_equal(
    raw$target$fractional_rank,
    rineq_by_source,
    tolerance = 2e-14
  )
  expect_equal(
    sum(raw$target$normalized_weight *
      raw$target$fractional_rank),
    0.5,
    tolerance = 2e-14
  )

  p <- data$weight / sum(data$weight)
  appendix_weight <- nrow(data) * p *
    (2 * expected_rank - 1)
  appendix_index <- sum(appendix_weight * data$health) /
    nrow(data)
  rineq_index <- rineq::ci(
    ineqvar = data$ses_tied,
    outcome = data$health,
    weights = data$weight,
    type = "CIg",
    method = "direct",
    df_correction = FALSE,
    rank_function = rineq::rank_gwt
  )$concentration_index

  expect_lt(
    abs(fit$results_overall$Estimate - appendix_index),
    2e-14
  )
  expect_lt(
    abs(fit$results_overall$Estimate - rineq_index),
    2e-14
  )
  expect_lt(
    abs(
      fit$results_overall$Estimate -
        (-0.0027796326505051161)
    ),
    2e-14
  )
  expect_equal(
    unname(raw$model$coefficients),
    c(
      -0.029423929122822073,
      0.010942885703556556,
      0.025487323928020923
    ),
    tolerance = 2e-14
  )
  expect_lt(
    abs(
      raw$index$regression_at_covariate_means -
        fit$results_overall$Estimate
    ),
    2e-14
  )
  expect_match(
    raw$design$point_weight_semantics,
    "sampling/design weight"
  )
  expect_match(
    raw$design$point_estimation_equation,
    "frequency expansion"
  )
})

test_that("KE phase 3 ranks ignore weight scale and within-tie row order", {
  data <- .phase3_ke_data()
  base <- .phase3_ke_fit(data)

  scaled_data <- data
  scaled_data$weight <- 10 * scaled_data$weight
  scaled <- .phase3_ke_fit(scaled_data)

  permuted_data <- data[order(data$permutation_key), , drop = FALSE]
  permuted <- .phase3_ke_fit(permuted_data)

  expect_equal(
    scaled$results_overall$Estimate,
    base$results_overall$Estimate,
    tolerance = 2e-14
  )
  expect_equal(
    scaled$raw$model$coefficients,
    base$raw$model$coefficients,
    tolerance = 2e-14
  )

  base_target <- base$raw$target
  base_target$source_id <- data$source_id[base_target$source_row]
  permuted_target <- permuted$raw$target
  permuted_target$source_id <-
    permuted_data$source_id[permuted_target$source_row]
  aligned <- merge(
    base_target[, c("source_id", "position", "target")],
    permuted_target[, c("source_id", "position", "target")],
    by = "source_id",
    suffixes = c("_base", "_permuted")
  )

  expect_equal(
    aligned$position_permuted,
    aligned$position_base,
    tolerance = 2e-14
  )
  expect_equal(
    aligned$target_permuted,
    aligned$target_base,
    tolerance = 2e-14
  )
  expect_equal(
    permuted$raw$model$coefficients,
    base$raw$model$coefficients,
    tolerance = 2e-14
  )
})

test_that("KE phase 3 WLS equals integer frequency expansion", {
  data <- .phase3_ke_data()
  weighted <- .phase3_ke_fit(data)
  expanded <- data[
    rep(seq_len(nrow(data)), data$weight),
    ,
    drop = FALSE
  ]
  expanded$expanded_id <- seq_len(nrow(expanded))
  expanded_fit <- suppressWarnings(ke_decomp(
    data = expanded,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses_tied",
    index_type = "rank",
    correction = "generalized",
    use_svy = FALSE,
    vce_method = "linearized",
    quiet = TRUE
  ))

  expect_equal(nrow(expanded), sum(data$weight))
  expect_equal(
    expanded_fit$results_overall$Estimate,
    weighted$results_overall$Estimate,
    tolerance = 2e-14
  )
  expect_equal(
    expanded_fit$raw$model$coefficients,
    weighted$raw$model$coefficients,
    tolerance = 2e-14
  )
  expect_equal(
    expanded_fit$model_metrics$R_Squared,
    weighted$model_metrics$R_Squared,
    tolerance = 2e-14
  )
})

test_that("KE phase 3 validates the level-weight formula in Appendix A", {
  data <- .phase3_ke_data()
  fit <- .phase3_ke_fit(data, "income", "level")
  p <- data$weight / sum(data$weight)
  mean_health <- sum(p * data$health)
  mean_income <- sum(p * data$income)
  relative_income <- data$income / mean_income
  expected_target <- relative_income * data$health - mean_health
  appendix_weight <- nrow(data) * p * (relative_income - 1)
  appendix_index <- sum(appendix_weight * data$health) /
    nrow(data)

  raw_target <- fit$raw$target
  expected_sorted <- expected_target[raw_target$source_row]
  expect_equal(
    raw_target$target,
    expected_sorted,
    tolerance = 2e-14
  )
  expect_lt(
    abs(fit$results_overall$Estimate - appendix_index),
    2e-14
  )
  expect_lt(
    abs(
      fit$results_overall$Estimate -
        (-0.0006805923199052794)
    ),
    2e-14
  )
  expect_equal(
    unname(fit$raw$model$coefficients),
    c(
      -0.010020441967182330,
      0.0079411602462004065,
      -0.013541637816091120
    ),
    tolerance = 2e-14
  )

  scaled_weight <- data
  scaled_weight$weight <- 10 * scaled_weight$weight
  scaled_income <- data
  scaled_income$income <- 10 * scaled_income$income
  weight_fit <- .phase3_ke_fit(
    scaled_weight, "income", "level"
  )
  income_fit <- .phase3_ke_fit(
    scaled_income, "income", "level"
  )

  expect_lt(
    abs(
      weight_fit$results_overall$Estimate -
        fit$results_overall$Estimate
    ),
    2e-14
  )
  expect_lt(
    abs(
      income_fit$results_overall$Estimate -
        fit$results_overall$Estimate
    ),
    2e-14
  )
  expect_equal(
    weight_fit$raw$model$coefficients,
    fit$raw$model$coefficients,
    tolerance = 2e-14
  )
  expect_equal(
    income_fit$raw$model$coefficients,
    fit$raw$model$coefficients,
    tolerance = 2e-14
  )
})

test_that("KE phase 3 handles the all-tied rank limit coherently", {
  data <- .phase3_ke_data()
  fit <- .phase3_ke_fit(data, "ses_all_tied", "rank")

  expect_equal(
    fit$raw$target$fractional_rank,
    rep(0.5, nrow(data)),
    tolerance = 2e-14
  )
  expect_equal(
    fit$results_overall$Estimate,
    0,
    tolerance = 2e-14
  )
  expect_gt(stats::var(fit$raw$target$target), 0)
  expect_equal(
    fit$raw$index$regression_at_covariate_means,
    0,
    tolerance = 2e-14
  )
})
