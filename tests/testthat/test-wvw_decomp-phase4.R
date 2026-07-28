.phase4_wvw_data <- function() {
  id <- seq_len(72L)
  tie_sizes <- c(6L, 9L, 11L, 14L, 17L, 15L)
  data <- data.frame(
    id = id,
    ses_unique = ((37L * id) %% 73L) + 1L,
    ses_tied = rep(c(1, 2, 4, 7, 10, 15), times = tie_sizes),
    ses_extreme = rep(c(0, 1), times = c(31L, 41L)),
    ses_all = 1,
    x1 = 1 + ((7L * id) %% 17L) / 4,
    x2 = 2 + ((11L * id) %% 19L) / 5,
    x3 = as.integer((id %% 5L) %in% c(0L, 1L)),
    weight = 1L + ((5L * id) %% 7L),
    permutation_key = ((29L * id) %% 73L)
  )
  data$noise <- (((23L * id) %% 29L) - 14) / 18
  data$y <- with(
    data,
    4.5 + 0.65 * x1 - 0.38 * x2 + 0.9 * x3 + noise
  )
  data$weight_zero <- data$weight
  data$weight_zero[c(1L, 8L, 17L, 29L, 46L, 64L)] <- 0
  data
}

.phase4_grouped_rank <- function(rank_variable, weight) {
  group_weight <- tapply(weight, rank_variable, sum)
  group_values <- as.numeric(names(group_weight))
  group_order <- order(group_values)
  group_weight <- group_weight[group_order]
  group_values <- group_values[group_order]
  group_rank <- (
    c(0, head(cumsum(group_weight), -1L)) +
      0.5 * group_weight
  ) / sum(group_weight)
  as.numeric(group_rank[match(rank_variable, group_values)])
}

.phase4_fit <- function(
    data, ses_var = NULL, precalc_rank_var = NULL,
    weight_var = "weight"
) {
  wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = c("x1", "x2", "x3"),
    ses_var = ses_var,
    precalc_rank_var = precalc_rank_var,
    correction = "standard",
    model_type = "ols",
    weight_var = weight_var,
    vce_method = "linearized",
    quiet = TRUE
  )
}

test_that("phase 4 grouped weighted ties match the frozen Stata result", {
  data <- .phase4_wvw_data()
  fit <- .phase4_fit(data, ses_var = "ses_tied")
  raw <- fit$raw
  rank_by_id <- raw$rank$fractional_rank[
    match(data$id, raw$rank$source_row)
  ]
  expected_rank <- .phase4_grouped_rank(data$ses_tied, data$weight)

  expect_equal(rank_by_id, expected_rank, tolerance = 1e-15)
  expect_equal(raw$index$rank_mean, 0.5, tolerance = 1e-15)
  expect_equal(
    raw$index$total,
    -0.007143594349996599,
    tolerance = 1e-12
  )
  expect_equal(
    raw$index$explained,
    -0.0013810195946953702,
    tolerance = 1e-12
  )
  expect_equal(
    raw$index$residual,
    -0.005762574755301229,
    tolerance = 1e-12
  )
  expect_equal(
    raw$determinants$contribution[
      match(c("x1", "x2", "x3"), raw$determinants$column)
    ],
    c(
      0.0024293383049555509,
      -0.00089019199701152242,
      -0.0029201659026393987
    ),
    tolerance = 1e-12
  )
})

test_that("phase 4 grouped ties are invariant to row permutations", {
  data <- .phase4_wvw_data()
  permuted <- data[order(data$permutation_key), , drop = FALSE]
  base <- .phase4_fit(data, ses_var = "ses_tied")$raw
  reordered <- .phase4_fit(permuted, ses_var = "ses_tied")$raw

  base_rank <- data.frame(
    id = data$id[base$rank$source_row],
    rank = base$rank$fractional_rank
  )
  reordered_rank <- data.frame(
    id = permuted$id[reordered$rank$source_row],
    rank = reordered$rank$fractional_rank
  )
  rank_pair <- merge(
    base_rank,
    reordered_rank,
    by = "id",
    suffixes = c("_base", "_permuted")
  )

  expect_equal(
    rank_pair$rank_permuted,
    rank_pair$rank_base,
    tolerance = 1e-15
  )
  expect_equal(
    c(
      reordered$index$total,
      reordered$index$explained,
      reordered$index$residual
    ),
    c(
      base$index$total,
      base$index$explained,
      base$index$residual
    ),
    tolerance = 1e-12
  )
  expect_equal(
    reordered$estimates$detailed,
    base$estimates$detailed,
    tolerance = 1e-12
  )
})

test_that("phase 4 canonical precalculated rank reproduces internal rank", {
  data <- .phase4_wvw_data()
  data$rank_precalculated <- .phase4_grouped_rank(
    data$ses_tied,
    data$weight
  )
  internal <- .phase4_fit(data, ses_var = "ses_tied")$raw
  supplied <- .phase4_fit(
    data,
    precalc_rank_var = "rank_precalculated"
  )$raw
  supplied_rank_by_id <- supplied$rank$fractional_rank[
    match(data$id, supplied$rank$source_row)
  ]

  expect_equal(
    supplied_rank_by_id,
    data$rank_precalculated,
    tolerance = 0
  )
  expect_identical(
    supplied$settings$rank_method,
    "precalculated rank supplied by user"
  )
  expect_equal(
    c(
      supplied$index$total,
      supplied$index$explained,
      supplied$index$residual
    ),
    c(
      internal$index$total,
      internal$index$explained,
      internal$index$residual
    ),
    tolerance = 1e-12
  )
  expect_equal(
    supplied$estimates$detailed,
    internal$estimates$detailed,
    tolerance = 1e-12
  )
})

test_that("phase 4 all-tied SES produces rank one half and zero inequality", {
  data <- .phase4_wvw_data()
  raw <- .phase4_fit(data, ses_var = "ses_all")$raw

  expect_equal(raw$rank$fractional_rank, rep(0.5, nrow(data)))
  expect_equal(raw$index$total, 0, tolerance = 1e-15)
  expect_equal(raw$index$explained, 0, tolerance = 1e-15)
  expect_equal(raw$index$residual, 0, tolerance = 1e-15)
  expect_equal(
    raw$determinants$contribution,
    rep(0, nrow(raw$determinants)),
    tolerance = 1e-15
  )
})

test_that("phase 4 zero-weight rows do not affect point estimates", {
  data <- .phase4_wvw_data()
  data$analysis_weight <- data$weight_zero
  retained <- .phase4_fit(
    data,
    ses_var = "ses_tied",
    weight_var = "analysis_weight"
  )$raw
  positive <- .phase4_fit(
    data[data$analysis_weight > 0, , drop = FALSE],
    ses_var = "ses_tied",
    weight_var = "analysis_weight"
  )$raw

  expect_equal(retained$sample$analytic_n, 66)
  expect_equal(positive$sample$analytic_n, 66)
  expect_equal(
    retained$sample$excluded_source_rows$zero_sampling_weight,
    c(1L, 8L, 17L, 29L, 46L, 64L)
  )
  expect_equal(
    c(
      retained$index$total,
      retained$index$explained,
      retained$index$residual
    ),
    c(
      positive$index$total,
      positive$index$explained,
      positive$index$residual
    ),
    tolerance = 1e-12
  )
  expect_equal(
    retained$estimates$detailed,
    positive$estimates$detailed,
    tolerance = 1e-12
  )
})

test_that("phase 4 precalculated ranks reject invalid fractional values", {
  data <- .phase4_wvw_data()
  data$rank_precalculated <- .phase4_grouped_rank(
    data$ses_tied,
    data$weight
  )

  data$rank_precalculated[1] <- Inf
  expect_error(
    .phase4_fit(data, precalc_rank_var = "rank_precalculated"),
    "Non-finite values detected",
    fixed = TRUE
  )

  data$rank_precalculated[1] <- -0.01
  expect_error(
    .phase4_fit(data, precalc_rank_var = "rank_precalculated"),
    "closed interval [0, 1]",
    fixed = TRUE
  )

  data$rank_precalculated[1] <- 1.01
  expect_error(
    .phase4_fit(data, precalc_rank_var = "rank_precalculated"),
    "closed interval [0, 1]",
    fixed = TRUE
  )
})
