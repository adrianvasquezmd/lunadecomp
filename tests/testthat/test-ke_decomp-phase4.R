.phase4_ke_data <- function() {
  n <- 72L
  source_id <- seq_len(n)
  tie_sizes <- c(7L, 9L, 11L, 13L, 15L, 17L)
  tie_group <- rep(seq_along(tie_sizes), times = tie_sizes)
  x1 <- ((11 * source_id) %% 27) - 13 + source_id / 19
  x2 <- cos(source_id / 5) + 0.17 * (source_id %% 4)
  weight <- ((5 * source_id) %% 7) + 1
  latent_health <- 0.34 +
    0.006 * x1 -
    0.038 * x2 +
    0.08 * tie_group +
    0.021 * sin(source_id / 4)
  health <- (latent_health - min(latent_health)) /
    (max(latent_health) - min(latent_health))
  health_binary <- as.integer(
    latent_health > stats::weighted.mean(latent_health, weight)
  )
  input_order <- c(
    seq(3, n, by = 3),
    seq(1, n, by = 3),
    seq(2, n, by = 3)
  )

  data.frame(
    source_id = source_id,
    health = health,
    health_binary = health_binary,
    health_complement = 1 - health_binary,
    health_affine = 2 + 3 * health,
    health_scale10 = 10 * health,
    ses_tied = c(1, 2, 4, 7, 11, 18)[tie_group],
    income = 380 + 10 * source_id + 28 * sin(source_id / 6),
    x1 = x1,
    x2 = x2,
    weight = weight,
    stringsAsFactors = FALSE
  )[input_order, , drop = FALSE]
}

.phase4_ke_fit <- function(
    data, outcome = "health", index_type = "rank",
    correction = "generalized", dep_min = NULL, dep_max = NULL) {
  suppressWarnings(ke_decomp(
    data = data,
    dep_var = outcome,
    indep_vars = c("x1", "x2"),
    ses_var = if (index_type == "rank") "ses_tied" else "income",
    index_type = index_type,
    correction = correction,
    dep_min = dep_min,
    dep_max = dep_max,
    use_svy = TRUE,
    weight_var = "weight",
    vce_method = "linearized",
    quiet = TRUE
  ))
}

test_that("KE phase 4 rank corrections obey their defining scalars", {
  data <- .phase4_ke_data()
  generalized <- .phase4_ke_fit(data)
  standard <- .phase4_ke_fit(
    data, correction = "standard"
  )
  erreygers <- .phase4_ke_fit(
    data,
    correction = "erreygers",
    dep_min = 0,
    dep_max = 1
  )
  wagstaff <- .phase4_ke_fit(
    data,
    correction = "wagstaff",
    dep_min = 0,
    dep_max = 1
  )
  mean_health <- generalized$raw$index$health_mean
  generalized_index <- generalized$results_overall$Estimate

  expect_lt(
    abs(
      standard$results_overall$Estimate -
        generalized_index / mean_health
    ),
    5e-14
  )
  expect_lt(
    abs(
      erreygers$results_overall$Estimate -
        4 * generalized_index
    ),
    5e-14
  )
  expect_lt(
    abs(
      wagstaff$results_overall$Estimate -
        generalized_index /
          (mean_health * (1 - mean_health))
    ),
    5e-14
  )
  expect_equal(
    c(
      generalized$raw$index$correction_scale,
      standard$raw$index$correction_scale,
      erreygers$raw$index$correction_scale,
      wagstaff$raw$index$correction_scale
    ),
    c(
      1,
      1 / mean_health,
      4,
      1 / (mean_health * (1 - mean_health))
    ),
    tolerance = 2e-14
  )

  # Frozen full-precision Stata conindex/composed-regression benchmarks.
  expect_lt(
    abs(generalized_index - 0.13155352552540198),
    1e-14
  )
  expect_lt(
    abs(
      standard$results_overall$Estimate -
        0.23639787436079196
    ),
    1e-14
  )
  expect_lt(
    abs(
      erreygers$results_overall$Estimate -
        0.52621410210160791
    ),
    1e-14
  )
  expect_lt(
    abs(
      wagstaff$results_overall$Estimate -
        0.53301828494391179
    ),
    1e-14
  )

  expect_equal(
    standard$raw$model$coefficients,
    generalized$raw$model$coefficients / mean_health,
    tolerance = 2e-14
  )
  expect_equal(
    erreygers$raw$model$coefficients,
    4 * generalized$raw$model$coefficients,
    tolerance = 2e-14
  )
})

test_that("KE phase 4 level corrections follow Kessels-Erreygers", {
  data <- .phase4_ke_data()
  generalized <- .phase4_ke_fit(
    data, index_type = "level"
  )
  standard <- .phase4_ke_fit(
    data,
    index_type = "level",
    correction = "standard"
  )
  erreygers <- .phase4_ke_fit(
    data,
    index_type = "level",
    correction = "erreygers",
    dep_min = 0,
    dep_max = 1
  )
  mean_health <- generalized$raw$index$health_mean

  expect_lt(
    abs(
      standard$results_overall$Estimate -
        generalized$results_overall$Estimate / mean_health
    ),
    5e-14
  )
  expect_lt(
    abs(
      erreygers$results_overall$Estimate -
        generalized$results_overall$Estimate
    ),
    5e-14
  )
  expect_equal(
    erreygers$raw$index$correction_scale,
    1,
    tolerance = 0
  )
  expect_lt(
    abs(
      generalized$results_overall$Estimate -
        0.062318707093639
    ),
    1e-14
  )
  expect_lt(
    abs(
      standard$results_overall$Estimate -
        0.11198491132040719
    ),
    1e-14
  )
})

test_that("KE phase 4 corrected rank indices match rineq", {
  data <- .phase4_ke_data()
  correction_types <- c(
    generalized = "CIg",
    standard = "CI",
    erreygers = "CIc",
    wagstaff = "CIw"
  )
  for (correction in names(correction_types)) {
    bounded <- correction %in% c("erreygers", "wagstaff")
    fit <- .phase4_ke_fit(
      data,
      correction = correction,
      dep_min = if (bounded) 0 else NULL,
      dep_max = if (bounded) 1 else NULL
    )
    reference <- rineq::ci(
      ineqvar = data$ses_tied,
      outcome = data$health,
      weights = data$weight,
      type = correction_types[[correction]],
      method = "direct",
      df_correction = FALSE,
      rank_function = rineq::rank_gwt
    )$concentration_index
    expect_lt(
      abs(fit$results_overall$Estimate - reference),
      1e-14
    )
  }
})

test_that("KE phase 4 requires valid theoretical bounded-outcome limits", {
  data <- .phase4_ke_data()

  expect_error(
    .phase4_ke_fit(data, correction = "erreygers"),
    "requires explicit theoretical health bounds"
  )
  expect_error(
    .phase4_ke_fit(
      data,
      correction = "wagstaff",
      dep_min = 0
    ),
    "requires explicit theoretical health bounds"
  )
  expect_error(
    .phase4_ke_fit(
      data,
      correction = "erreygers",
      dep_min = 0.1,
      dep_max = 1
    ),
    "outside the declared theoretical bounds"
  )
  expect_error(
    .phase4_ke_fit(
      data,
      correction = "erreygers",
      dep_min = 1,
      dep_max = 0
    ),
    "dep_min >= dep_max"
  )
  expect_error(
    .phase4_ke_fit(
      data,
      index_type = "level",
      correction = "wagstaff",
      dep_min = 0,
      dep_max = 1
    ),
    "defined only for index_type = 'rank'"
  )
})

test_that("KE phase 4 distinguishes index invariance from coefficient sensitivity", {
  data <- .phase4_ke_data()
  rank_erreygers <- .phase4_ke_fit(
    data,
    correction = "erreygers",
    dep_min = 0,
    dep_max = 1
  )
  affine_erreygers <- .phase4_ke_fit(
    data,
    outcome = "health_affine",
    correction = "erreygers",
    dep_min = 2,
    dep_max = 5
  )
  scaled_generalized <- .phase4_ke_fit(
    data,
    outcome = "health_scale10",
    correction = "generalized"
  )
  scaled_standard <- .phase4_ke_fit(
    data,
    outcome = "health_scale10",
    correction = "standard"
  )
  base_generalized <- .phase4_ke_fit(data)
  base_standard <- .phase4_ke_fit(
    data,
    correction = "standard"
  )

  expect_lt(
    abs(
      affine_erreygers$results_overall$Estimate -
        rank_erreygers$results_overall$Estimate
    ),
    5e-14
  )
  expect_gt(
    max(abs(
      affine_erreygers$raw$model$coefficients -
        rank_erreygers$raw$model$coefficients
    )),
    1e-8
  )
  expect_lt(
    abs(
      scaled_generalized$results_overall$Estimate -
        10 * base_generalized$results_overall$Estimate
    ),
    5e-14
  )
  expect_equal(
    scaled_generalized$raw$model$coefficients,
    10 * base_generalized$raw$model$coefficients,
    tolerance = 2e-14
  )
  expect_lt(
    abs(
      scaled_standard$results_overall$Estimate -
        base_standard$results_overall$Estimate
    ),
    5e-14
  )
  expect_equal(
    scaled_standard$raw$model$coefficients,
    base_standard$raw$model$coefficients,
    tolerance = 2e-14
  )
})

test_that("KE phase 4 mirror corrections flip the aggregate index", {
  data <- .phase4_ke_data()
  for (correction in c("erreygers", "wagstaff")) {
    original <- .phase4_ke_fit(
      data,
      outcome = "health_binary",
      correction = correction,
      dep_min = 0,
      dep_max = 1
    )
    complement <- .phase4_ke_fit(
      data,
      outcome = "health_complement",
      correction = correction,
      dep_min = 0,
      dep_max = 1
    )
    expect_lt(
      abs(
        complement$results_overall$Estimate +
          original$results_overall$Estimate
      ),
      5e-14
    )
    expect_gt(
      max(abs(
        complement$raw$model$coefficients +
          original$raw$model$coefficients
      )),
      1e-8
    )
  }
})
