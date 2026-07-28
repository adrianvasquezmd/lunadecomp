.phase6_wvw_data <- function() {
  id <- seq_len(96L)
  data <- data.frame(
    id = id,
    ses = rep(c(1, 2, 4, 7, 10, 14, 19, 25), each = 12L),
    x1 = 1 + ((7L * id) %% 17L) / 4,
    x2 = 2 + ((11L * id) %% 19L) / 5,
    x3 = as.integer((id %% 5L) %in% c(0L, 1L)),
    region = 1L + ((3L * id) %% 4L),
    education = 1L + ((5L * id) %% 3L),
    weight = 1L + ((5L * id) %% 7L)
  )
  data$noise <- (((23L * id) %% 29L) - 14) / 20
  region_effect <- c(0, 0.45, -0.30, 0.70)[data$region]
  education_effect <- c(0, 0.28, -0.18)[data$education]
  data$y <- with(
    data,
    5.5 + 0.62 * x1 - 0.34 * x2 + 0.82 * x3 +
      region_effect + education_effect + noise
  )
  for (level in 1:4) {
    data[[paste0("region", level)]] <-
      as.integer(data$region == level)
  }
  data
}

.phase6_wvw_fit <- function(
    data, predictors = c("x1", "x2", "x3", "region"),
    groupings = list(), correction = "standard", relax = FALSE
) {
  wvw_decomp(
    data = data,
    dep_var = "y",
    indep_vars = predictors,
    ses_var = "ses",
    groupings = groupings,
    correction = correction,
    model_type = "ols",
    weight_var = "weight",
    vce_method = "linearized",
    relax = relax,
    quiet = TRUE
  )
}

.phase6_factor_data <- function(
    levels = c(1, 2, 3, 4)
) {
  data <- .phase6_wvw_data()
  data$region <- factor(data$region, levels = levels)
  data$education <- factor(data$education, levels = 1:3)
  data
}

test_that("phase 6 factors match frozen Stata decomposition values", {
  base1 <- .phase6_wvw_fit(.phase6_factor_data())$raw
  base2 <- .phase6_wvw_fit(
    .phase6_factor_data(c(2, 1, 3, 4))
  )$raw
  two_factors <- .phase6_wvw_fit(
    .phase6_factor_data(),
    c("x1", "x2", "x3", "region", "education")
  )$raw

  expect_equal(
    unlist(base1$index[c("total", "explained", "residual")]),
    c(
      total = -0.0037421393930548064,
      explained = -0.0019850811585981516,
      residual = -0.0017570582344563755
    ),
    tolerance = 1e-12
  )
  expect_equal(
    unlist(two_factors$index[c("total", "explained", "residual")]),
    c(
      total = -0.0037421393930548064,
      explained = -0.0016196182948779491,
      residual = -0.002122521098176578
    ),
    tolerance = 1e-12
  )

  factor_block <- function(raw, term) {
    sum(
      raw$determinants$contribution[
        raw$determinants$source_term == term
      ]
    )
  }
  expect_equal(
    factor_block(base1, "region"),
    0.000068997014880893139,
    tolerance = 1e-12
  )
  expect_equal(
    factor_block(base2, "region"),
    factor_block(base1, "region"),
    tolerance = 1e-12
  )
  expect_equal(
    factor_block(two_factors, "education"),
    0.00035222436740586956,
    tolerance = 1e-12
  )
})

test_that("phase 6 factor blocks are invariant to the omitted category", {
  base1 <- .phase6_wvw_fit(.phase6_factor_data())$raw
  base2 <- .phase6_wvw_fit(
    .phase6_factor_data(c(2, 1, 3, 4))
  )$raw

  expect_equal(
    base2$model$fitted,
    base1$model$fitted,
    tolerance = 1e-12
  )
  expect_equal(
    unlist(base2$index[c("total", "explained", "residual")]),
    unlist(base1$index[c("total", "explained", "residual")]),
    tolerance = 1e-12
  )

  contribution_by_term <- function(raw, term) {
    sum(
      raw$determinants$contribution[
        raw$determinants$source_term == term
      ]
    )
  }
  expect_equal(
    contribution_by_term(base2, "region"),
    contribution_by_term(base1, "region"),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(
    base2$model$coefficients[
      grepl("^region", names(base2$model$coefficients))
    ],
    base1$model$coefficients[
      grepl("^region", names(base1$model$coefficients))
    ]
  )))
})

test_that("phase 6 formula factors equal explicit manual dummies", {
  factor_fit <- .phase6_wvw_fit(.phase6_factor_data())$raw
  manual_fit <- .phase6_wvw_fit(
    .phase6_wvw_data(),
    c("x1", "x2", "x3", "region2", "region3", "region4")
  )$raw

  expect_equal(
    manual_fit$model$fitted,
    factor_fit$model$fitted,
    tolerance = 1e-12
  )
  expect_equal(
    manual_fit$model$coefficients,
    factor_fit$model$coefficients,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(
    manual_fit$estimates$detailed,
    factor_fit$estimates$detailed,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_equal(
    unlist(manual_fit$index[c("total", "explained", "residual")]),
    unlist(factor_fit$index[c("total", "explained", "residual")]),
    tolerance = 1e-12
  )
})

test_that("phase 6 custom groups form a complete non-overlapping partition", {
  data <- .phase6_factor_data()
  predictors <- c("x1", "x2", "x3", "region", "education")
  full_fit <- .phase6_wvw_fit(
    data,
    predictors,
    groupings = list(
      Continuous = c("x1", "x2", "x3"),
      Categorical = c("region", "education")
    )
  )
  partial_fit <- .phase6_wvw_fit(
    data,
    predictors,
    groupings = list(Context = c("region", "education"))
  )
  full <- full_fit$raw
  partial <- partial_fit$raw

  for (raw in list(full, partial)) {
    assignment <- unlist(raw$groups, use.names = FALSE)
    expect_equal(
      sum(raw$estimates$grouped),
      raw$index$explained,
      tolerance = 1e-12
    )
    expect_length(assignment, ncol(raw$model$matrix) - 1L)
    expect_identical(anyDuplicated(assignment), 0L)
    expect_setequal(
      assignment,
      setdiff(colnames(raw$model$matrix), "(Intercept)")
    )
  }
  expect_match(
    partial_fit$diagnostics$groupings_uncovered,
    "retained automatically",
    fixed = TRUE
  )
})

test_that("phase 6 rejects invalid custom grouping specifications", {
  data <- .phase6_factor_data()
  predictors <- c("x1", "x2", "x3", "region")

  expect_error(
    .phase6_wvw_fit(
      data,
      predictors,
      groupings = list(Bad = "does_not_exist")
    ),
    "unknown grouping token",
    fixed = TRUE
  )
  expect_error(
    .phase6_wvw_fit(
      data,
      predictors,
      groupings = list(A = "region", B = "region2")
    ),
    "assigned to more than one group",
    fixed = TRUE
  )

  ambiguous <- data
  ambiguous$region_extra <-
    ((13L * ambiguous$id) %% 23L) / 7
  expect_error(
    .phase6_wvw_fit(
      ambiguous,
      c(predictors, "region_extra"),
      groupings = list(Ambiguous = "reg")
    ),
    "ambiguous prefix",
    fixed = TRUE
  )
  expect_error(
    .phase6_wvw_fit(
      data,
      predictors,
      groupings = list(A = "x1", A = "x2")
    ),
    "unique, non-empty names",
    fixed = TRUE
  )
})

test_that("phase 6 missing required values define one common sample", {
  data <- .phase6_factor_data()
  data$y[1L] <- NA_real_
  data$x1[2L] <- NA_real_
  data$region[3L] <- NA
  data$ses[4L] <- NA_real_
  data$weight[5L] <- NA_real_
  data$x2[6L] <- NA_real_

  raw <- .phase6_wvw_fit(data)$raw

  expect_identical(raw$sample$analytic_n, 90L)
  expect_setequal(
    raw$sample$excluded_source_rows[["missing_required_values"]],
    1:6
  )
  expect_setequal(raw$sample$analytic_source_rows, 7:96)
  expect_equal(
    unlist(raw$index[c("total", "explained", "residual")]),
    c(
      total = -0.0015191160108172664,
      explained = -0.0012472004278683076,
      residual = -0.00027191558294982119
    ),
    tolerance = 1e-12
  )
})

test_that("phase 6 zero weights are excluded and audited", {
  data <- .phase6_factor_data()
  data$weight[c(2L, 7L, 31L)] <- 0
  with_zero_fit <- .phase6_wvw_fit(data)
  with_zero <- with_zero_fit$raw
  without_zero <- .phase6_wvw_fit(
    data[data$weight > 0, , drop = FALSE]
  )$raw

  expect_identical(with_zero$sample$analytic_n, 93L)
  expect_setequal(
    with_zero$sample$excluded_source_rows$zero_sampling_weight,
    c(2L, 7L, 31L)
  )
  expect_equal(
    unlist(with_zero$index[c("total", "explained", "residual")]),
    unlist(without_zero$index[c("total", "explained", "residual")]),
    tolerance = 1e-12
  )
  expect_match(
    with_zero_fit$diagnostics$zero_weights,
    "were excluded",
    fixed = TRUE
  )
})

test_that("phase 6 validates the complete sampling-weight domain", {
  data <- .phase6_factor_data()

  negative <- data
  negative$weight[1L] <- -1
  expect_error(
    .phase6_wvw_fit(negative),
    "cannot be negative",
    fixed = TRUE
  )

  nonfinite <- data
  nonfinite$weight[1L] <- Inf
  expect_error(
    .phase6_wvw_fit(nonfinite),
    "Non-finite values detected",
    fixed = TRUE
  )

  all_zero <- data
  all_zero$weight[] <- 0
  expect_error(
    .phase6_wvw_fit(all_zero),
    "No positive sampling weight remains",
    fixed = TRUE
  )

  character_weight <- data
  character_weight$weight <- as.character(character_weight$weight)
  expect_equal(
    unlist(
      .phase6_wvw_fit(character_weight)$raw$index[
        c("total", "explained", "residual")
      ]
    ),
    unlist(
      .phase6_wvw_fit(data)$raw$index[
        c("total", "explained", "residual")
      ]
    ),
    tolerance = 1e-12
  )

  invalid_character <- character_weight
  invalid_character$weight[1L] <- "invalid"
  expect_error(
    .phase6_wvw_fit(invalid_character),
    "Numeric conversion introduced NA",
    fixed = TRUE
  )
})

test_that("phase 6 exact collinearity always stops the decomposition", {
  data <- .phase6_factor_data()
  data$x4 <- 2 * data$x1

  expect_error(
    .phase6_wvw_fit(
      data,
      c("x1", "x2", "x3", "region", "x4")
    ),
    "Exact collinearity detected",
    fixed = TRUE
  )
  expect_error(
    .phase6_wvw_fit(
      data,
      c("x1", "x2", "x3", "region", "x4"),
      relax = TRUE
    ),
    "Exact collinearity detected",
    fixed = TRUE
  )
})

test_that("phase 6 rejects malformed predictors and non-finite matrices", {
  data <- .phase6_factor_data()

  expect_error(
    .phase6_wvw_fit(
      data,
      c("x1", "x2", "x3", "region", "not_a_variable")
    ),
    "variables not found in data: not_a_variable",
    fixed = TRUE
  )
  expect_error(
    .phase6_wvw_fit(data, c("x1", "x1", "x3", "region")),
    "duplicated variable names",
    fixed = TRUE
  )

  one_level <- data
  one_level$region <- factor("only")
  expect_error(
    .phase6_wvw_fit(one_level),
    "fewer than two observed levels",
    fixed = TRUE
  )

  nonfinite <- data
  nonfinite$x1[1L] <- Inf
  expect_error(
    .phase6_wvw_fit(nonfinite),
    "Non-finite values detected in model-matrix columns",
    fixed = TRUE
  )
})

test_that("phase 6 standard correction requires a positive outcome mean", {
  data <- .phase6_factor_data()
  weighted_mean <- stats::weighted.mean(data$y, data$weight)
  data$y <- data$y - weighted_mean

  expect_error(
    .phase6_wvw_fit(data),
    "requires a strictly positive mean",
    fixed = TRUE
  )

  generalized <- .phase6_wvw_fit(
    data,
    correction = "generalized"
  )
  expect_true(is.finite(generalized$raw$index$total))
  expect_equal(
    generalized$raw$index$total,
    generalized$raw$index$explained +
      generalized$raw$index$residual,
    tolerance = 1e-12
  )
})

test_that("phase 6 a zero-mean determinant has undefined display CI only", {
  data <- .phase6_factor_data()
  data$x_centered <- data$x1 -
    stats::weighted.mean(data$x1, data$weight)
  fit <- .phase6_wvw_fit(
    data,
    c("x_centered", "x2", "x3", "region")
  )$raw
  determinant <- fit$determinants[
    fit$determinants$column == "x_centered",
    ,
    drop = FALSE
  ]

  expect_true(is.na(determinant$concentration_index))
  expect_true(is.finite(determinant$contribution))
  expect_equal(
    fit$index$total,
    fit$index$explained + fit$index$residual,
    tolerance = 1e-12
  )
})

test_that("phase 6 transformations and interactions remain outside the API", {
  data <- .phase6_factor_data()

  expect_error(
    .phase6_wvw_fit(
      data,
      c("x1", "x2:x3", "region")
    ),
    "variables not found in data: x2:x3",
    fixed = TRUE
  )
  expect_error(
    .phase6_wvw_fit(
      data,
      c("x1", "I(x2^2)", "region")
    ),
    "variables not found in data: I(x2^2)",
    fixed = TRUE
  )
})
