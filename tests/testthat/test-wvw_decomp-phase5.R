.phase5_wvw_data <- function() {
  id <- seq_len(80L)
  data <- data.frame(
    id = id,
    ses = rep(c(1, 2, 4, 7, 10, 14, 19, 25), each = 10L),
    x1 = 1 + ((7L * id) %% 17L) / 4,
    x2 = 2 + ((11L * id) %% 19L) / 5,
    x3 = as.integer((id %% 5L) %in% c(0L, 1L)),
    weight = 1L + ((5L * id) %% 7L)
  )
  data$noise <- (((23L * id) %% 29L) - 14) / 20
  data$y_positive <- with(
    data,
    6 + 0.58 * x1 - 0.31 * x2 + 0.75 * x3 + noise
  )
  bounded_fraction <- ((17L * id) %% 83L) / 82
  bounded_fraction[1L] <- 0
  bounded_fraction[80L] <- 1
  data$y_bounded <- 10 + 20 * bounded_fraction
  data$y_bounded_short <- 40 - data$y_bounded
  interior_fraction <- 0.15 +
    0.70 * (((19L * id) %% 83L) / 82)
  data$y_bounded_interior <- 10 + 20 * interior_fraction
  data$y_binary <- as.integer(((13L * id) %% 17L) < 7L)
  data$y_binary_short <- 1L - data$y_binary
  data
}

.phase5_fit <- function(
    data, outcome, correction, lower = NULL, upper = NULL,
    shortfall = FALSE
) {
  wvw_decomp(
    data = data,
    dep_var = outcome,
    indep_vars = c("x1", "x2", "x3"),
    ses_var = "ses",
    correction = correction,
    model_type = "ols",
    dep_min = lower,
    dep_max = upper,
    is_shortfall = shortfall,
    weight_var = "weight",
    vce_method = "linearized",
    quiet = TRUE
  )
}

test_that("phase 5 four corrections match frozen Stata totals", {
  data <- .phase5_wvw_data()
  fits <- list(
    standard = .phase5_fit(
      data, "y_positive", "standard"
    ),
    generalized = .phase5_fit(
      data, "y_positive", "generalized"
    ),
    erreygers = .phase5_fit(
      data, "y_bounded", "erreygers", 10, 30
    ),
    wagstaff = .phase5_fit(
      data, "y_bounded", "wagstaff", 10, 30
    )
  )
  expected <- list(
    standard = c(
      total = -0.0023070089519601094,
      explained = -0.0013939696346889263,
      residual = -0.00091303931727116316
    ),
    generalized = c(
      total = -0.015825000000001616,
      explained = -0.0095619782706997756,
      residual = -0.0062630217293002085
    ),
    erreygers = c(
      total = 0.0693592797256085,
      explained = 0.019080317988605995,
      residual = 0.05027896173700408
    ),
    wagstaff = c(
      total = 0.069362055668137926,
      explained = 0.019081081633591391,
      residual = 0.050280974034546719
    )
  )

  for (correction in names(fits)) {
    observed <- unlist(
      fits[[correction]]$raw$index[
        c("total", "explained", "residual")
      ]
    )
    expect_equal(
      unname(observed),
      unname(expected[[correction]]),
      tolerance = 1e-12,
      info = correction
    )
    expect_equal(
      fits[[correction]]$raw$index$total,
      fits[[correction]]$raw$index$explained +
        fits[[correction]]$raw$index$residual,
      tolerance = 1e-12,
      info = correction
    )
  }
})

test_that("phase 5 correction scale applies to every component", {
  data <- .phase5_wvw_data()
  outcome <- data$y_bounded
  mean_y <- stats::weighted.mean(outcome, data$weight)
  expected_scale <- c(
    standard = 1 / mean_y,
    generalized = 1,
    erreygers = 4 / 20,
    wagstaff = 20 / ((30 - mean_y) * (mean_y - 10))
  )
  fits <- lapply(
    names(expected_scale),
    function(correction) {
      .phase5_fit(
        data,
        "y_bounded",
        correction,
        if (correction %in% c("erreygers", "wagstaff")) 10 else NULL,
        if (correction %in% c("erreygers", "wagstaff")) 30 else NULL
      )
    }
  )
  names(fits) <- names(expected_scale)

  for (correction in names(fits)) {
    raw <- fits[[correction]]$raw
    determinant <- raw$determinants[
      raw$determinants$column != "(Intercept)",
      ,
      drop = FALSE
    ]
    expected_contribution <- expected_scale[[correction]] *
      determinant$marginal_effect *
      determinant$generalized_concentration_covariance

    expect_equal(
      raw$index$correction_scale,
      expected_scale[[correction]],
      tolerance = 1e-15,
      info = correction
    )
    expect_equal(
      determinant$contribution,
      expected_contribution,
      tolerance = 1e-14,
      info = correction
    )
    expect_equal(
      raw$index$total,
      expected_scale[[correction]] *
        raw$index$generalized_concentration_covariance_y,
      tolerance = 1e-14,
      info = correction
    )
  }
})

test_that("phase 5 binary corrections match frozen Stata values", {
  data <- .phase5_wvw_data()
  expected <- c(
    standard = 0.041918945312499838,
    generalized = 0.016767578125000075,
    erreygers = 0.067070312500000298,
    wagstaff = 0.069864908854166657
  )

  for (correction in names(expected)) {
    bounded <- correction %in% c("erreygers", "wagstaff")
    fit <- .phase5_fit(
      data,
      "y_binary",
      correction,
      if (bounded) 0 else NULL,
      if (bounded) 1 else NULL
    )
    expect_equal(
      fit$raw$index$total,
      expected[[correction]],
      tolerance = 1e-12,
      info = correction
    )
  }
})

test_that("phase 5 theoretical bounds are explicit and enforced", {
  data <- .phase5_wvw_data()

  expect_error(
    .phase5_fit(data, "y_bounded", "erreygers"),
    "require explicit theoretical dep_min and dep_max",
    fixed = TRUE
  )
  expect_error(
    .phase5_fit(data, "y_bounded", "wagstaff", 10, NULL),
    "require explicit theoretical dep_min and dep_max",
    fixed = TRUE
  )
  expect_error(
    .phase5_fit(data, "y_bounded", "erreygers", 11, 30),
    "fall outside the declared theoretical bounds",
    fixed = TRUE
  )
  expect_error(
    .phase5_fit(data, "y_bounded", "wagstaff", 10, 29),
    "fall outside the declared theoretical bounds",
    fixed = TRUE
  )
  expect_error(
    .phase5_fit(data, "y_bounded", "erreygers", 10, Inf),
    "finite numeric scalars",
    fixed = TRUE
  )
  expect_error(
    .phase5_fit(data, "y_bounded", "erreygers", 30, 10),
    "dep_min < dep_max",
    fixed = TRUE
  )
  expect_error(
    .phase5_fit(
      data,
      "y_bounded",
      "erreygers",
      10,
      30,
      shortfall = NA
    ),
    "is_shortfall must be exactly TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("phase 5 complements obey the intended orientation properties", {
  data <- .phase5_wvw_data()

  for (correction in c("erreygers", "wagstaff")) {
    attainment <- .phase5_fit(
      data,
      "y_bounded",
      correction,
      10,
      30,
      shortfall = FALSE
    )
    shortfall <- .phase5_fit(
      data,
      "y_bounded_short",
      correction,
      10,
      30,
      shortfall = TRUE
    )
    expect_equal(
      unname(unlist(
        shortfall$raw$index[c("total", "explained", "residual")]
      )),
      -unname(unlist(
        attainment$raw$index[c("total", "explained", "residual")]
      )),
      tolerance = 1e-12,
      info = correction
    )
    expect_true(grepl(
      "SHORTFALL",
      shortfall$diagnostics$orientation,
      fixed = TRUE
    ))
    expect_true(grepl(
      "ATTAINMENT",
      attainment$diagnostics$orientation,
      fixed = TRUE
    ))
  }

  flag_true <- .phase5_fit(
    data, "y_bounded_short", "wagstaff", 10, 30, TRUE
  )
  flag_false <- .phase5_fit(
    data, "y_bounded_short", "wagstaff", 10, 30, FALSE
  )
  expect_equal(
    flag_true$raw$estimates$detailed,
    flag_false$raw$estimates$detailed,
    tolerance = 0
  )
  expect_equal(
    flag_true$raw$estimates$overall,
    flag_false$raw$estimates$overall,
    tolerance = 0
  )
})

test_that("phase 5 records rineq observed-range restriction", {
  skip_if_not_installed("rineq")
  data <- .phase5_wvw_data()
  fit <- .phase5_fit(
    data,
    "y_bounded_interior",
    "erreygers",
    10,
    30
  )
  standardized <- (data$y_bounded_interior - 10) / 20
  native <- rineq::ci(
    ineqvar = data$ses,
    outcome = standardized,
    weights = data$weight,
    type = "CIc",
    method = "direct",
    rank_function = rineq::rank_gwt
  )$concentration_index
  generalized <- rineq::ci(
    ineqvar = data$ses,
    outcome = data$y_bounded_interior,
    weights = data$weight,
    type = "CIg",
    method = "direct",
    rank_function = rineq::rank_gwt
  )$concentration_index
  aligned <- 4 * generalized / 20
  observed_range <- diff(range(data$y_bounded_interior))

  expect_equal(fit$raw$index$total, aligned, tolerance = 1e-12)
  expect_gt(abs(native - fit$raw$index$total), 1e-4)
  expect_equal(
    native,
    fit$raw$index$total * 20 / observed_range,
    tolerance = 1e-12
  )
})
