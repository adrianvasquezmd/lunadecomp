.phase7_wvw_data <- function() {
  id <- seq_len(320L)
  data <- data.frame(
    id = id,
    ses = rep(
      c(
        1, 2, 4, 7, 10, 14, 19, 25,
        32, 40, 49, 60, 72, 85, 99, 115
      ),
      each = 20L
    ),
    x1 = (((7L * id) %% 31L) - 15) / 8,
    x2 = (((11L * id) %% 37L) - 18) / 10,
    x3 = as.integer((id %% 5L) %in% c(0L, 1L)),
    region = 1L + ((3L * id) %% 4L),
    education = 1L + ((5L * id) %% 3L),
    weight = 1L + ((5L * id) %% 7L)
  )
  region_effect <- c(0, 0.42, -0.33, 0.58)[data$region]
  education_effect <- c(0, 0.25, -0.21)[data$education]
  eta_logit <- with(
    data,
    -0.28 + 0.56 * x1 - 0.38 * x2 + 0.61 * x3 +
      region_effect + education_effect
  )
  eta_probit <- with(
    data,
    -0.18 + 0.34 * x1 - 0.27 * x2 + 0.39 * x3 +
      0.70 * region_effect + 0.75 * education_effect
  )
  u_logit <- (((43L * id) %% 331L) + 0.5) / 331
  u_probit <- (((97L * id) %% 337L) + 0.5) / 337
  data$y_logit <- as.integer(
    u_logit < stats::plogis(eta_logit)
  )
  data$y_probit <- as.integer(
    u_probit < stats::pnorm(eta_probit)
  )
  for (level in 1:4) {
    data[[paste0("region", level)]] <-
      as.integer(data$region == level)
  }
  data
}

.phase7_wvw_fit <- function(
    model = "logit",
    correction = "standard",
    region_levels = c(1, 2, 3, 4),
    manual = FALSE,
    weight_multiplier = 1
) {
  data <- .phase7_wvw_data()
  data$weight <- weight_multiplier * data$weight
  predictors <- if (manual) {
    c("x1", "x2", "x3", "region2", "region3", "region4")
  } else {
    data$region <- factor(data$region, levels = region_levels)
    c("x1", "x2", "x3", "region")
  }
  wvw_decomp(
    data = data,
    dep_var = paste0("y_", model),
    indep_vars = predictors,
    ses_var = "ses",
    correction = correction,
    model_type = model,
    weight_var = "weight",
    vce_method = "linearized",
    quiet = TRUE
  )
}

test_that("phase 7 weighted logit matches frozen Stata margins", {
  raw <- .phase7_wvw_fit("logit")$raw
  expected_coefficients <- c(
    `(Intercept)` = -0.42084451114335975,
    x1 = 0.51835163663138895,
    x2 = -0.74405858391154289,
    x3 = 0.85777135041047448,
    region2 = 0.25240351449029180,
    region3 = -0.67435413091665863,
    region4 = 0.70126215000309167
  )
  expected_ame <- c(
    `(Intercept)` = 0,
    x1 = 0.10019081917714709,
    x2 = -0.14381711905522176,
    x3 = 0.16734030204437861,
    region2 = 0.050634893418915639,
    region3 = -0.13126440004681833,
    region4 = 0.13887221851368076
  )

  expect_equal(
    raw$model$coefficients,
    expected_coefficients,
    tolerance = 1e-8
  )
  expect_equal(
    stats::setNames(
      raw$determinants$marginal_effect,
      raw$determinants$column
    ),
    expected_ame,
    tolerance = 1e-8
  )
  expect_equal(
    unlist(raw$index[c(
      "total", "explained", "residual",
      "nonlinear_approximation_error"
    )]),
    c(
      total = 0.003357345813954624,
      explained = -0.0058484445926349944,
      residual = 0.0057566393699797751,
      nonlinear_approximation_error =
        0.0034491510366098429
    ),
    tolerance = 1e-9
  )
})

test_that("phase 7 weighted probit matches frozen Stata margins", {
  raw <- .phase7_wvw_fit("probit")$raw
  expected_coefficients <- c(
    `(Intercept)` = -0.35501874013255708,
    x1 = 0.41354628739302884,
    x2 = -0.32245419238878859,
    x3 = 0.51075549224431094,
    region2 = 0.32678328640874754,
    region3 = -0.23935792093594310,
    region4 = 0.39742332412758186
  )
  expected_ame <- c(
    `(Intercept)` = 0,
    x1 = 0.13571745870406154,
    x2 = -0.10582289062574871,
    x3 = 0.16933082865982824,
    region2 = 0.10978236311910072,
    region3 = -0.078285676594498976,
    region4 = 0.13325821337813293
  )

  expect_equal(
    raw$model$coefficients,
    expected_coefficients,
    tolerance = 1e-7
  )
  expect_equal(
    stats::setNames(
      raw$determinants$marginal_effect,
      raw$determinants$column
    ),
    expected_ame,
    tolerance = 1e-7
  )
  expect_equal(
    unlist(raw$index[c(
      "total", "explained", "residual",
      "nonlinear_approximation_error"
    )]),
    c(
      total = -0.038132087227414332,
      explained = -0.0055582433843455346,
      residual = -0.033903959687732785,
      nonlinear_approximation_error =
        0.0013301158446639844
    ),
    tolerance = 1e-7
  )
})

test_that("phase 7 continuous and binary AMEs have the stated estimands", {
  for (model in c("logit", "probit")) {
    raw <- .phase7_wvw_fit(model)$raw
    matrix_x <- raw$model$matrix
    coefficients <- raw$model$coefficients
    weights <- raw$design$normalized_weights
    linear_predictor <- as.vector(matrix_x %*% coefficients)
    density <- if (model == "logit") {
      stats::dlogis(linear_predictor)
    } else {
      stats::dnorm(linear_predictor)
    }
    observed_ame <- stats::setNames(
      raw$determinants$marginal_effect,
      raw$determinants$column
    )

    for (term in c("x1", "x2")) {
      expect_equal(
        observed_ame[[term]],
        coefficients[[term]] * sum(weights * density),
        tolerance = 1e-14,
        info = paste(model, term)
      )
    }

    matrix_one <- matrix_x
    matrix_zero <- matrix_x
    matrix_one[, "x3"] <- 1
    matrix_zero[, "x3"] <- 0
    inverse_link <- if (model == "logit") {
      stats::plogis
    } else {
      stats::pnorm
    }
    expected_x3 <- sum(weights * (
      inverse_link(as.vector(matrix_one %*% coefficients)) -
        inverse_link(as.vector(matrix_zero %*% coefficients))
    ))
    expect_equal(
      observed_ame[["x3"]],
      expected_x3,
      tolerance = 1e-14,
      info = model
    )
  }
})

test_that("phase 7 factor AMEs are changes from the omitted level", {
  for (model in c("logit", "probit")) {
    raw <- .phase7_wvw_fit(model)$raw
    matrix_x <- raw$model$matrix
    coefficients <- raw$model$coefficients
    weights <- raw$design$normalized_weights
    inverse_link <- if (model == "logit") {
      stats::plogis
    } else {
      stats::pnorm
    }
    factor_columns <- paste0("region", 2:4)

    for (column in factor_columns) {
      matrix_level <- matrix_x
      matrix_reference <- matrix_x
      matrix_level[, factor_columns] <- 0
      matrix_reference[, factor_columns] <- 0
      matrix_level[, column] <- 1
      expected <- sum(weights * (
        inverse_link(as.vector(matrix_level %*% coefficients)) -
          inverse_link(
            as.vector(matrix_reference %*% coefficients)
          )
      ))
      observed <- raw$determinants$marginal_effect[
        raw$determinants$column == column
      ]
      expect_equal(
        observed,
        expected,
        tolerance = 1e-14,
        info = paste(model, column)
      )
    }
  }
})

test_that("phase 7 factor blocks are invariant to reference level", {
  for (model in c("logit", "probit")) {
    base1 <- .phase7_wvw_fit(model)$raw
    base2 <- .phase7_wvw_fit(
      model,
      region_levels = c(2, 1, 3, 4)
    )$raw
    contribution <- function(raw) {
      sum(
        raw$determinants$contribution[
          raw$determinants$source_term == "region"
        ]
      )
    }

    expect_equal(
      base2$model$fitted,
      base1$model$fitted,
      tolerance = 1e-8,
      info = model
    )
    expect_equal(
      unlist(base2$index[c(
        "total", "explained", "residual",
        "nonlinear_approximation_error"
      )]),
      unlist(base1$index[c(
        "total", "explained", "residual",
        "nonlinear_approximation_error"
      )]),
      tolerance = 1e-8,
      info = model
    )
    expect_equal(
      contribution(base2),
      contribution(base1),
      tolerance = 1e-7,
      info = model
    )
  }
})

test_that("phase 7 one-hot dummies use independent toggle AMEs", {
  factor_fit <- .phase7_wvw_fit("logit")$raw
  manual_fit <- .phase7_wvw_fit(
    "logit",
    manual = TRUE
  )$raw

  expect_equal(
    manual_fit$model$coefficients,
    factor_fit$model$coefficients,
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  expect_equal(
    manual_fit$model$fitted,
    factor_fit$model$fitted,
    tolerance = 1e-10
  )
  expect_gt(
    max(abs(
      manual_fit$determinants$marginal_effect[-1] -
        factor_fit$determinants$marginal_effect[-1]
    )),
    1e-4
  )
  expect_gt(
    abs(
      manual_fit$index$explained -
        factor_fit$index$explained
    ),
    1e-6
  )
})

test_that("phase 7 binary corrections match frozen Stata values", {
  expected <- list(
    logit = data.frame(
      correction = c(
        "standard", "generalized", "erreygers", "wagstaff"
      ),
      total = c(
        0.003357345813954624,
        0.0016655991304432209,
        0.0066623965217728837,
        0.0066628006570598735
      ),
      explained = c(
        -0.0058484445926349944,
        -0.002901447979367984,
        -0.011605791917471936,
        -0.011606495914904688
      ),
      residual = c(
        0.0057566393699797751,
        0.0028559028650133309,
        0.011423611460053324,
        0.011424304406575011
      )
    ),
    probit = data.frame(
      correction = c(
        "standard", "generalized", "erreygers", "wagstaff"
      ),
      total = c(
        -0.038132087227414332,
        -0.018561179530478162,
        -0.074244718121912648,
        -0.074296813353566013
      ),
      explained = c(
        -0.0055582433843455346,
        -0.0027055312423800307,
        -0.010822124969520123,
        -0.010829718521243804
      ),
      residual = c(
        -0.033903959687732785,
        -0.016503095642393297,
        -0.066012382569573189,
        -0.066058701424960398
      )
    )
  )

  for (model in names(expected)) {
    for (row in seq_len(nrow(expected[[model]]))) {
      target <- expected[[model]][row, ]
      raw <- .phase7_wvw_fit(
        model,
        correction = target$correction
      )$raw
      expect_equal(
        c(
          total = raw$index$total,
          explained = raw$index$explained,
          residual = raw$index$residual
        ),
        unlist(target[c("total", "explained", "residual")]),
        tolerance = 1e-7,
        info = paste(model, target$correction)
      )
      expect_equal(
        raw$index$total,
        raw$index$explained +
          raw$index$residual +
          raw$index$nonlinear_approximation_error,
        tolerance = 1e-12,
        info = paste(model, target$correction)
      )
    }
  }
})

test_that("phase 7 nonlinear points are invariant to weight scale", {
  for (model in c("logit", "probit")) {
    base <- .phase7_wvw_fit(model)$raw
    scaled <- .phase7_wvw_fit(
      model,
      weight_multiplier = 10
    )$raw

    expect_equal(
      scaled$model$coefficients,
      base$model$coefficients,
      tolerance = 1e-8,
      info = model
    )
    expect_equal(
      scaled$determinants$marginal_effect,
      base$determinants$marginal_effect,
      tolerance = 1e-8,
      info = model
    )
    expect_equal(
      unlist(scaled$index),
      unlist(base$index),
      tolerance = 1e-8,
      info = model
    )
  }
})

test_that("phase 7 rejects invalid or perfectly separated outcomes", {
  data <- .phase7_wvw_data()
  data$region <- factor(data$region)

  nonbinary <- data
  nonbinary$y_logit[1L] <- 2
  expect_error(
    wvw_decomp(
      nonbinary,
      "y_logit",
      c("x1", "x2"),
      ses_var = "ses",
      model_type = "logit",
      quiet = TRUE
    ),
    "strictly binary",
    fixed = TRUE
  )

  constant <- data
  constant$y_logit[] <- 1
  expect_error(
    wvw_decomp(
      constant,
      "y_logit",
      c("x1", "x2"),
      ses_var = "ses",
      model_type = "logit",
      quiet = TRUE
    ),
    "has no variation",
    fixed = TRUE
  )

  separated <- data
  separated$y_logit <- as.integer(separated$x1 > 0)
  expect_error(
    wvw_decomp(
      separated,
      "y_logit",
      c("x1", "x2"),
      ses_var = "ses",
      model_type = "logit",
      quiet = TRUE
    ),
    "perfect separation detected",
    fixed = TRUE
  )
  separated$y_probit <- separated$y_logit
  expect_error(
    wvw_decomp(
      separated,
      "y_probit",
      c("x1", "x2"),
      ses_var = "ses",
      model_type = "probit",
      quiet = TRUE
    ),
    "perfect separation detected",
    fixed = TRUE
  )
})
