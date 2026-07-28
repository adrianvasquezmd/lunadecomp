phase2_ke_data <- function() {
  n <- 36L
  id <- seq_len(n)
  ses <- 25 + ((17 * id) %% 53) + id / 1000
  income <- 300 + 11 * id + 20 * sin(id / 5)
  x1 <- ((7 * id) %% 19) - 9 + id / 13
  x2 <- cos(id / 4) + 0.15 * (id %% 3)
  ses_scaled <- (ses - min(ses)) / (max(ses) - min(ses))
  health <- 0.28 +
    0.009 * x1 -
    0.035 * x2 +
    0.16 * ses_scaled +
    0.018 * sin(id / 3)
  order <- c(seq(2, n, by = 2), seq(1, n, by = 2))

  data.frame(
    id = id[order],
    health = health[order],
    ses = ses[order],
    income = income[order],
    x1 = x1[order],
    x2 = x2[order],
    stringsAsFactors = FALSE
  )
}

manual_phase2_ke <- function(data, index_type) {
  n <- nrow(data)
  health_mean <- mean(data$health)
  if (index_type == "rank") {
    position <- (rank(data$ses) - 0.5) / n
    target <- 2 * position * data$health - health_mean
  } else {
    position <- data$income / mean(data$income)
    target <- position * data$health - health_mean
  }
  matrix <- cbind("(Intercept)" = 1, x1 = data$x1, x2 = data$x2)
  model <- stats::lm.fit(matrix, target)
  coefficients <- model$coefficients
  names(coefficients) <- colnames(matrix)
  fitted <- as.vector(matrix %*% coefficients)
  residual <- target - fitted

  list(
    position = position,
    target = target,
    coefficients = coefficients,
    fitted = fitted,
    residual = residual,
    index = mean(target),
    r_squared =
      1 - sum(residual^2) / sum((target - mean(target))^2)
  )
}

test_that("KE phase 2 rank core equals its independent equations", {
  data <- phase2_ke_data()
  manual <- manual_phase2_ke(data, "rank")
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))

  order <- order(data$ses)
  expect_equal(
    fit$raw$target$position,
    manual$position[order],
    tolerance = 2e-14
  )
  expect_equal(
    fit$raw$target$target,
    manual$target[order],
    tolerance = 2e-14
  )
  expect_equal(
    fit$raw$model$coefficients,
    manual$coefficients,
    tolerance = 1e-14
  )
  expect_equal(
    fit$raw$target$fitted,
    manual$fitted[order],
    tolerance = 1e-14
  )
  expect_equal(
    fit$raw$target$residual,
    manual$residual[order],
    tolerance = 1e-14
  )
  expect_lt(
    abs(fit$results_overall$Estimate - manual$index),
    5e-14
  )
  expect_lt(
    abs(fit$model_metrics$R_Squared - manual$r_squared),
    5e-14
  )
})

test_that("KE phase 2 level core equals its independent equations", {
  data <- phase2_ke_data()
  manual <- manual_phase2_ke(data, "level")
  fit <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "income",
    index_type = "level",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expect_equal(
    fit$raw$target$position,
    manual$position,
    tolerance = 2e-14
  )
  expect_equal(
    fit$raw$target$target,
    manual$target,
    tolerance = 2e-14
  )
  expect_equal(
    fit$raw$model$coefficients,
    manual$coefficients,
    tolerance = 1e-14
  )
  expect_equal(
    fit$raw$target$fitted,
    manual$fitted,
    tolerance = 1e-14
  )
  expect_equal(
    fit$raw$target$residual,
    manual$residual,
    tolerance = 1e-14
  )
  expect_lt(
    abs(fit$results_overall$Estimate - manual$index),
    5e-14
  )
  expect_lt(
    abs(fit$model_metrics$R_Squared - manual$r_squared),
    5e-14
  )
})

test_that("KE phase 2 rank core is invariant to order and affine SES", {
  data <- phase2_ke_data()
  base <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))
  reversed <- suppressWarnings(ke_decomp(
    data = data[nrow(data):1, ],
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))
  affine_data <- data
  affine_data$ses <- 7 + 3 * affine_data$ses
  affine <- suppressWarnings(ke_decomp(
    data = affine_data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))

  for (variant in list(reversed, affine)) {
    expect_equal(
      variant$results_overall$Estimate,
      base$results_overall$Estimate,
      tolerance = 2e-14
    )
    expect_equal(
      variant$raw$model$coefficients,
      base$raw$model$coefficients,
      tolerance = 1e-14
    )
    expect_equal(
      variant$raw$target$position,
      base$raw$target$position,
      tolerance = 2e-14
    )
    expect_equal(
      variant$raw$target$target,
      base$raw$target$target,
      tolerance = 2e-14
    )
  }
})

test_that("KE phase 2 level core is invariant to order and SES scale", {
  data <- phase2_ke_data()
  base <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "income",
    index_type = "level",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))
  reversed <- suppressWarnings(ke_decomp(
    data = data[nrow(data):1, ],
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "income",
    index_type = "level",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))
  scaled_data <- data
  scaled_data$income <- 10 * scaled_data$income
  scaled <- suppressWarnings(ke_decomp(
    data = scaled_data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "income",
    index_type = "level",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expect_lt(
    abs(
      reversed$results_overall$Estimate -
        base$results_overall$Estimate
    ),
    5e-14
  )
  expect_equal(
    reversed$raw$model$coefficients,
    base$raw$model$coefficients,
    tolerance = 1e-14
  )
  expect_lt(
    abs(
      scaled$results_overall$Estimate -
        base$results_overall$Estimate
    ),
    5e-14
  )
  expect_equal(
    scaled$raw$model$coefficients,
    base$raw$model$coefficients,
    tolerance = 1e-14
  )
  expect_equal(
    scaled$raw$target$position,
    base$raw$target$position,
    tolerance = 2e-14
  )
  expect_equal(
    scaled$raw$target$target,
    base$raw$target$target,
    tolerance = 2e-14
  )
})

test_that("KE phase 2 predictor rescaling preserves fitted decomposition", {
  data <- phase2_ke_data()
  base <- suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))
  scaled_data <- data
  scaled_data$x1 <- 100 * scaled_data$x1
  scaled <- suppressWarnings(ke_decomp(
    data = scaled_data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = "rank",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  ))

  expect_lt(
    abs(
      scaled$results_overall$Estimate -
        base$results_overall$Estimate
    ),
    5e-14
  )
  expect_equal(
    100 * scaled$raw$model$coefficients[["x1"]],
    base$raw$model$coefficients[["x1"]],
    tolerance = 1e-14
  )
  expect_equal(
    scaled$raw$target$fitted,
    base$raw$target$fitted,
    tolerance = 1e-14
  )
  expect_equal(
    scaled$raw$target$residual,
    base$raw$target$residual,
    tolerance = 1e-14
  )
  expect_lt(
    abs(
      scaled$model_metrics$R_Squared -
        base$model_metrics$R_Squared
    ),
    5e-14
  )
})
