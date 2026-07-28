make_fairlie_phase7_data <- function() {
  n_per_group <- 60L
  id <- seq_len(2L * n_per_group)
  group <- rep(0:1, each = n_per_group)
  x1 <- as.numeric(scale(
    seq(-2.5, 2.5, length.out = length(id)) +
      0.35 * sin(id / 5) - 0.4 * group
  ))
  x2 <- as.integer(((id * 7L + group) %% 13L) < 6L)
  x3 <- as.integer(((id * 5L + 2L * group) %% 17L) < 7L)
  y <- as.integer(
    ((id * 11L + group * 3L + x2 * 2L) %% 19L) <
      (8L + group + x3)
  )
  data.frame(
    id = id,
    group = group,
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    weight = 0.5 + (id %% 17L) / 5
  )
}

fit_fairlie_phase7 <- function(
    data,
    model = "logit",
    reference = "group0",
    variables = c("x1", "x2", "x3"),
    group_levels = NULL,
    weight_var = NULL
) {
  f_decomp(
    data = data,
    dep_var = "y",
    group_var = "group",
    group_levels = group_levels,
    indep_vars = variables,
    ref_method = reference,
    model_type = model,
    randomize_order = FALSE,
    weight_var = weight_var,
    vce_method = "linearized",
    reps = 1L,
    quiet = TRUE
  )
}

test_that("Phase 7 exposes an exact row-level analytic-sample audit", {
  data <- make_fairlie_phase7_data()
  data$group <- ifelse(data$group == 0, 10, 40)
  data$y[c(2L, 62L)] <- NA
  data$group[c(3L, 63L)] <- c(20, 50)
  data$weight[c(4L, 64L)] <- 0

  fit <- fit_fairlie_phase7(
    data,
    group_levels = c(10, 40),
    weight_var = "weight"
  )

  expect_equal(fit$summary_stats$N, 114)
  expect_equal(fit$summary_stats$N0, 57)
  expect_equal(fit$summary_stats$N1, 57)
  expect_identical(fit$raw$sample_flow$missing_rows, c(2L, 62L))
  expect_identical(
    fit$raw$sample_flow$group_filtered_rows,
    c(3L, 63L)
  )
  expect_identical(
    fit$raw$sample_flow$zero_weight_rows,
    c(4L, 64L)
  )
  expect_identical(
    fit$raw$sample_flow$analytic_rows,
    fit$raw$analytic_sample
  )
  expect_false(any(2:4 %in% fit$raw$analytic_sample))
  expect_false(any(62:64 %in% fit$raw$analytic_sample))
})

test_that("Phase 7 requires outcome variation only in fitted equations", {
  data <- make_fairlie_phase7_data()
  required_constant <- data
  required_constant$y[required_constant$group == 0] <- 0
  expect_error(
    fit_fairlie_phase7(required_constant, reference = "group0"),
    "group 0 model cannot be estimated"
  )

  unused_constant <- data
  unused_constant$y[unused_constant$group == 1] <- 0
  fit <- fit_fairlie_phase7(unused_constant, reference = "group0")
  expect_equal(fit$raw$fitted_models, "group0")
  expect_true(all(is.finite(fit$raw$reference_coefficients)))
})

test_that("Phase 7 rejects perfect separation but accepts sparse events", {
  data <- make_fairlie_phase7_data()
  separated <- data
  separated$separator <- separated$id %% 2L
  separated$y[separated$group == 0] <-
    separated$separator[separated$group == 0]
  for (model in c("logit", "probit")) {
    expect_error(
      fit_fairlie_phase7(
        separated,
        model = model,
        variables = c("x1", "x2", "x3", "separator")
      ),
      "perfect separation"
    )
  }

  sparse <- data
  sparse$y[sparse$group == 0] <- 0
  sparse$y[sparse$id %in% c(1L, 31L)] <- 1
  expect_type(fit_fairlie_phase7(sparse), "list")
  expect_type(
    fit_fairlie_phase7(sparse, model = "probit"),
    "list"
  )
})

test_that("Phase 7 follows Stata omission behavior for aliased columns", {
  data <- make_fairlie_phase7_data()
  data$x_duplicate <- 2 * data$x1
  variables <- c("x1", "x2", "x3", "x_duplicate")

  for (model in c("logit", "probit")) {
    fit <- fit_fairlie_phase7(
      data,
      model = model,
      variables = variables
    )
    expect_match(fit$diagnostics$vif, "exact collinearity")
    expect_match(fit$diagnostics$aliased, "aliased")
    expect_equal(fit$raw$reference_coefficients[["x_duplicate"]], 0)
  }

  data$within_group_constant <- ifelse(data$group == 0, 0, data$x1)
  fit_constant <- fit_fairlie_phase7(
    data,
    variables = c("x1", "x2", "x3", "within_group_constant")
  )
  expect_match(
    fit_constant$diagnostics$zero_variance,
    "zero within-group variance"
  )
  expect_equal(
    unname(fit_constant$raw$reference_coefficients[
      "within_group_constant"
    ]),
    0
  )
})

test_that("Phase 7 preserves strict and explicit group-selection contracts", {
  data <- make_fairlie_phase7_data()
  invalid <- data
  invalid$group[c(4L, 64L)] <- 2
  expect_error(
    fit_fairlie_phase7(invalid),
    "exactly two distinct groups"
  )

  selected <- data
  selected$group <- ifelse(selected$group == 0, 10, 40)
  selected$group[c(4L, 64L)] <- c(20, 50)
  fit <- fit_fairlie_phase7(
    selected,
    group_levels = c(10, 40)
  )
  expect_equal(fit$summary_stats$n_group_filtered, 2)
  expect_identical(
    fit$summary_stats$selected_group_levels,
    c("10", "40")
  )
})
