make_fairlie_phase1_data <- function() {
  set.seed(20260726)
  n_per_group <- 64L
  n <- 2L * n_per_group
  id <- seq_len(n)
  group <- rep(0:1, each = n_per_group)
  x1 <- as.numeric(scale(
    seq(-2.4, 2.4, length.out = n) +
      0.35 * sin(id / 4) - 0.45 * group
  ))
  x2 <- as.numeric(scale(
    cos(id / 7) + ((id * 13L) %% 31L) / 10 +
      0.30 * group
  ))
  x3 <- as.numeric(scale(
    sin(id / 9) - ((id * 7L) %% 19L) / 12 -
      0.20 * group
  ))
  eta <- -0.30 + 0.70 * x1 - 0.45 * x2 + 0.35 * x3 + 0.40 * group
  data.frame(
    id = id,
    group = group,
    y_binary = rbinom(n, size = 1L, prob = plogis(eta)),
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
}

fit_fairlie_phase1 <- function(
    data,
    model = "logit",
    reference = "group0",
    reps = 1,
    seed = 20260726
) {
  f_decomp(
    data = data,
    dep_var = "y_binary",
    group_var = "group",
    indep_vars = c("x1", "x2", "x3"),
    favored_group = 0,
    ref_method = reference,
    model_type = model,
    randomize_order = FALSE,
    vce_method = "linearized",
    reps = reps,
    seed = seed,
    quiet = TRUE
  )
}

test_that("f_decomp exposes unrounded Phase 1 audit results", {
  fit <- fit_fairlie_phase1(make_fairlie_phase1_data())

  expect_true(all(c(
    "estimates", "standard_errors", "vcov",
    "overall_estimates", "detailed_estimates",
    "reference_coefficients", "reference_vcov",
    "reference_predictions", "analytic_sample",
    "decomposition_blocks", "matching"
  ) %in% names(fit$raw)))

  public_overall <- setNames(
    fit$results_overall$Estimate,
    fit$results_overall$Term
  )
  expect_equal(
    public_overall,
    fit$raw$overall_estimates,
    tolerance = 0
  )

  public_detail <- setNames(
    fit$results_detailed_explained$Estimate,
    sub(
      " \\(Continuous\\)$",
      "",
      fit$results_detailed_explained$Term
    )
  )
  raw_detail <- fit$raw$detailed_estimates
  names(raw_detail) <- sub("^Exp_", "", names(raw_detail))
  expect_equal(
    public_detail[names(raw_detail)],
    raw_detail,
    tolerance = 0
  )
  expect_false(
    identical(
      fit$raw$overall_estimates[["Explained"]],
      round(fit$raw$overall_estimates[["Explained"]], 6)
    )
  )
})

test_that("f_decomp estimates only models required by the reference", {
  data <- make_fairlie_phase1_data()

  fit_group0 <- fit_fairlie_phase1(data, reference = "group0")
  expect_identical(fit_group0$raw$fitted_models, "group0")
  expect_null(fit_group0$raw$group_coefficients$group_1)
  expect_setequal(
    names(fit_group0$models_coefficients),
    c("Group_0", "Reference")
  )

  fit_group1 <- fit_fairlie_phase1(data, reference = "group1")
  expect_identical(fit_group1$raw$fitted_models, "group1")
  expect_null(fit_group1$raw$group_coefficients$group_0)
  expect_setequal(
    names(fit_group1$models_coefficients),
    c("Group_1", "Reference")
  )

  fit_neumark <- fit_fairlie_phase1(data, reference = "neumark")
  expect_identical(fit_neumark$raw$fitted_models, "pooled")
  expect_null(fit_neumark$raw$group_coefficients$group_0)
  expect_null(fit_neumark$raw$group_coefficients$group_1)
  expect_identical(names(fit_neumark$models_coefficients), "Reference")

  fit_reimers <- fit_fairlie_phase1(data, reference = "reimers")
  expect_setequal(fit_reimers$raw$fitted_models, c("group0", "group1"))
  expect_false(is.null(fit_reimers$raw$group_coefficients$group_0))
  expect_false(is.null(fit_reimers$raw$group_coefficients$group_1))
})

test_that("unused group separation does not invalidate a group0 reference", {
  data <- make_fairlie_phase1_data()
  data$y_binary[data$group == 1] <- 1

  fit <- expect_no_error(
    fit_fairlie_phase1(data, reference = "group0")
  )
  expect_identical(fit$raw$fitted_models, "group0")
  expect_null(fit$raw$group_coefficients$group_1)
})

test_that("Fairlie details are no longer rescaled", {
  data <- make_fairlie_phase1_data()
  data <- data[-which(data$group == 1)[seq_len(11)], ]

  fit <- fit_fairlie_phase1(
    data,
    reference = "group0",
    reps = 7,
    seed = 1827
  )

  expect_false(fit$raw$matching$rescaled)
  expect_equal(
    sum(fit$raw$detailed_estimates),
    fit$raw$matching$detailed_sum,
    tolerance = 0
  )
  expect_equal(
    fit$raw$matching$difference,
    fit$raw$matching$detailed_sum -
      fit$raw$matching$full_sample_explained,
    tolerance = 0
  )
  expect_gt(abs(fit$raw$matching$difference), 1e-6)
})

test_that("Phase 1 deterministic point estimates match frozen Stata values", {
  data <- make_fairlie_phase1_data()
  expected <- list(
    logit_group0 = c(
      Explained = -0.2812006990439595,
      Exp_x1 = -0.3084654466939129,
      Exp_x2 = 0.0025891651793874936,
      Exp_x3 = 0.024675582470565945
    ),
    logit_group1 = c(
      Explained = -0.16114237428726341,
      Exp_x1 = -0.18242991083633614,
      Exp_x2 = 0.0059920921812124618,
      Exp_x3 = 0.015295444367860308
    ),
    probit_group0 = c(
      Explained = -0.29013934013107301,
      Exp_x1 = -0.32000847871840599,
      Exp_x2 = 0.0023872665350154051,
      Exp_x3 = 0.027481872052317603
    ),
    probit_group1 = c(
      Explained = -0.15923139170019029,
      Exp_x1 = -0.18135372962565227,
      Exp_x2 = 0.0058523478945496678,
      Exp_x3 = 0.01626999003091233
    )
  )

  for (scenario in names(expected)) {
    model <- sub("_group[01]$", "", scenario)
    reference <- sub("^[^_]+_", "", scenario)
    fit <- fit_fairlie_phase1(
      data,
      model = model,
      reference = reference
    )
    actual <- c(
      Explained = fit$raw$overall_estimates[["Explained"]],
      fit$raw$detailed_estimates
    )
    expect_equal(
      actual,
      expected[[scenario]],
      tolerance = if (model == "logit") 1e-10 else 5e-8,
      info = scenario
    )
  }
})
