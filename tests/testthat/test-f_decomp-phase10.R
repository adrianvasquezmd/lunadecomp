make_fairlie_phase10_data <- function() {
  set.seed(20260727)
  n_per_group <- 48L
  id <- seq_len(2L * n_per_group)
  group <- rep(0:1, each = n_per_group)
  x1 <- as.numeric(scale(
    seq(-2.6, 2.6, length.out = length(id)) +
      0.42 * sin(id / 5) - 0.38 * group
  ))
  x2 <- as.integer(((id * 7L + 3L * group) %% 17L) < (7L + group))
  x3 <- as.numeric(scale(
    cos(id / 8) + ((id * 11L) %% 23L) / 9 + 0.25 * group
  ))
  eta <- -0.45 + 0.68 * x1 - 0.37 * x2 + 0.31 * x3 +
    0.42 * group
  weight <- 0.65 + ((id * 13L) %% 29L) / 20 + 0.22 * group
  data.frame(
    id = id,
    group = group,
    y = stats::rbinom(length(id), 1, stats::plogis(eta)),
    x1 = x1,
    x2 = x2,
    x3 = x3,
    x12 = x1 * x2,
    weight = weight
  )
}

fit_fairlie_phase10 <- function(
    data,
    reference,
    model = "logit",
    method = "linearized",
    weight_var = NULL,
    group_levels = c(0, 1),
    indep_vars = c("x1", "x2", "x3"),
    groupings = list(),
    randomize_order = FALSE,
    reps = 1L,
    seed = 101010L
) {
  f_decomp(
    data = data,
    dep_var = "y",
    group_var = "group",
    group_levels = group_levels,
    indep_vars = indep_vars,
    groupings = groupings,
    ref_method = reference,
    model_type = model,
    randomize_order = randomize_order,
    weight_var = weight_var,
    vce_method = method,
    reps = reps,
    boot_reps = 20L,
    seed = seed,
    quiet = TRUE
  )
}

test_that("Reimers satisfies its coefficient and covariance identities", {
  data <- make_fairlie_phase10_data()

  for (model in c("logit", "probit")) {
    fit <- fit_fairlie_phase10(data, "reimers", model = model)
    b0 <- fit$raw$group_coefficients$group_0
    b1 <- fit$raw$group_coefficients$group_1
    V0 <- fit$raw$group_vcov$group_0
    V1 <- fit$raw$group_vcov$group_1

    expect_equal(
      fit$raw$reference_coefficients,
      0.5 * b0 + 0.5 * b1,
      tolerance = 0
    )
    expect_equal(
      fit$raw$reference_vcov,
      0.25 * V0 + 0.25 * V1,
      tolerance = 0
    )
    expect_equal(
      fit$raw$reference_structure$group_0_coefficient_weight,
      0.5
    )
    expect_true(
      fit$raw$reference_structure$
        weight_fixed_across_sampling_replicates
    )
    expect_false(
      fit$raw$reference_structure$fairlie_native_comparator
    )
  }
})

test_that("Cotton uses the full-sample weighted group share", {
  data <- make_fairlie_phase10_data()
  fit <- fit_fairlie_phase10(
    data,
    "cotton",
    model = "logit",
    weight_var = "weight"
  )
  p0 <- with(data, sum(weight[group == 0]) / sum(weight))
  b0 <- fit$raw$group_coefficients$group_0
  b1 <- fit$raw$group_coefficients$group_1
  V0 <- fit$raw$group_vcov$group_0
  V1 <- fit$raw$group_vcov$group_1

  expect_equal(
    fit$raw$reference_coefficients,
    p0 * b0 + (1 - p0) * b1,
    tolerance = 1e-15
  )
  expect_equal(
    fit$raw$reference_vcov,
    p0^2 * V0 + (1 - p0)^2 * V1,
    tolerance = 1e-15
  )
  expect_equal(
    fit$raw$reference_structure$group_0_coefficient_weight,
    p0,
    tolerance = 0
  )
  expect_identical(
    fit$raw$reference_structure$reference_weight_inference,
    "conditional on the full-sample reference weight"
  )
})

test_that("Reimers and Cotton coincide when group weights are equal", {
  data <- make_fairlie_phase10_data()
  reimers <- fit_fairlie_phase10(data, "reimers")
  cotton <- fit_fairlie_phase10(data, "cotton")

  expect_equal(
    reimers$raw$reference_coefficients,
    cotton$raw$reference_coefficients,
    tolerance = 0
  )
  expect_equal(
    reimers$raw$reference_vcov,
    cotton$raw$reference_vcov,
    tolerance = 0
  )
  expect_equal(
    reimers$raw$overall_estimates,
    cotton$raw$overall_estimates,
    tolerance = 0
  )
  expect_equal(
    reimers$raw$detailed_estimates,
    cotton$raw$detailed_estimates,
    tolerance = 0
  )
})

test_that("Reimers and Cotton are orientation symmetric", {
  data <- make_fairlie_phase10_data()

  for (reference in c("reimers", "cotton")) {
    forward <- fit_fairlie_phase10(data, reference)
    reversed <- fit_fairlie_phase10(
      data,
      reference,
      group_levels = c(1, 0)
    )

    expect_equal(
      forward$raw$reference_coefficients,
      reversed$raw$reference_coefficients,
      tolerance = 1e-14
    )
    expect_equal(
      forward$raw$overall_estimates[["Group_0"]],
      reversed$raw$overall_estimates[["Group_1"]],
      tolerance = 1e-14
    )
    expect_equal(
      forward$raw$overall_estimates[["Group_1"]],
      reversed$raw$overall_estimates[["Group_0"]],
      tolerance = 1e-14
    )
    expect_equal(
      forward$raw$overall_estimates[
        c("Difference", "Explained", "Unexplained")
      ],
      -reversed$raw$overall_estimates[
        c("Difference", "Explained", "Unexplained")
      ],
      tolerance = 1e-14
    )
  }
})

test_that("Cotton reference weight stays fixed inside jackknife replicates", {
  data <- make_fairlie_phase10_data()
  fit <- fit_fairlie_phase10(
    data,
    "cotton",
    method = "jackknife"
  )
  deleted_first <- data[-1, , drop = FALSE]
  family <- stats::quasibinomial(link = "logit")
  model0 <- stats::glm(
    y ~ x1 + x2 + x3,
    data = deleted_first[deleted_first$group == 0, ],
    family = family,
    control = stats::glm.control(epsilon = 1e-12, maxit = 100)
  )
  model1 <- stats::glm(
    y ~ x1 + x2 + x3,
    data = deleted_first[deleted_first$group == 1, ],
    family = family,
    control = stats::glm.control(epsilon = 1e-12, maxit = 100)
  )
  X <- stats::model.matrix(~x1 + x2 + x3, data = deleted_first)
  fixed_weight <- 0.5
  beta_fixed <- fixed_weight * stats::coef(model0) +
    (1 - fixed_weight) * stats::coef(model1)
  probability_fixed <- stats::plogis(as.vector(X %*% beta_fixed))
  expected_fixed <- mean(
    probability_fixed[deleted_first$group == 0]
  ) - mean(probability_fixed[deleted_first$group == 1])

  replicate_weight <- mean(deleted_first$group == 0)
  beta_reestimated <- replicate_weight * stats::coef(model0) +
    (1 - replicate_weight) * stats::coef(model1)
  probability_reestimated <- stats::plogis(
    as.vector(X %*% beta_reestimated)
  )
  expected_reestimated <- mean(
    probability_reestimated[deleted_first$group == 0]
  ) - mean(probability_reestimated[deleted_first$group == 1])

  observed <- unname(
    fit$raw$replication$replicates[1, "Explained"]
  )
  expect_equal(observed, expected_fixed, tolerance = 1e-14)
  expect_gt(abs(observed - expected_reestimated), 1e-8)
})

test_that("Cotton is invariant to a global weight rescaling", {
  data <- make_fairlie_phase10_data()
  scaled <- data
  scaled$weight <- 1000 * scaled$weight
  original <- fit_fairlie_phase10(
    data,
    "cotton",
    model = "probit",
    weight_var = "weight",
    reps = 25L
  )
  rescaled <- fit_fairlie_phase10(
    scaled,
    "cotton",
    model = "probit",
    weight_var = "weight",
    reps = 25L
  )

  expect_equal(
    original$raw$reference_coefficients,
    rescaled$raw$reference_coefficients,
    tolerance = 1e-9
  )
  expect_equal(
    original$raw$overall_estimates,
    rescaled$raw$overall_estimates,
    tolerance = 1e-9
  )
  expect_equal(
    original$raw$detailed_estimates,
    rescaled$raw$detailed_estimates,
    tolerance = 1e-9
  )
})

test_that("deterministic estimates are invariant to analytic-row order", {
  data <- make_fairlie_phase10_data()
  reordered <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  original <- fit_fairlie_phase10(data, "group0")
  reversed_rows <- fit_fairlie_phase10(reordered, "group0")

  expect_equal(
    original$raw$overall_estimates,
    reversed_rows$raw$overall_estimates,
    tolerance = 1e-14
  )
  expect_equal(
    original$raw$detailed_estimates,
    reversed_rows$raw$detailed_estimates,
    tolerance = 1e-14
  )
})

test_that("seeded stochastic paths are reproducible and restore RNG state", {
  data <- make_fairlie_phase10_data()
  set.seed(440011)
  state_before <- .Random.seed
  first <- fit_fairlie_phase10(
    data,
    "reimers",
    randomize_order = TRUE,
    reps = 80L,
    seed = 440012L
  )
  state_after <- .Random.seed
  second <- fit_fairlie_phase10(
    data,
    "reimers",
    randomize_order = TRUE,
    reps = 80L,
    seed = 440012L
  )

  expect_identical(state_after, state_before)
  expect_identical(.Random.seed, state_before)
  expect_equal(
    first$raw$detailed_estimates,
    second$raw$detailed_estimates,
    tolerance = 0
  )
  expect_equal(
    first$raw$matching$detailed_monte_carlo_se,
    second$raw$matching$detailed_monte_carlo_se,
    tolerance = 0
  )
})

test_that("interaction contract rejects formulas and permits one joint block", {
  data <- make_fairlie_phase10_data()

  expect_error(
    fit_fairlie_phase10(
      data,
      "group0",
      indep_vars = c("x1", "x2", "x1:x2")
    ),
    "Interaction terms detected"
  )

  joint <- fit_fairlie_phase10(
    data,
    "group0",
    indep_vars = c("x1", "x2", "x12"),
    groupings = list(interacted_system = c("x1", "x2", "x12"))
  )
  expect_identical(
    joint$raw$interaction_handling$formula_interactions,
    "rejected"
  )
  expect_false(
    joint$raw$interaction_handling$separate_hierarchical_attribution
  )
  expect_identical(
    names(joint$raw$decomposition_blocks),
    "interacted_system"
  )
  expect_equal(
    sum(joint$raw$detailed_estimates),
    joint$raw$overall_estimates[["Explained"]],
    tolerance = 1e-14
  )
})

test_that("overlap diagnostic is descriptive and never called propensity", {
  data <- make_fairlie_phase10_data()
  fit <- fit_fairlie_phase10(data, "group0")
  overlap <- fit$raw$predicted_outcome_overlap

  expect_identical(
    overlap$quantity,
    "reference-model predicted outcome probability"
  )
  expect_false(overlap$causal_propensity_score)
  expect_false(overlap$trimming_applied)
  expect_true(all(overlap$off_overlap_unweighted >= 0))
  expect_true(all(overlap$off_overlap_unweighted <= 1))
  expect_false(any(grepl(
    "propensity",
    unlist(fit$diagnostics),
    ignore.case = TRUE
  )))
})
