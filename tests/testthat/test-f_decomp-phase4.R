make_fairlie_phase4_data <- function() {
  set.seed(20260728)
  n_per_group <- 80L
  n <- 2L * n_per_group
  id <- seq_len(n)
  group <- rep(0:1, each = n_per_group)
  x_cont <- as.numeric(scale(
    seq(-2.5, 2.5, length.out = n) +
      0.50 * sin(id / 6) - 0.42 * group
  ))
  binary <- as.integer(
    ((id * 7L + group * 3L) %% 11L) < (4L + group)
  )
  education_code <- (id * 5L + group * 2L) %% 9L
  education <- factor(
    ifelse(
      education_code < 3,
      "A",
      ifelse(education_code < 6, "B", "C")
    ),
    levels = c("A", "B", "C")
  )
  edu_b <- as.integer(education == "B")
  edu_c <- as.integer(education == "C")
  eta <- -0.35 + 0.58 * x_cont + 0.37 * binary +
    0.28 * edu_b + 0.63 * edu_c + 0.45 * group
  data.frame(
    id = id,
    group = group,
    y_binary = rbinom(n, size = 1L, prob = plogis(eta)),
    x_cont = x_cont,
    binary = binary,
    education = education,
    edu_b = edu_b,
    edu_c = edu_c
  )
}

phase4_grouping <- function(order = c(
    "x_cont",
    "binary",
    "education"
)) {
  blocks <- list(
    x_cont = "x_cont",
    binary = "binary",
    education = "education"
  )
  blocks[order]
}

fit_fairlie_phase4 <- function(
    data,
    groupings = phase4_grouping(),
    randomize_order = FALSE,
    reps = 1,
    seed = 410101,
    indep_vars = c("x_cont", "binary", "education")
) {
  f_decomp(
    data = data,
    dep_var = "y_binary",
    group_var = "group",
    group_levels = c(0, 1),
    indep_vars = indep_vars,
    groupings = groupings,
    ref_method = "group0",
    model_type = "logit",
    randomize_order = randomize_order,
    vce_method = "linearized",
    reps = reps,
    seed = seed,
    quiet = TRUE
  )
}

canonical_phase4_details <- function(fit) {
  values <- fit$raw$detailed_estimates
  names(values) <- sub("^Exp_", "", names(values))
  values[c("x_cont", "binary", "education")]
}

test_that("a multilevel factor is exchanged as one complete block", {
  data <- make_fairlie_phase4_data()
  fit <- fit_fairlie_phase4(data)

  expect_identical(
    names(fit$raw$decomposition_blocks),
    c("x_cont", "binary", "education")
  )
  expect_identical(
    unname(fit$raw$decomposition_blocks$education),
    c("educationB", "educationC")
  )
  expect_equal(
    sum(fit$raw$detailed_estimates),
    fit$raw$overall_estimates[["Explained"]],
    tolerance = 1e-14
  )
})

test_that("automatic factor coding equals grouped manual dummies", {
  data <- make_fairlie_phase4_data()
  factor_fit <- fit_fairlie_phase4(data)
  dummy_fit <- fit_fairlie_phase4(
    data,
    indep_vars = c("x_cont", "binary", "edu_b", "edu_c"),
    groupings = list(
      x_cont = "x_cont",
      binary = "binary",
      education = c("edu_b", "edu_c")
    )
  )

  expect_equal(
    factor_fit$raw$overall_estimates,
    dummy_fit$raw$overall_estimates,
    tolerance = 0
  )
  expect_equal(
    canonical_phase4_details(factor_fit),
    canonical_phase4_details(dummy_fit),
    tolerance = 0
  )
  factor_beta <- factor_fit$raw$reference_coefficients
  names(factor_beta)[names(factor_beta) == "educationB"] <- "edu_b"
  names(factor_beta)[names(factor_beta) == "educationC"] <- "edu_c"
  expect_equal(
    factor_beta[names(dummy_fit$raw$reference_coefficients)],
    dummy_fit$raw$reference_coefficients,
    tolerance = 0
  )
})

test_that("alternative fixed paths telescope but can redistribute details", {
  data <- make_fairlie_phase4_data()
  forward <- fit_fairlie_phase4(data)
  reverse <- fit_fairlie_phase4(
    data,
    groupings = phase4_grouping(
      c("education", "binary", "x_cont")
    )
  )

  expect_equal(
    forward$raw$overall_estimates[["Explained"]],
    reverse$raw$overall_estimates[["Explained"]],
    tolerance = 0
  )
  expect_equal(
    sum(forward$raw$detailed_estimates),
    forward$raw$overall_estimates[["Explained"]],
    tolerance = 1e-14
  )
  expect_equal(
    sum(reverse$raw$detailed_estimates),
    reverse$raw$overall_estimates[["Explained"]],
    tolerance = 1e-14
  )
  expect_gt(
    max(abs(
      canonical_phase4_details(forward) -
        canonical_phase4_details(reverse)
    )),
    1e-4
  )
})

test_that("random ordering converges to the exact permutation average", {
  data <- make_fairlie_phase4_data()
  permutations <- list(
    c("x_cont", "binary", "education"),
    c("x_cont", "education", "binary"),
    c("binary", "x_cont", "education"),
    c("binary", "education", "x_cont"),
    c("education", "x_cont", "binary"),
    c("education", "binary", "x_cont")
  )
  fixed_values <- vapply(
    permutations,
    function(order) {
      canonical_phase4_details(fit_fairlie_phase4(
        data,
        groupings = phase4_grouping(order)
      ))
    },
    numeric(3)
  )
  exact_average <- rowMeans(fixed_values)

  randomized <- fit_fairlie_phase4(
    data,
    randomize_order = TRUE,
    reps = 6000,
    seed = 412121
  )
  observed <- canonical_phase4_details(randomized)
  monte_carlo_se <- randomized$raw$matching$detailed_monte_carlo_se
  names(monte_carlo_se) <- sub("^Exp_", "", names(monte_carlo_se))
  monte_carlo_se <- monte_carlo_se[names(observed)]

  expect_true(all(
    abs(observed - exact_average) <=
      4.5 * monte_carlo_se + 1e-12
  ))
  expect_true(all(monte_carlo_se > 0))
})

test_that("invalid grouping specifications fail before decomposition", {
  data <- make_fairlie_phase4_data()

  expect_error(
    fit_fairlie_phase4(
      data,
      groupings = list(A = "x_cont", B = "x_cont")
    ),
    "cannot belong to more than one block"
  )
  expect_error(
    fit_fairlie_phase4(
      data,
      groupings = list(A = "not_a_model_term")
    ),
    "did not match"
  )
  expect_error(
    fit_fairlie_phase4(
      data,
      groupings = list(
        Education_B = "educationB",
        Education_C = "educationC",
        Other = c("x_cont", "binary")
      )
    ),
    "multi-column factor term"
  )
  expect_error(
    fit_fairlie_phase4(
      data,
      groupings = list("x_cont", "binary")
    ),
    "unique, non-empty name"
  )
})
