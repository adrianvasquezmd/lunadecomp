make_fairlie_phase2_data <- function() {
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
  eta <- -0.30 + 0.70 * x1 - 0.45 * x2 + 0.35 * x3 +
    0.40 * group
  base <- data.frame(
    id = id,
    group = group,
    y_binary = rbinom(n, size = 1L, prob = plogis(eta)),
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
  base$group_label <- ifelse(
    base$group == 0,
    "advantaged",
    "disadvantaged"
  )
  middle <- base[seq_len(16), , drop = FALSE]
  middle$id <- max(base$id) + seq_len(nrow(middle))
  middle$group <- 2
  middle$group_label <- "middle"
  middle$y_binary <- as.integer(seq_len(nrow(middle)) %% 3 == 0)
  middle$x1 <- middle$x1 + 0.17
  middle$x2 <- middle$x2 - 0.11
  middle$x3 <- middle$x3 + 0.08
  rbind(base, middle)
}

fit_fairlie_phase2 <- function(
    data,
    model = "logit",
    levels = c("advantaged", "disadvantaged"),
    reference = "pooled",
    anchor = "favored"
) {
  f_decomp(
    data = data,
    dep_var = "y_binary",
    group_var = "group_label",
    group_levels = levels,
    indep_vars = c("x1", "x2", "x3"),
    ref_method = reference,
    pooled_anchor = anchor,
    model_type = model,
    randomize_order = FALSE,
    vce_method = "linearized",
    reps = 1,
    seed = 20260726,
    quiet = TRUE
  )
}

test_that("group_levels selects, orders, and audits a multilevel comparison", {
  data <- make_fairlie_phase2_data()
  fit <- fit_fairlie_phase2(data)

  expect_equal(fit$summary_stats$N, 128)
  expect_equal(fit$summary_stats$N0, 64)
  expect_equal(fit$summary_stats$N1, 64)
  expect_equal(fit$summary_stats$n_group_filtered, 16)
  expect_identical(
    fit$summary_stats$selected_group_levels,
    c("advantaged", "disadvantaged")
  )
  expect_identical(fit$raw$group_mapping$selection_mode, "group_levels")
  expect_identical(fit$raw$group_mapping$excluded, "middle")
  expect_identical(
    fit$raw$group_mapping$selected,
    c("advantaged", "disadvantaged")
  )
  expect_equal(length(fit$raw$analytic_sample), 128)
  expect_false(any(data$id[data$group_label == "middle"] %in%
    data$id[fit$raw$analytic_sample]))

  reversed <- fit_fairlie_phase2(
    data,
    levels = c("disadvantaged", "advantaged"),
    anchor = "disadvantaged"
  )
  expect_equal(
    fit$raw$overall_estimates[["Difference"]],
    -reversed$raw$overall_estimates[["Difference"]],
    tolerance = 0
  )
  expect_equal(
    fit$raw$overall_estimates[["Explained"]],
    -reversed$raw$overall_estimates[["Explained"]],
    tolerance = 1e-14
  )
})

test_that("pooled anchors expose the intended group-indicator reference", {
  data <- make_fairlie_phase2_data()
  fits <- lapply(
    c("favored", "disadvantaged", "centered"),
    function(anchor) {
      fit_fairlie_phase2(data, anchor = anchor)
    }
  )
  names(fits) <- c("favored", "disadvantaged", "centered")

  expect_equal(fits$favored$raw$pooled_anchor$reference_value, 0)
  expect_equal(fits$disadvantaged$raw$pooled_anchor$reference_value, 1)
  expect_equal(fits$centered$raw$pooled_anchor$reference_value, 0.5)
  expect_equal(fits$centered$raw$pooled_anchor$disadvantaged_share, 0.5)
  expect_true(all(vapply(
    fits,
    function(fit) isTRUE(fit$raw$pooled_anchor$applicable),
    logical(1)
  )))

  explained <- vapply(
    fits,
    function(fit) fit$raw$overall_estimates[["Explained"]],
    numeric(1)
  )
  expect_gt(diff(range(explained)), 1e-3)
})

test_that("Phase 2 pooled and Neumark points match frozen Stata values", {
  data <- make_fairlie_phase2_data()
  expected <- list(
    logit = c(
      neumark = -0.26467347836553029,
      favored = -0.21086965091970869,
      disadvantaged = -0.21539019246433369,
      centered = -0.21456764532214581
    ),
    probit = c(
      neumark = -0.26747405401914781,
      favored = -0.21183683811720599,
      disadvantaged = -0.21593995013521980,
      centered = -0.21523559327914971
    )
  )

  for (model in names(expected)) {
    for (reference_name in names(expected[[model]])) {
      fit <- if (reference_name == "neumark") {
        fit_fairlie_phase2(
          data,
          model = model,
          reference = "neumark"
        )
      } else {
        fit_fairlie_phase2(
          data,
          model = model,
          anchor = reference_name
        )
      }
      expect_equal(
        fit$raw$overall_estimates[["Explained"]],
        unname(expected[[model]][[reference_name]]),
        tolerance = if (model == "logit") 5e-9 else 5e-8,
        info = paste(model, reference_name)
      )
    }
  }
})

test_that("Reimers and Cotton remain available outside Stata's native scope", {
  data <- make_fairlie_phase2_data()

  for (reference in c("reimers", "cotton")) {
    fit_favored <- fit_fairlie_phase2(
      data,
      reference = reference,
      anchor = "favored"
    )
    fit_centered <- fit_fairlie_phase2(
      data,
      reference = reference,
      anchor = "centered"
    )
    expect_false(fit_favored$raw$pooled_anchor$applicable)
    expect_true(is.na(fit_favored$raw$pooled_anchor$effective))
    expect_equal(
      fit_favored$raw$overall_estimates,
      fit_centered$raw$overall_estimates,
      tolerance = 0
    )
    expect_setequal(
      fit_favored$raw$fitted_models,
      c("group0", "group1")
    )
  }
})

test_that("group_levels validates ambiguous or invalid declarations", {
  data <- make_fairlie_phase2_data()

  expect_error(
    f_decomp(
      data = data,
      dep_var = "y_binary",
      group_var = "group_label",
      group_levels = c("advantaged", "disadvantaged"),
      favored_group = "advantaged",
      indep_vars = c("x1", "x2", "x3"),
      reps = 1,
      quiet = TRUE
    ),
    "either group_levels or favored_group"
  )
  expect_error(
    fit_fairlie_phase2(
      data,
      levels = c("advantaged", "not_observed")
    ),
    "not observed"
  )
  expect_error(
    fit_fairlie_phase2(
      data,
      levels = c("advantaged", "advantaged")
    ),
    "two distinct"
  )
  expect_error(
    fit_fairlie_phase2(data, anchor = "unsupported"),
    "pooled_anchor must be one of"
  )
})
