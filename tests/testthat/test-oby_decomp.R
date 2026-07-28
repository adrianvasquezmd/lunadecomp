make_oby_normalization_data <- function() {
  n <- 600L
  id <- seq_len(n)
  group <- rep(c(0, 1), each = n / 2)
  x1 <- (((id * 17L) %% 101L) - 50) / 10 + 0.6 * group
  x2 <- (((id * 29L) %% 97L) - 48) / 8 - 0.4 * group
  error <- (((id * 61L) %% 83L) - 41) / 13

  insurance_score <- (id * 37L) %% 100L
  insurance_yes <- as.integer(
    insurance_score < ifelse(group == 0, 68L, 42L)
  )
  insurance <- factor(
    ifelse(insurance_yes == 1L, "yes", "no"),
    levels = c("no", "yes")
  )

  education_score <- (id * 53L) %% 100L
  education <- ifelse(
    group == 0,
    ifelse(
      education_score < 25L,
      "low",
      ifelse(education_score < 67L, "middle", "high")
    ),
    ifelse(
      education_score < 48L,
      "low",
      ifelse(education_score < 83L, "middle", "high")
    )
  )
  education <- factor(
    education,
    levels = c("low", "middle", "high")
  )
  education_middle <- as.integer(education == "middle")
  education_high <- as.integer(education == "high")

  y <- 8 +
    1.50 * x1 -
    0.80 * x2 +
    1.20 * insurance_yes +
    0.70 * education_middle +
    1.80 * education_high -
    1.50 * group +
    group * (
      0.25 * x1 -
        0.50 * insurance_yes +
        0.30 * education_middle -
        0.20 * education_high
    ) +
    error

  data.frame(
    group = group,
    y = y,
    x1 = x1,
    x2 = x2,
    insurance = insurance,
    education = education
  )
}

make_oby_survey_design_data <- function() {
  n_strata <- 4L
  n_psu_per_stratum <- 4L
  n_per_psu <- 8L
  n <- n_strata * n_psu_per_stratum * n_per_psu
  id <- seq_len(n)
  strata <- ((id - 1L) %/% (n_psu_per_stratum * n_per_psu)) + 1L
  psu_reused <- (((id - 1L) %/% n_per_psu) %% n_psu_per_stratum) + 1L
  psu_unique <- (strata - 1L) * n_psu_per_stratum + psu_reused
  within_psu <- ((id - 1L) %% n_per_psu) + 1L
  group <- as.integer(within_psu %% 2L == 0L)
  x1 <- (((id * 17L) %% 101L) - 50) / 13 + 0.4 * group
  x2 <- (((id * 29L) %% 97L) - 48) / 11 - 0.3 * group
  error <- (((id * 61L) %% 83L) - 41) / 17
  y <- 6 + 1.2 * x1 - 0.7 * x2 - group +
    0.15 * psu_reused + error
  weight <- 5 + ((id * 13L + strata) %% 17L)
  psu_lonely <- psu_reused
  psu_lonely[strata == 1L] <- 1L

  data.frame(
    group = group,
    y = y,
    x1 = x1,
    x2 = x2,
    weight = weight,
    strata = strata,
    psu_reused = psu_reused,
    psu_unique = psu_unique,
    psu_lonely = psu_lonely
  )
}

test_that("oby_decomp works for OLS decomposition", {
  fit <- oby_decomp(
    data = lunadecomp_example,
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    model_type = "ols",
    ref_method = "pooled",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_s3_class(fit$results_overall, "tbl_df")
  expect_true(all(c("Group_0", "Group_1", "Difference", "Explained", "Unexplained") %in% fit$results_overall$Term))
  expect_true(is.list(fit))
  expect_equal(
    fit$results_overall$Estimate,
    unname(fit$raw$estimates[fit$results_overall$Term]),
    tolerance = 0
  )
  expect_equal(
    fit$results_overall$Std_Error,
    unname(fit$raw$standard_errors[fit$results_overall$Term]),
    tolerance = 0
  )
  expect_equal(
    colnames(fit$raw$vcov),
    names(fit$raw$estimates)
  )
  expect_false("Exp_(Intercept)" %in% names(fit$raw$estimates))
})

test_that("oby_decomp works for linear probability model", {
  fit <- oby_decomp(
    data = lunadecomp_example,
    dep_var = "y_binary",
    group_var = "group",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    model_type = "lpm",
    ref_method = "pooled",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_s3_class(fit$results_overall, "tbl_df")
  expect_true(all(c("Group_0", "Group_1", "Difference", "Explained", "Unexplained") %in% fit$results_overall$Term))
})

test_that("oby_decomp works for logit decomposition", {
  fit <- oby_decomp(
    data = lunadecomp_example,
    dep_var = "y_binary",
    group_var = "group",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    model_type = "logit",
    ref_method = "pooled",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_s3_class(fit$results_overall, "tbl_df")
  expect_true(all(c("Group_0", "Group_1", "Difference", "Explained", "Unexplained") %in% fit$results_overall$Term))
})

make_oby_pooled_anchor_data <- function() {
  n <- 900L
  id <- seq_len(n)
  group <- rep(c(0, 1), each = n / 2)
  x1 <- (((id * 37L) %% 211L) - 105) / 34 + 0.75 * group
  x2 <- (((id * 71L + 13L) %% 307L) - 153) / 53 - 0.35 * group
  eta <- -1.55 + 0.55 * x1 - 0.40 * x2 + 1.10 * group
  uniform_draw <- (((id * 193L + 41L) %% 1009L) + 0.5) / 1009

  data.frame(
    group = group,
    y = as.integer(uniform_draw < stats::plogis(eta)),
    x1 = x1,
    x2 = x2
  )
}

test_that("pooled_anchor selects the nonlinear pooled reference point", {
  validation_data <- make_oby_pooled_anchor_data()
  common_args <- list(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    ref_method = "pooled",
    vce_method = "linearized",
    quiet = TRUE
  )

  for (model_type in c("logit", "probit")) {
    fit_favored <- do.call(
      oby_decomp,
      c(common_args, list(
        model_type = model_type,
        pooled_anchor = "favored"
      ))
    )
    fit_disadvantaged <- do.call(
      oby_decomp,
      c(common_args, list(
        model_type = model_type,
        pooled_anchor = "disadvantaged"
      ))
    )
    fit_centered <- do.call(
      oby_decomp,
      c(common_args, list(
        model_type = model_type,
        pooled_anchor = "centered"
      ))
    )

    overall_names <- c("Group_0", "Group_1", "Difference")
    expect_equal(
      fit_favored$raw$estimates[overall_names],
      fit_disadvantaged$raw$estimates[overall_names],
      tolerance = 1e-12
    )
    expect_equal(
      fit_favored$raw$estimates[overall_names],
      fit_centered$raw$estimates[overall_names],
      tolerance = 1e-12
    )

    delta <- unname(
      fit_favored$models$reference$coefficients[
        ".oby_pooled_group"
      ]
    )
    group_mean <- mean(validation_data$group)
    favored_intercept <- unname(
      fit_favored$raw$reference_coefficients["(Intercept)"]
    )
    expect_equal(
      unname(
        fit_disadvantaged$raw$reference_coefficients["(Intercept)"]
      ),
      favored_intercept + delta,
      tolerance = 1e-9
    )
    expect_equal(
      unname(fit_centered$raw$reference_coefficients["(Intercept)"]),
      favored_intercept + group_mean * delta,
      tolerance = 1e-9
    )
    expect_equal(
      fit_favored$raw$reference_coefficients[c("x1", "x2")],
      fit_disadvantaged$raw$reference_coefficients[c("x1", "x2")],
      tolerance = 1e-9
    )
    expect_equal(
      fit_favored$raw$reference_coefficients[c("x1", "x2")],
      fit_centered$raw$reference_coefficients[c("x1", "x2")],
      tolerance = 1e-9
    )

    explained <- c(
      favored = fit_favored$raw$estimates[["Explained"]],
      disadvantaged =
        fit_disadvantaged$raw$estimates[["Explained"]],
      centered = fit_centered$raw$estimates[["Explained"]]
    )
    expect_gt(max(explained) - min(explained), 1e-5)

    for (fit in list(
      fit_favored,
      fit_disadvantaged,
      fit_centered
    )) {
      expect_equal(
        fit$raw$estimates[["Explained"]] +
          fit$raw$estimates[["Unexplained"]],
        fit$raw$estimates[["Difference"]],
        tolerance = 1e-12
      )
      expect_true(fit$raw$pooled_anchor$applicable)
      expect_true(fit$raw$pooled_anchor$material_nonlinear)
    }

    expect_identical(
      fit_favored$raw$pooled_anchor$effective,
      "favored"
    )
    expect_identical(
      fit_disadvantaged$raw$pooled_anchor$effective,
      "disadvantaged"
    )
    expect_identical(
      fit_centered$raw$pooled_anchor$effective,
      "centered"
    )
    expect_equal(
      fit_centered$raw$pooled_anchor$reference_value,
      group_mean,
      tolerance = 1e-15
    )
  }
})

test_that("pooled_anchor is decomposition-invariant for linear models", {
  validation_data <- make_oby_pooled_anchor_data()
  fit_anchor <- function(anchor) {
    oby_decomp(
      data = validation_data,
      dep_var = "y",
      group_var = "group",
      indep_vars = c("x1", "x2"),
      ref_method = "pooled",
      pooled_anchor = anchor,
      model_type = "lpm",
      vce_method = "linearized",
      quiet = TRUE
    )
  }
  fits <- lapply(
    c("favored", "disadvantaged", "centered"),
    fit_anchor
  )

  expect_equal(
    fits[[1]]$raw$estimates,
    fits[[2]]$raw$estimates,
    tolerance = 1e-12
  )
  expect_equal(
    fits[[1]]$raw$estimates,
    fits[[3]]$raw$estimates,
    tolerance = 1e-12
  )
  expect_equal(
    fits[[1]]$raw$vcov,
    fits[[2]]$raw$vcov,
    tolerance = 1e-11
  )
  expect_equal(
    fits[[1]]$raw$vcov,
    fits[[3]]$raw$vcov,
    tolerance = 1e-11
  )
  expect_false(fits[[1]]$raw$pooled_anchor$material_nonlinear)
})

test_that("pooled_anchor is ignored outside the pooled reference", {
  validation_data <- make_oby_pooled_anchor_data()
  fit_favored <- oby_decomp(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    ref_method = "neumark",
    pooled_anchor = "favored",
    model_type = "logit",
    vce_method = "linearized",
    quiet = TRUE
  )
  fit_disadvantaged <- oby_decomp(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    ref_method = "neumark",
    pooled_anchor = "disadvantaged",
    model_type = "logit",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_equal(
    fit_favored$raw$estimates,
    fit_disadvantaged$raw$estimates,
    tolerance = 0
  )
  expect_equal(
    fit_favored$raw$vcov,
    fit_disadvantaged$raw$vcov,
    tolerance = 0
  )
  expect_false(fit_disadvantaged$raw$pooled_anchor$applicable)
  expect_true(is.na(fit_disadvantaged$raw$pooled_anchor$effective))
})

test_that("pooled_anchor validates its public interface", {
  expect_error(
    oby_decomp(
      data = lunadecomp_example,
      dep_var = "y_binary",
      group_var = "group",
      indep_vars = c("age", "rural"),
      pooled_anchor = "stata_fixed",
      model_type = "logit",
      quiet = TRUE
    ),
    "pooled_anchor must be one of"
  )
})

test_that("centered pooled anchor supports replicate VCE", {
  common_args <- list(
    data = lunadecomp_example,
    dep_var = "y_binary",
    group_var = "group",
    indep_vars = c("age", "rural"),
    ref_method = "pooled",
    pooled_anchor = "centered",
    model_type = "logit",
    seed = 42,
    relax = TRUE,
    quiet = TRUE
  )
  fit_jackknife <- do.call(
    oby_decomp,
    c(common_args, list(vce_method = "jackknife"))
  )
  fit_bootstrap <- do.call(
    oby_decomp,
    c(
      common_args,
      list(vce_method = "bootstrap", boot_reps = 8)
    )
  )

  expect_true(all(is.finite(fit_jackknife$results_overall$Std_Error)))
  expect_true(all(is.finite(fit_bootstrap$results_overall$Std_Error)))
  expect_identical(
    fit_jackknife$raw$pooled_anchor$effective,
    "centered"
  )
  expect_identical(
    fit_bootstrap$raw$pooled_anchor$effective,
    "centered"
  )
})

test_that("nonlinear linearized VCE uses Stata-compatible ML sandwiches", {
  n <- 600L
  id <- seq_len(n)
  group <- rep(c(0, 1), each = n / 2)
  x1 <- (((id * 37L) %% 211L) - 105) / 31 + 0.25 * group
  x2 <- (((id * 71L + 13L) %% 307L) - 153) / 46 - 0.15 * group
  eta <- -0.20 + 0.40 * x1 - 0.30 * x2 - 0.35 * group
  uniform_draw <- (((id * 193L + 41L) %% 1009L) + 0.5) / 1009
  validation_data <- data.frame(
    group = group,
    y = as.integer(uniform_draw < stats::plogis(eta)),
    x1 = x1,
    x2 = x2
  )

  manual_ml_vcov <- function(model_type) {
    group_data <- validation_data[validation_data$group == 0, ]
    X <- stats::model.matrix(~ x1 + x2, group_data)
    family <- stats::binomial(link = model_type)
    model <- stats::glm.fit(
      X,
      group_data$y,
      family = family,
      control = stats::glm.control(epsilon = 1e-16, maxit = 100)
    )
    mu <- pmin(pmax(model$fitted.values, 1e-15), 1 - 1e-15)
    eta <- model$linear.predictors
    variance <- mu * (1 - mu)
    q <- if (model_type == "logit") {
      variance
    } else {
      stats::dnorm(eta)
    }
    residual <- group_data$y - mu
    score <- X * as.numeric(residual * q / variance)
    hessian_weights <- q^2 / variance
    if (model_type == "probit") {
      q_prime <- -eta * q
      hessian_weights <- hessian_weights -
        residual * q_prime / variance +
        residual * q^2 * (1 - 2 * mu) / variance^2
    }
    bread <- MASS::ginv(crossprod(X, X * hessian_weights))
    n_group <- nrow(group_data)
    (n_group / (n_group - 1)) *
      bread %*% crossprod(score) %*% bread
  }

  for (model_type in c("logit", "probit")) {
    fit <- oby_decomp(
      data = validation_data,
      dep_var = "y",
      group_var = "group",
      indep_vars = c("x1", "x2"),
      ref_method = "group1",
      model_type = model_type,
      vce_method = "linearized",
      quiet = TRUE
    )
    terms <- paste0("b0:", c("(Intercept)", "x1", "x2"))
    expect_equal(
      unname(fit$raw$base_vcov[terms, terms]),
      unname(manual_ml_vcov(model_type)),
      tolerance = 1e-11
    )
  }

  probit_data <- validation_data[validation_data$group == 0, ]
  X <- stats::model.matrix(~ x1 + x2, probit_data)
  model <- stats::glm.fit(
    X,
    probit_data$y,
    family = stats::binomial(link = "probit"),
    control = stats::glm.control(epsilon = 1e-16, maxit = 100)
  )
  mu <- model$fitted.values
  fisher_bread <- MASS::ginv(crossprod(
    X,
    X * as.numeric(stats::dnorm(model$linear.predictors)^2 / (mu * (1 - mu)))
  ))
  fisher_vcov <- (nrow(X) / (nrow(X) - 1)) *
    fisher_bread %*%
    crossprod(X * as.numeric(
      (probit_data$y - mu) * stats::dnorm(model$linear.predictors) /
        (mu * (1 - mu))
    )) %*%
    fisher_bread
  expect_gt(
    max(abs(manual_ml_vcov("probit") - fisher_vcov)),
    1e-8
  )
})

test_that("group_levels selects and orders two levels from a multilevel group", {
  group_sizes <- c(
    lowest = 90L,
    lower_middle = 70L,
    middle = 80L,
    upper_middle = 75L,
    highest = 110L
  )
  group <- rep(names(group_sizes), times = group_sizes)
  id <- seq_along(group)
  x1 <- (((id * 31L) %% 127L) - 63) / 19
  x2 <- sin(id / 13)
  group_effect <- unname(c(
    lowest = -1.4,
    lower_middle = -0.7,
    middle = 0,
    upper_middle = 0.8,
    highest = 1.6
  )[group])
  validation_data <- data.frame(
    group = group,
    y = 4 + 0.7 * x1 - 0.4 * x2 + group_effect,
    x1 = x1,
    x2 = x2,
    stringsAsFactors = FALSE
  )

  fit <- oby_decomp(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    group_levels = c("highest", "lowest"),
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )

  manual_data <- validation_data[
    validation_data$group %in% c("highest", "lowest"),
  ]
  manual_data$group <- ifelse(manual_data$group == "highest", 0, 1)
  manual_fit <- oby_decomp(
    data = manual_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_equal(fit$raw$estimates, manual_fit$raw$estimates, tolerance = 0)
  expect_equal(fit$raw$vcov, manual_fit$raw$vcov, tolerance = 0)
  expect_identical(fit$summary_stats$favored, "highest")
  expect_identical(fit$summary_stats$disadvantaged, "lowest")
  expect_identical(fit$summary_stats$group_selection, "group_levels")
  expect_equal(fit$summary_stats$N0, unname(group_sizes[["highest"]]))
  expect_equal(fit$summary_stats$N1, unname(group_sizes[["lowest"]]))
  expect_equal(
    fit$summary_stats$filtered_observations,
    sum(group_sizes) - group_sizes[["highest"]] - group_sizes[["lowest"]]
  )
  expect_setequal(
    fit$summary_stats$filtered_groups,
    c("lower_middle", "middle", "upper_middle")
  )
  expect_identical(
    fit$raw$group_mapping$selected,
    c("highest", "lowest")
  )
  expect_identical(
    fit$raw$analytic_sample,
    which(validation_data$group %in% c("highest", "lowest"))
  )

  reverse_fit <- oby_decomp(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    group_levels = c("lowest", "highest"),
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )
  expect_equal(
    fit$raw$estimates[["Group_0"]],
    reverse_fit$raw$estimates[["Group_1"]],
    tolerance = 1e-12
  )
  expect_equal(
    fit$raw$estimates[["Group_1"]],
    reverse_fit$raw$estimates[["Group_0"]],
    tolerance = 1e-12
  )
  expect_equal(
    fit$raw$estimates[["Difference"]],
    -reverse_fit$raw$estimates[["Difference"]],
    tolerance = 1e-12
  )

  two_group_data <- validation_data[
    validation_data$group %in% c("highest", "lowest"),
  ]
  legacy_fit <- oby_decomp(
    data = two_group_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    favored_group = "highest",
    ref_method = "pooled",
    quiet = TRUE
  )
  expect_equal(
    fit$raw$estimates,
    legacy_fit$raw$estimates,
    tolerance = 0
  )
})

test_that("group_levels validates explicit group comparisons", {
  validation_data <- data.frame(
    group = rep(c("low", "middle", "high"), each = 20),
    y = seq_len(60),
    x = rep(seq_len(20), 3)
  )
  common_args <- list(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = "x",
    quiet = TRUE
  )

  expect_error(
    do.call(
      oby_decomp,
      c(common_args, list(group_levels = "high"))
    ),
    "ordered vector of length two"
  )
  expect_error(
    do.call(
      oby_decomp,
      c(common_args, list(group_levels = c("high", "high")))
    ),
    "two distinct values"
  )
  expect_error(
    do.call(
      oby_decomp,
      c(common_args, list(group_levels = c("high", "missing")))
    ),
    "not observed"
  )
  expect_error(
    do.call(
      oby_decomp,
      c(
        common_args,
        list(
          group_levels = c("high", "low"),
          favored_group = "high"
        )
      )
    ),
    "not both"
  )
  expect_error(
    do.call(
      oby_decomp,
      c(common_args, list(favored_group = "high"))
    ),
    "Use group_levels"
  )
})

test_that("oby_decomp exposes the approved VCE interface", {
  expect_false("robust" %in% names(formals(oby_decomp)))

  expect_error(
    oby_decomp(
      data = lunadecomp_example,
      dep_var = "y_continuous",
      group_var = "group",
      indep_vars = c("age", "rural"),
      vce_method = "analytic",
      quiet = TRUE
    ),
    "vce_method must be one of"
  )

  expect_error(
    oby_decomp(
      data = lunadecomp_example,
      dep_var = "y_continuous",
      group_var = "group",
      indep_vars = c("age", "rural"),
      vce_method = "svy_linearized",
      quiet = TRUE
    ),
    "vce_method must be one of"
  )

  fit <- oby_decomp(
    data = lunadecomp_example,
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("age", "rural"),
    use_svy = TRUE,
    weight_var = "weight",
    strata_var = "strata",
    psu_var = "psu",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_identical(fit$raw$vce_method, "linearized")
  expect_true(fit$raw$survey_mode)
})

test_that("oby_decomp validates sampling weights and survey degrees of freedom", {
  validation_data <- lunadecomp_example

  fit_weighted <- oby_decomp(
    data = validation_data,
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("age", "rural"),
    weight_var = "weight",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_equal(
    fit_weighted$summary_stats$DF,
    fit_weighted$summary_stats$PSUs -
      fit_weighted$summary_stats$Strata
  )

  zero_weight_data <- validation_data
  zero_weight_data$weight[1] <- 0
  fit_zero <- oby_decomp(
    data = zero_weight_data,
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("age", "rural"),
    weight_var = "weight",
    quiet = TRUE
  )
  expect_identical(
    fit_zero$summary_stats$N,
    nrow(zero_weight_data)
  )

  negative_weight_data <- validation_data
  negative_weight_data$weight[1] <- -1
  expect_error(
    oby_decomp(
      data = negative_weight_data,
      dep_var = "y_continuous",
      group_var = "group",
      indep_vars = c("age", "rural"),
      weight_var = "weight",
      quiet = TRUE
    ),
    "non-negative"
  )

  infinite_weight_data <- validation_data
  infinite_weight_data$weight[1] <- Inf
  expect_error(
    oby_decomp(
      data = infinite_weight_data,
      dep_var = "y_continuous",
      group_var = "group",
      indep_vars = c("age", "rural"),
      weight_var = "weight",
      quiet = TRUE
    ),
    "finite"
  )

  character_weight_data <- validation_data
  character_weight_data$weight_text <- as.character(
    character_weight_data$weight
  )
  expect_error(
    oby_decomp(
      data = character_weight_data,
      dep_var = "y_continuous",
      group_var = "group",
      indep_vars = c("age", "rural"),
      weight_var = "weight_text",
      quiet = TRUE
    ),
    "numeric sampling-weight"
  )
})

test_that("sampling-weight scale and Cotton reference are invariant", {
  validation_data <- lunadecomp_example
  scaled_data <- validation_data
  scaled_data$weight <- 100 * scaled_data$weight

  common_args <- list(
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("age", "rural"),
    weight_var = "weight",
    ref_method = "cotton",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )
  fit <- do.call(
    oby_decomp,
    c(list(data = validation_data), common_args)
  )
  fit_scaled <- do.call(
    oby_decomp,
    c(list(data = scaled_data), common_args)
  )

  expect_equal(
    fit$raw$estimates,
    fit_scaled$raw$estimates,
    tolerance = 1e-12
  )
  expect_equal(
    fit$raw$vcov,
    fit_scaled$raw$vcov,
    tolerance = 1e-12
  )

  favored_share <- with(
    validation_data,
    sum(weight[group == 0]) / sum(weight)
  )
  expected_reference <-
    favored_share * fit$raw$group_coefficients$group_0 +
    (1 - favored_share) * fit$raw$group_coefficients$group_1
  expect_equal(
    fit$raw$reference_coefficients,
    expected_reference,
    tolerance = 1e-12
  )
})

test_that("survey PSU nesting and lonely-stratum rules are explicit", {
  validation_data <- make_oby_survey_design_data()
  common_args <- list(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    weight_var = "weight",
    strata_var = "strata",
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )

  fit_unique <- do.call(
    oby_decomp,
    c(common_args, list(psu_var = "psu_unique"))
  )
  fit_reused <- do.call(
    oby_decomp,
    c(common_args, list(psu_var = "psu_reused"))
  )
  expect_equal(
    fit_unique$raw$estimates,
    fit_reused$raw$estimates,
    tolerance = 0
  )
  expect_equal(
    fit_unique$raw$vcov,
    fit_reused$raw$vcov,
    tolerance = 1e-12
  )
  expect_equal(fit_unique$summary_stats$PSUs, 16)
  expect_equal(fit_reused$summary_stats$PSUs, 16)
  expect_equal(fit_reused$summary_stats$DF, 12)

  lonely_args <- c(common_args, list(psu_var = "psu_lonely"))
  fit_adjust <- do.call(
    oby_decomp,
    c(lonely_args, list(lonely_psu = "adjust"))
  )
  fit_certainty <- do.call(
    oby_decomp,
    c(lonely_args, list(lonely_psu = "certainty"))
  )
  fit_remove <- do.call(
    oby_decomp,
    c(lonely_args, list(lonely_psu = "remove"))
  )
  fit_average <- do.call(
    oby_decomp,
    c(lonely_args, list(lonely_psu = "average"))
  )

  expect_equal(
    fit_certainty$raw$vcov,
    fit_remove$raw$vcov,
    tolerance = 0
  )
  expect_true(all(is.finite(fit_adjust$raw$vcov)))
  expect_true(all(is.finite(fit_average$raw$vcov)))
  expect_gt(
    max(abs(fit_adjust$raw$vcov - fit_certainty$raw$vcov)),
    1e-12
  )
  expect_gt(
    max(abs(fit_average$raw$vcov - fit_certainty$raw$vcov)),
    1e-12
  )
  expect_error(
    do.call(
      oby_decomp,
      c(lonely_args, list(lonely_psu = "fail"))
    ),
    "single PSU"
  )
})

test_that("reported inference uses stable normal and survey-t tails", {
  validation_data <- make_oby_survey_design_data()

  ordinary_fit <- oby_decomp(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    level = 0.90,
    quiet = TRUE
  )
  ordinary <- ordinary_fit$results_overall
  ordinary_stat <- ordinary$Estimate / ordinary$Std_Error
  ordinary_critical <- stats::qnorm(0.95)
  expect_true(is.infinite(ordinary_fit$summary_stats$DF))
  expect_equal(ordinary$Statistic, ordinary_stat, tolerance = 0)
  expect_equal(
    ordinary$P_Value,
    2 * stats::pnorm(abs(ordinary_stat), lower.tail = FALSE),
    tolerance = 1e-15
  )
  expect_equal(
    ordinary$Conf_Low,
    ordinary$Estimate - ordinary_critical * ordinary$Std_Error,
    tolerance = 1e-15
  )
  expect_equal(
    ordinary$Conf_High,
    ordinary$Estimate + ordinary_critical * ordinary$Std_Error,
    tolerance = 1e-15
  )

  survey_fit <- oby_decomp(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    weight_var = "weight",
    strata_var = "strata",
    psu_var = "psu_unique",
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    level = 0.90,
    quiet = TRUE
  )
  survey <- survey_fit$results_overall
  survey_df <- survey_fit$summary_stats$DF
  survey_stat <- survey$Estimate / survey$Std_Error
  survey_critical <- stats::qt(0.95, df = survey_df)
  expect_equal(survey_df, 12)
  expect_equal(
    survey$P_Value,
    2 * stats::pt(
      abs(survey_stat),
      df = survey_df,
      lower.tail = FALSE
    ),
    tolerance = 1e-15
  )
  expect_equal(
    survey$Conf_Low,
    survey$Estimate - survey_critical * survey$Std_Error,
    tolerance = 1e-15
  )
  expect_equal(
    survey$Conf_High,
    survey$Estimate + survey_critical * survey$Std_Error,
    tolerance = 1e-15
  )
  expect_true(all(survey$P_Value > 0))
})

test_that("oby_decomp supports jackknife and bootstrap VCE", {
  common_args <- list(
    data = lunadecomp_example,
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("age", "rural"),
    model_type = "ols",
    ref_method = "pooled",
    seed = 123,
    quiet = TRUE
  )

  fit_jk <- do.call(
    oby_decomp,
    c(common_args, list(vce_method = "jackknife"))
  )
  fit_boot <- do.call(
    oby_decomp,
    c(common_args, list(vce_method = "bootstrap", boot_reps = 10))
  )

  expect_true(all(is.finite(fit_jk$results_overall$Std_Error)))
  expect_true(all(is.finite(fit_boot$results_overall$Std_Error)))
  expect_identical(fit_jk$raw$vce_method, "jackknife")
  expect_identical(fit_boot$raw$vce_method, "bootstrap")
})

test_that("replicate VCE uses Stata-compatible scales and centering", {
  validation_data <- make_oby_survey_design_data()
  common_args <- list(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2"),
    model_type = "ols",
    ref_method = "pooled",
    quiet = TRUE
  )

  set.seed(271800)
  rng_before_jackknife <- .Random.seed
  ordinary_jk <- do.call(
    oby_decomp,
    c(common_args, list(vce_method = "jackknife"))
  )
  expect_identical(.Random.seed, rng_before_jackknife)
  ordinary_boot <- do.call(
    oby_decomp,
    c(common_args, list(
      vce_method = "bootstrap",
      boot_reps = 12,
      seed = 271801
    ))
  )
  survey_args <- c(
    common_args,
    list(
      weight_var = "weight",
      strata_var = "strata",
      psu_var = "psu_unique"
    )
  )
  survey_jk <- do.call(
    oby_decomp,
    c(survey_args, list(vce_method = "jackknife"))
  )
  survey_boot <- do.call(
    oby_decomp,
    c(survey_args, list(
      vce_method = "bootstrap",
      boot_reps = 12,
      seed = 271803
    ))
  )

  expect_equal(ordinary_jk$summary_stats$DF, nrow(validation_data) - 1)
  expect_true(is.infinite(ordinary_boot$summary_stats$DF))
  expect_equal(survey_jk$summary_stats$DF, 12)
  expect_equal(survey_boot$summary_stats$DF, 12)

  expect_identical(
    ordinary_jk$raw$replication$center_method,
    "replicate mean"
  )
  expect_identical(
    ordinary_boot$raw$replication$center_method,
    "replicate mean"
  )
  expect_identical(
    ordinary_boot$raw$replication$engine,
    "boot::boot"
  )
  expect_identical(
    ordinary_boot$raw$replication$simulation,
    "ordinary"
  )
  expect_identical(
    ordinary_boot$raw$replication$statistic_type,
    "indices"
  )
  expect_identical(
    ordinary_boot$raw$replication$fit_weight_type,
    "frequency"
  )
  expect_identical(
    survey_jk$raw$replication$center_method,
    "stratum-specific replicate mean"
  )
  expect_identical(
    survey_boot$raw$replication$center_method,
    "replicate mean"
  )
  expect_identical(
    ordinary_jk$raw$replication$engine,
    "resample::jackknife"
  )
  expect_true(
    ordinary_jk$raw$replication$ordinary_jackknife
  )
  expect_identical(
    ordinary_jk$raw$replication$simulation,
    "delete-one"
  )
  expect_identical(
    ordinary_jk$raw$replication$statistic_type,
    "retained row indices"
  )
  expect_identical(
    ordinary_jk$raw$replication$fit_weight_type,
    "JK1 replicate factor"
  )
  expect_equal(
    nrow(ordinary_jk$raw$replication$replicates),
    nrow(validation_data)
  )
  expect_equal(
    ordinary_jk$raw$replication$requested_replicates,
    nrow(validation_data)
  )
  expect_identical(
    survey_jk$raw$replication$engine,
    "survey::withReplicates"
  )
  expect_identical(
    survey_boot$raw$replication$engine,
    "survey::withReplicates"
  )
  expect_identical(
    survey_boot$raw$replication$survey_replicate_type,
    "Rao-Wu rescaled bootstrap"
  )

  expect_equal(
    ordinary_jk$raw$replication$scale,
    (nrow(validation_data) - 1) / nrow(validation_data),
    tolerance = 0
  )
  expect_equal(
    ordinary_boot$raw$replication$scale,
    1 / 11,
    tolerance = 0
  )
  expect_equal(
    survey_jk$raw$replication$scale,
    1,
    tolerance = 0
  )
  expect_equal(
    survey_boot$raw$replication$scale,
    1 / 12,
    tolerance = 0
  )
  expect_false(
    isTRUE(all.equal(
      survey_boot$raw$replication$scale,
      survey_boot$raw$replication$generator_scale
    ))
  )
  expect_equal(
    survey_boot$raw$replication$generator_scale,
    1 / 11,
    tolerance = 0
  )

  reconstruct_vcov <- function(fit) {
    replication <- fit$raw$replication
    differences <- if (is.matrix(replication$center)) {
      replication$replicates - replication$center
    } else {
      sweep(
        replication$replicates,
        2,
        replication$center,
        "-"
      )
    }
    replication$scale * crossprod(
      differences * sqrt(replication$rscales)
    )
  }
  for (fit in list(
    ordinary_jk, ordinary_boot, survey_jk, survey_boot
  )) {
    expect_equal(
      reconstruct_vcov(fit),
      unname(fit$raw$vcov),
      tolerance = 1e-14,
      ignore_attr = TRUE
    )
    expect_equal(fit$raw$replication$failed_replicates, 0)
  }

  ordinary_boot_repeat <- do.call(
    oby_decomp,
    c(common_args, list(
      vce_method = "bootstrap",
      boot_reps = 12,
      seed = 271801
    ))
  )
  expect_equal(
    ordinary_boot_repeat$raw$vcov,
    ordinary_boot$raw$vcov,
    tolerance = 0
  )
})

test_that("Rao-Wu survey bootstrap rejects singleton PSU strata", {
  validation_data <- make_oby_survey_design_data()
  validation_data$singleton_strata <- ifelse(
    validation_data$strata == 1,
    "singleton",
    paste0(
      "other_",
      validation_data$strata,
      "_",
      validation_data$psu_unique
    )
  )
  validation_data$singleton_psu <- ifelse(
    validation_data$strata == 1,
    "one",
    validation_data$psu_unique
  )

  expect_error(
    oby_decomp(
      data = validation_data,
      dep_var = "y",
      group_var = "group",
      indep_vars = c("x1", "x2"),
      weight_var = "weight",
      strata_var = "singleton_strata",
      psu_var = "singleton_psu",
      vce_method = "bootstrap",
      boot_reps = 8,
      seed = 271804,
      quiet = TRUE
    ),
    "Rao-Wu survey bootstrap requires at least two PSUs in every stratum"
  )
})

test_that("oby_decomp handles and groups multi-level categorical predictors", {
  fit <- oby_decomp(
    data = lunadecomp_example,
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("age", "rural", "education"),
    groupings = list(
      Residence = "rural",
      Education = "education"
    ),
    model_type = "ols",
    ref_method = "pooled",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(all(c(
    "ruralrural",
    "educationsecondary",
    "educationhigher"
  ) %in% fit$raw$model_terms))

  education_exp_terms <- c(
    "Exp_educationsecondary",
    "Exp_educationhigher"
  )
  education_unexp_terms <- c(
    "Unexp_educationsecondary",
    "Unexp_educationhigher"
  )

  grouped_exp <- fit$results_grouped_explained[
    fit$results_grouped_explained$Term == "Education",
  ]
  grouped_unexp <- fit$results_grouped_unexplained[
    fit$results_grouped_unexplained$Term == "Education",
  ]

  expect_equal(
    grouped_exp$Estimate,
    sum(fit$raw$estimates[education_exp_terms]),
    tolerance = 1e-12
  )
  expect_equal(
    grouped_exp$Std_Error,
    sqrt(sum(fit$raw$vcov[education_exp_terms, education_exp_terms])),
    tolerance = 1e-12
  )
  expect_equal(
    grouped_unexp$Estimate,
    sum(fit$raw$estimates[education_unexp_terms]),
    tolerance = 1e-12
  )
  expect_equal(
    grouped_unexp$Std_Error,
    sqrt(sum(fit$raw$vcov[education_unexp_terms, education_unexp_terms])),
    tolerance = 1e-12
  )
})

test_that("oby_decomp normalizes selected or all factor predictors", {
  validation_data <- make_oby_normalization_data()
  common_args <- list(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2", "insurance", "education"),
    groupings = list(
      insurance = "insurance",
      education = "education"
    ),
    model_type = "ols",
    ref_method = "pooled",
    vce_method = "linearized",
    quiet = TRUE
  )

  fit_none <- do.call(oby_decomp, common_args)
  fit_all <- do.call(
    oby_decomp,
    c(common_args, list(normalize = "all"))
  )
  fit_named <- do.call(
    oby_decomp,
    c(common_args, list(normalize = c("insurance", "education")))
  )
  fit_education <- do.call(
    oby_decomp,
    c(common_args, list(normalize = "education"))
  )

  expect_identical(
    fit_all$raw$normalization$factors,
    c("insurance", "education")
  )
  expect_true(fit_all$raw$normalization$active)
  expect_equal(
    fit_all$raw$estimates,
    fit_named$raw$estimates,
    tolerance = 0
  )
  expect_equal(
    fit_all$raw$vcov,
    fit_named$raw$vcov,
    tolerance = 0
  )
  expect_equal(
    fit_all$raw$estimates[c(
      "Group_0", "Group_1", "Difference", "Explained", "Unexplained"
    )],
    fit_none$raw$estimates[c(
      "Group_0", "Group_1", "Difference", "Explained", "Unexplained"
    )],
    tolerance = 1e-12
  )

  expect_true(all(c(
    "insuranceno", "insuranceyes",
    "educationlow", "educationmiddle", "educationhigh"
  ) %in% fit_all$raw$model_terms))
  expect_false("insuranceno" %in% fit_education$raw$model_terms)
  expect_true("educationlow" %in% fit_education$raw$model_terms)

  for (coefs in c(
    fit_all$raw$group_coefficients,
    list(reference = fit_all$raw$reference_coefficients)
  )) {
    expect_equal(
      sum(coefs[c("insuranceno", "insuranceyes")]),
      0,
      tolerance = 1e-12
    )
    expect_equal(
      sum(coefs[c(
        "educationlow", "educationmiddle", "educationhigh"
      )]),
      0,
      tolerance = 1e-12
    )
  }

  for (means in fit_all$raw$covariate_means) {
    expect_equal(
      sum(means[c("insuranceno", "insuranceyes")]),
      1,
      tolerance = 1e-12
    )
    expect_equal(
      sum(means[c(
        "educationlow", "educationmiddle", "educationhigh"
      )]),
      1,
      tolerance = 1e-12
    )
  }

  insurance_terms <- paste0(
    "Unexp_",
    c("insuranceno", "insuranceyes")
  )
  grouped_insurance <- fit_all$results_grouped_unexplained[
    fit_all$results_grouped_unexplained$Term == "insurance",
  ]
  expect_equal(
    grouped_insurance$Estimate,
    sum(fit_all$raw$estimates[insurance_terms]),
    tolerance = 1e-12
  )
  expect_equal(
    grouped_insurance$Std_Error,
    sqrt(sum(fit_all$raw$vcov[
      insurance_terms,
      insurance_terms
    ])),
    tolerance = 1e-12
  )
})

test_that("normalized categorical results are invariant to omitted levels", {
  validation_data <- make_oby_normalization_data()
  changed_base <- validation_data
  changed_base$insurance <- stats::relevel(
    changed_base$insurance,
    ref = "yes"
  )
  changed_base$education <- stats::relevel(
    changed_base$education,
    ref = "high"
  )

  fit_original <- oby_decomp(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2", "insurance", "education"),
    normalize = "all",
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )
  fit_changed <- oby_decomp(
    data = changed_base,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2", "insurance", "education"),
    normalize = "all",
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )

  common_terms <- intersect(
    names(fit_original$raw$estimates),
    names(fit_changed$raw$estimates)
  )
  expect_setequal(
    names(fit_original$raw$estimates),
    names(fit_changed$raw$estimates)
  )
  expect_equal(
    fit_original$raw$estimates[common_terms],
    fit_changed$raw$estimates[common_terms],
    tolerance = 1e-12
  )
  expect_equal(
    fit_original$raw$standard_errors[common_terms],
    fit_changed$raw$standard_errors[common_terms],
    tolerance = 1e-12
  )
  expect_equal(
    fit_original$raw$vcov[common_terms, common_terms],
    fit_changed$raw$vcov[common_terms, common_terms],
    tolerance = 1e-12
  )
})

test_that("oby_decomp normalizes factor-by-continuous interactions", {
  validation_data <- make_oby_normalization_data()
  common_args <- list(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c(
      "x2",
      "insurance * x1",
      "education"
    ),
    normalize = "insurance",
    groupings = list(
      insurance = "insurance",
      insurance_x1 = "insurance:x1",
      insurance_block = c(
        "insurance",
        "x1",
        "insurance:x1"
      )
    ),
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )

  fit <- do.call(oby_decomp, common_args)
  fit_unnormalized <- do.call(
    oby_decomp,
    within(common_args, rm(normalize))
  )
  fit_expanded <- oby_decomp(
    data = validation_data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c(
      "x2",
      "insurance",
      "x1",
      "education",
      "insurance:x1"
    ),
    normalize = "insurance",
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )

  interaction_terms <- c(
    "insuranceno:x1",
    "insuranceyes:x1"
  )
  expect_true(all(interaction_terms %in% fit$raw$model_terms))
  expect_identical(
    fit$raw$normalization$blocks$insurance$
      interactions$`insurance:x1`$continuous,
    "x1"
  )
  expect_equal(
    fit$raw$estimates,
    fit_expanded$raw$estimates,
    tolerance = 0
  )
  expect_equal(
    fit$raw$vcov,
    fit_expanded$raw$vcov,
    tolerance = 0
  )
  expect_equal(
    fit$raw$estimates[c(
      "Group_0", "Group_1", "Difference", "Explained", "Unexplained"
    )],
    fit_unnormalized$raw$estimates[c(
      "Group_0", "Group_1", "Difference", "Explained", "Unexplained"
    )],
    tolerance = 1e-12
  )

  for (coefs in c(
    fit$raw$group_coefficients,
    list(reference = fit$raw$reference_coefficients)
  )) {
    expect_equal(
      sum(coefs[c("insuranceno", "insuranceyes")]),
      0,
      tolerance = 1e-12
    )
    expect_equal(
      sum(coefs[interaction_terms]),
      0,
      tolerance = 1e-12
    )
  }
  for (means in fit$raw$covariate_means) {
    expect_equal(
      sum(means[c("insuranceno", "insuranceyes")]),
      1,
      tolerance = 1e-12
    )
    expect_equal(
      sum(means[interaction_terms]),
      means[["x1"]],
      tolerance = 1e-12
    )
  }

  grouped_interaction <- fit$results_grouped_unexplained[
    fit$results_grouped_unexplained$Term == "insurance_x1",
  ]
  raw_interaction_terms <- paste0("Unexp_", interaction_terms)
  expect_equal(
    grouped_interaction$Estimate,
    sum(fit$raw$estimates[raw_interaction_terms]),
    tolerance = 1e-12
  )
  expect_equal(
    grouped_interaction$Std_Error,
    sqrt(sum(fit$raw$vcov[
      raw_interaction_terms,
      raw_interaction_terms
    ])),
    tolerance = 1e-12
  )
})

test_that("normalized interactions are invariant to the omitted factor level", {
  validation_data <- make_oby_normalization_data()
  changed_base <- validation_data
  changed_base$insurance <- stats::relevel(
    changed_base$insurance,
    ref = "yes"
  )

  common_args <- list(
    dep_var = "y",
    group_var = "group",
    indep_vars = c(
      "x2",
      "insurance * x1",
      "education"
    ),
    normalize = "insurance",
    ref_method = "pooled",
    model_type = "ols",
    vce_method = "linearized",
    quiet = TRUE
  )
  fit_original <- do.call(
    oby_decomp,
    c(list(data = validation_data), common_args)
  )
  fit_changed <- do.call(
    oby_decomp,
    c(list(data = changed_base), common_args)
  )

  common_terms <- intersect(
    names(fit_original$raw$estimates),
    names(fit_changed$raw$estimates)
  )
  expect_setequal(
    names(fit_original$raw$estimates),
    names(fit_changed$raw$estimates)
  )
  expect_equal(
    fit_original$raw$estimates[common_terms],
    fit_changed$raw$estimates[common_terms],
    tolerance = 1e-12
  )
  expect_equal(
    fit_original$raw$standard_errors[common_terms],
    fit_changed$raw$standard_errors[common_terms],
    tolerance = 1e-12
  )
  expect_equal(
    fit_original$raw$vcov[common_terms, common_terms],
    fit_changed$raw$vcov[common_terms, common_terms],
    tolerance = 1e-12
  )
})

test_that("oby_decomp validates normalization specifications", {
  common_args <- list(
    data = lunadecomp_example,
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("age", "rural", "education"),
    quiet = TRUE
  )

  expect_error(
    do.call(
      oby_decomp,
      c(common_args, list(normalize = c("all", "education")))
    ),
    "must be used by itself"
  )
  expect_error(
    do.call(
      oby_decomp,
      c(common_args, list(normalize = "missing_factor"))
    ),
    "not included in indep_vars"
  )
  expect_error(
    do.call(
      oby_decomp,
      c(common_args, list(normalize = "age"))
    ),
    "only include factor or character"
  )
  expect_error(
    do.call(
      oby_decomp,
      c(common_args, list(normalize = TRUE))
    ),
    "normalize must be NULL"
  )

  interaction_data <- make_oby_normalization_data()
  expect_error(
    oby_decomp(
      data = interaction_data,
      dep_var = "y",
      group_var = "group",
      indep_vars = "insurance:x1",
      normalize = "insurance",
      quiet = TRUE
    ),
    "must enter the model with its main effect"
  )
  expect_error(
    oby_decomp(
      data = interaction_data,
      dep_var = "y",
      group_var = "group",
      indep_vars = c("insurance", "insurance:x1"),
      normalize = "insurance",
      quiet = TRUE
    ),
    "requires the main effect of continuous variable"
  )
  expect_error(
    oby_decomp(
      data = interaction_data,
      dep_var = "y",
      group_var = "group",
      indep_vars = "insurance * education",
      normalize = "all",
      quiet = TRUE
    ),
    "supports interactions between one normalized factor and one continuous"
  )
})

test_that("replicate VCE normalizes factors within every replicate", {
  common_args <- list(
    data = lunadecomp_example,
    dep_var = "y_continuous",
    group_var = "group",
    indep_vars = c("rural * age", "education"),
    normalize = "all",
    model_type = "ols",
    ref_method = "pooled",
    seed = 321,
    quiet = TRUE
  )

  fit_jk <- do.call(
    oby_decomp,
    c(common_args, list(vce_method = "jackknife"))
  )
  fit_boot <- do.call(
    oby_decomp,
    c(common_args, list(
      vce_method = "bootstrap",
      boot_reps = 10
    ))
  )

  normalized_terms <- c(
    "Unexp_ruralurban", "Unexp_ruralrural",
    "Unexp_ruralurban:age", "Unexp_ruralrural:age",
    "Unexp_educationprimary", "Unexp_educationsecondary",
    "Unexp_educationhigher"
  )
  expect_true(all(normalized_terms %in% names(fit_jk$raw$estimates)))
  expect_true(all(normalized_terms %in% names(fit_boot$raw$estimates)))
  expect_true(all(is.finite(
    fit_jk$raw$standard_errors[normalized_terms]
  )))
  expect_true(all(is.finite(
    fit_boot$raw$standard_errors[normalized_terms]
  )))
})

test_that("oby_decomp rejects invalid model_type", {
  expect_error(
    oby_decomp(
      data = lunadecomp_example,
      dep_var = "y_binary",
      group_var = "group",
      indep_vars = c("age", "rural"),
      model_type = "linear",
      quiet = TRUE
    ),
    "model_type must be one of"
  )
})
