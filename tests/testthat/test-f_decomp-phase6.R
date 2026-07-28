make_fairlie_phase6_data <- function() {
  set.seed(20260730)
  n_per_group <- 60L
  n <- 2L * n_per_group
  id <- seq_len(n)
  group <- rep(0:1, each = n_per_group)
  x_cont <- as.numeric(scale(
    seq(-2.8, 2.8, length.out = n) +
      0.44 * sin(id / 6) - 0.48 * group
  ))
  binary <- as.integer(
    ((id * 7L + group * 3L) %% 13L) < (5L + group)
  )
  education_code <- (id * 5L + group * 2L) %% 12L
  education <- factor(
    ifelse(
      education_code < 4,
      "A",
      ifelse(education_code < 8, "B", "C")
    ),
    levels = c("A", "B", "C")
  )
  eta <- -0.32 + 0.61 * x_cont + 0.39 * binary +
    0.31 * (education == "B") + 0.67 * (education == "C") +
    0.43 * group
  weight <- 0.40 + ((id * 19L) %% 31L) / 8 +
    0.16 * (x_cont - min(x_cont))
  data.frame(
    id = id,
    group = group,
    y_binary = rbinom(n, size = 1L, prob = plogis(eta)),
    x_cont = x_cont,
    binary = binary,
    education = education,
    weight = weight
  )
}

fit_fairlie_phase6 <- function(
    data,
    model = "logit",
    reference = "group0",
    weight_var = NULL,
    reps = 1L,
    seed = 610107L
) {
  f_decomp(
    data = data,
    dep_var = "y_binary",
    group_var = "group",
    group_levels = c(0, 1),
    indep_vars = c("x_cont", "binary", "education"),
    ref_method = reference,
    pooled_anchor = "favored",
    model_type = model,
    randomize_order = FALSE,
    weight_var = weight_var,
    vce_method = "linearized",
    reps = reps,
    seed = seed,
    quiet = TRUE
  )
}

manual_phase6_robust_vcov <- function(
    x,
    y,
    beta,
    weights,
    model,
    observed = TRUE
) {
  eta <- as.vector(x %*% beta)
  mu <- if (model == "logit") {
    stats::plogis(eta)
  } else {
    stats::pnorm(eta)
  }
  mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
  variance <- mu * (1 - mu)
  derivative <- if (model == "logit") {
    variance
  } else {
    stats::dnorm(eta)
  }
  residual <- y - mu
  score_values <- weights * residual * derivative / variance
  hessian_weights <- weights * derivative^2 / variance
  if (model == "probit" && observed) {
    derivative_prime <- -eta * derivative
    hessian_weights <- weights * (
      derivative^2 / variance -
        residual * derivative_prime / variance +
        residual * derivative^2 * (1 - 2 * mu) / variance^2
    )
  }
  bread <- MASS::ginv(crossprod(x, x * hessian_weights))
  covariance <- bread %*% crossprod(x * score_values) %*% bread
  covariance <- (nrow(x) / (nrow(x) - 1)) * covariance
  dimnames(covariance) <- list(colnames(x), colnames(x))
  covariance
}

test_that("linearized group VCE is the Stata-compatible ML sandwich", {
  data <- make_fairlie_phase6_data()
  x <- stats::model.matrix(
    ~x_cont + binary + education,
    data = data
  )

  for (model in c("logit", "probit")) {
    fit <- fit_fairlie_phase6(data, model = model)
    keep <- data$group == 0
    expected <- manual_phase6_robust_vcov(
      x[keep, , drop = FALSE],
      data$y_binary[keep],
      fit$raw$reference_coefficients,
      rep(1, sum(keep)),
      model
    )
    expect_equal(
      fit$raw$reference_vcov,
      expected,
      tolerance = 1e-12
    )
  }
})

test_that("weighted linearized VCE uses pweight scores and is scale invariant", {
  data <- make_fairlie_phase6_data()
  data$rescaled_weight <- 1000 * data$weight
  x <- stats::model.matrix(
    ~x_cont + binary + education,
    data = data
  )
  keep <- data$group == 0

  for (model in c("logit", "probit")) {
    fit <- fit_fairlie_phase6(
      data,
      model = model,
      weight_var = "weight",
      reps = 17L
    )
    rescaled <- fit_fairlie_phase6(
      data,
      model = model,
      weight_var = "rescaled_weight",
      reps = 17L
    )
    expected <- manual_phase6_robust_vcov(
      x[keep, , drop = FALSE],
      data$y_binary[keep],
      fit$raw$reference_coefficients,
      data$weight[keep],
      model
    )
    expect_equal(
      fit$raw$reference_vcov,
      expected,
      tolerance = 1e-11
    )
    expect_equal(
      fit$raw$reference_vcov,
      rescaled$raw$reference_vcov,
      tolerance = 5e-10
    )
  }
})

test_that("pooled reference covariance includes the group-control equation", {
  data <- make_fairlie_phase6_data()
  fit <- fit_fairlie_phase6(
    data,
    reference = "pooled"
  )
  base_x <- stats::model.matrix(
    ~x_cont + binary + education,
    data = data
  )
  group_term <- fit$raw$pooled_anchor$model_term
  x <- cbind(
    base_x[, "(Intercept)", drop = FALSE],
    stats::setNames(data.frame(data$group), group_term),
    base_x[, setdiff(colnames(base_x), "(Intercept)"), drop = FALSE]
  )
  x <- as.matrix(x)
  beta <- fit$raw$reference_model_coefficients[colnames(x)]
  expected <- manual_phase6_robust_vcov(
    x,
    data$y_binary,
    beta,
    rep(1, nrow(data)),
    "logit"
  )

  expect_equal(
    fit$raw$reference_model_vcov[colnames(x), colnames(x)],
    expected,
    tolerance = 1e-12
  )
})

test_that("detailed VCE reproduces Fairlie block gradients", {
  data <- make_fairlie_phase6_data()
  fit <- fit_fairlie_phase6(data)
  x <- stats::model.matrix(
    ~x_cont + binary + education,
    data = data
  )
  x0 <- x[data$group == 0, , drop = FALSE]
  x1 <- x[data$group == 1, , drop = FALSE]
  beta <- fit$raw$reference_coefficients
  covariance <- fit$raw$reference_vcov
  p0 <- stats::plogis(as.vector(x0 %*% beta))
  p1 <- stats::plogis(as.vector(x1 %*% beta))
  left <- x0[order(p0), , drop = FALSE]
  matched1 <- x1[order(p1), , drop = FALSE]
  expected_variance <- numeric(
    length(fit$raw$decomposition_blocks)
  )
  names(expected_variance) <- paste0(
    "Exp_",
    names(fit$raw$decomposition_blocks)
  )

  for (block in seq_along(fit$raw$decomposition_blocks)) {
    right <- left
    columns <- match(
      fit$raw$decomposition_blocks[[block]],
      colnames(x)
    )
    right[, columns] <- matched1[, columns, drop = FALSE]
    eta_left <- as.vector(left %*% beta)
    eta_right <- as.vector(right %*% beta)
    derivative_left <- stats::plogis(eta_left) *
      (1 - stats::plogis(eta_left))
    derivative_right <- stats::plogis(eta_right) *
      (1 - stats::plogis(eta_right))
    gradient <- colMeans(
      derivative_left * left - derivative_right * right
    )
    expected_variance[[block]] <- as.numeric(
      t(gradient) %*% covariance %*% gradient
    )
    left <- right
  }

  expect_equal(
    diag(fit$raw$detailed_vcov),
    expected_variance,
    tolerance = 1e-14
  )
  expect_true(all(
    fit$raw$detailed_vcov[row(fit$raw$detailed_vcov) !=
      col(fit$raw$detailed_vcov)] == 0
  ))
  expect_false(fit$raw$detailed_vcov_complete)
})

test_that("ordinary linearized inference is normal asymptotic", {
  fit <- fit_fairlie_phase6(make_fairlie_phase6_data())
  detail <- fit$results_detailed_explained
  expected_p <- 2 * stats::pnorm(
    abs(detail$Statistic),
    lower.tail = FALSE
  )
  critical <- stats::qnorm(0.975)

  expect_true(is.infinite(fit$summary_stats$DF))
  expect_equal(detail$P_Value, expected_p, tolerance = 0)
  expect_equal(
    detail$Conf_Low,
    detail$Estimate - critical * detail$Std_Error,
    tolerance = 0
  )
  expect_equal(
    detail$Conf_High,
    detail$Estimate + critical * detail$Std_Error,
    tolerance = 0
  )
  expect_identical(
    fit$raw$linearized_vce$coefficient_covariance,
    "Huber-White maximum-likelihood sandwich"
  )
  expect_identical(
    fit$raw$linearized_vce$finite_sample_correction,
    "N/(N-1) by fitted Logit/Probit equation"
  )
  expect_identical(
    fit$raw$linearized_vce$inference_distribution,
    "normal"
  )
  expect_false("robust" %in% names(formals(f_decomp)))
})

test_that("Probit linearization uses observed rather than Fisher bread", {
  data <- make_fairlie_phase6_data()
  fit <- fit_fairlie_phase6(data, model = "probit")
  x <- stats::model.matrix(
    ~x_cont + binary + education,
    data = data
  )
  keep <- data$group == 0
  fisher <- manual_phase6_robust_vcov(
    x[keep, , drop = FALSE],
    data$y_binary[keep],
    fit$raw$reference_coefficients,
    rep(1, sum(keep)),
    "probit",
    observed = FALSE
  )

  expect_gt(max(abs(fit$raw$reference_vcov - fisher)), 1e-8)
  expect_identical(
    fit$raw$linearized_vce$probit_bread,
    "observed information"
  )
})
