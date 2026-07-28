phase9_survey_data <- function() {
  set.seed(920901)
  data <- expand.grid(
    observation = seq_len(4),
    psu_within_stratum = seq_len(4),
    strata = seq_len(4)
  )
  data$psu <- interaction(
    data$strata,
    data$psu_within_stratum,
    drop = TRUE
  )
  data$id <- seq_len(nrow(data))
  data$group <- as.integer(
    (data$observation + data$psu_within_stratum + data$strata) %% 2
  )
  data$x1 <- stats::rnorm(nrow(data))
  data$x2 <- as.integer(
    (data$observation + 2 * data$psu_within_stratum) %% 3 == 0
  )
  eta <- -0.5 + 0.7 * data$x1 - 0.35 * data$x2 +
    0.45 * data$group + 0.08 * data$strata
  data$y <- stats::rbinom(nrow(data), 1, stats::plogis(eta))
  data$weight <- 0.8 + 0.15 * data$observation +
    0.1 * data$psu_within_stratum + 0.05 * data$strata
  data
}

fit_phase9_survey <- function(method, boot_reps = 12) {
  f_decomp(
    data = phase9_survey_data(),
    dep_var = "y",
    group_var = "group",
    group_levels = c(0, 1),
    indep_vars = c("x1", "x2"),
    ref_method = "group0",
    model_type = "logit",
    randomize_order = FALSE,
    weight_var = "weight",
    psu_var = "psu",
    strata_var = "strata",
    vce_method = method,
    reps = 3,
    boot_reps = boot_reps,
    seed = 920902,
    quiet = TRUE
  )
}

test_that("declared survey variables activate survey mode", {
  fit <- fit_phase9_survey("jackknife")

  expect_true(fit$raw$survey_mode)
  expect_true(fit$raw$complex_survey_design)
  expect_identical(
    fit$raw$matching$sampling_method,
    "ranked_pps_with_replacement"
  )
  expect_equal(fit$summary_stats$Strata, 4)
  expect_equal(fit$summary_stats$PSUs, 16)
  expect_equal(fit$summary_stats$DF, 12)
  expect_identical(
    fit$raw$survey_design_scope$variance_structure,
    "ultimate cluster"
  )
  expect_false(
    fit$raw$survey_design_scope$first_stage_pps_resampling
  )
})

test_that("complex survey linearization is rejected rather than understated", {
  data <- phase9_survey_data()
  expect_error(
    f_decomp(
      data = data,
      dep_var = "y",
      group_var = "group",
      indep_vars = c("x1", "x2"),
      psu_var = "psu",
      strata_var = "strata",
      vce_method = "linearized",
      quiet = TRUE
    ),
    "complete Taylor-linearized VCE is not available"
  )
})

test_that("group models retain the full design before domain selection", {
  data <- phase9_survey_data()
  fit <- fit_phase9_survey("jackknife")
  design <- survey::svydesign(
    ids = ~psu,
    strata = ~strata,
    weights = ~weight,
    data = data,
    nest = TRUE
  )
  domain_design <- design[data$group == 0, ]
  manual_model <- survey::svyglm(
    y ~ x1 + x2,
    design = domain_design,
    family = stats::quasibinomial(link = "logit"),
    control = stats::glm.control(epsilon = 1e-16, maxit = 100)
  )

  expect_equal(
    unname(fit$raw$group_coefficients$group_0),
    unname(stats::coef(manual_model)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(fit$raw$group_vcov$group_0),
    unname(stats::vcov(manual_model)),
    tolerance = 1e-9
  )
})

test_that("stratified survey jackknife uses JKn centering and scale", {
  fit <- fit_phase9_survey("jackknife")
  replication <- fit$raw$replication

  expect_identical(replication$survey_replicate_type, "JKn")
  expect_identical(
    replication$center_method,
    "stratum-specific replicate mean"
  )
  expect_true(is.matrix(replication$center))
  expect_equal(replication$scale, 1)
  expect_equal(replication$failed_replicates, 0)

  deviations <- replication$replicates - replication$center
  manual_vcov <- replication$scale * crossprod(
    deviations * sqrt(replication$rscales)
  )
  expect_equal(
    unname(replication$vcov),
    unname(manual_vcov),
    tolerance = 1e-12
  )
})

test_that("survey bootstrap uses Rao-Wu rescaled PSU weights", {
  bootstrap_reps <- 12
  data <- phase9_survey_data()
  fit <- fit_phase9_survey("bootstrap", bootstrap_reps)
  replication <- fit$raw$replication

  expect_identical(
    replication$survey_replicate_type,
    "Rao-Wu rescaled bootstrap"
  )
  expect_equal(replication$scale, 1 / bootstrap_reps)
  expect_equal(
    replication$generator_scale,
    1 / (bootstrap_reps - 1)
  )
  expect_identical(replication$center_method, "replicate mean")

  factors <- sweep(
    replication$replicate_weights,
    1,
    data$weight,
    "/"
  )
  for (rep_index in seq_len(ncol(factors))) {
    for (stratum in unique(data$strata)) {
      rows <- data$strata == stratum
      psu_factors <- vapply(
        split(factors[rows, rep_index], as.character(data$psu[rows])),
        function(values) unique(values)[[1]],
        numeric(1)
      )
      # Four sampled PSUs imply three draws, rescaled by 4/3.
      expect_equal(sum(psu_factors), 4, tolerance = 1e-12)
      expect_equal(
        psu_factors / (4 / 3),
        round(psu_factors / (4 / 3)),
        tolerance = 1e-12
      )
    }
  }
})

test_that("Rao-Wu bootstrap rejects singleton PSU strata explicitly", {
  data <- phase9_survey_data()
  data <- data[
    data$strata != 1 | data$psu_within_stratum == 1,
    ,
    drop = FALSE
  ]

  expect_error(
    f_decomp(
      data = data,
      dep_var = "y",
      group_var = "group",
      indep_vars = c("x1", "x2"),
      weight_var = "weight",
      psu_var = "psu",
      strata_var = "strata",
      vce_method = "bootstrap",
      boot_reps = 8,
      reps = 2,
      seed = 920904,
      quiet = TRUE
    ),
    "requires at least two PSUs in every stratum"
  )
})

test_that("survey bootstrap can treat singleton PSUs as certainty with svrep", {
  bootstrap_reps <- 12
  data <- phase9_survey_data()
  data <- data[
    data$strata != 1 | data$psu_within_stratum == 1,
    ,
    drop = FALSE
  ]

  fit <- f_decomp(
    data = data,
    dep_var = "y",
    group_var = "group",
    group_levels = c(0, 1),
    indep_vars = c("x1", "x2"),
    ref_method = "group0",
    model_type = "logit",
    randomize_order = FALSE,
    weight_var = "weight",
    psu_var = "psu",
    strata_var = "strata",
    vce_method = "bootstrap",
    reps = 3,
    boot_reps = bootstrap_reps,
    bootstrap_singleton = "certainty",
    seed = 920905,
    quiet = TRUE
  )

  replication <- fit$raw$replication
  expect_identical(
    replication$engine,
    "svrep::as_bootstrap_design"
  )
  expect_identical(
    replication$survey_replicate_type,
    "Rao-Wu-Yue-Beaumont bootstrap with singleton certainty"
  )
  expect_true(replication$bootstrap_singleton$applied)
  expect_equal(replication$bootstrap_singleton$n_singleton_strata, 1)
  expect_identical(
    replication$bootstrap_singleton$first_stage_variance,
    "zero for singleton strata"
  )
  expect_equal(replication$scale, 1 / bootstrap_reps)
  expect_equal(replication$generator_scale, 1 / bootstrap_reps)

  factors <- sweep(
    replication$replicate_weights,
    1,
    data$weight,
    "/"
  )
  singleton_rows <- data$strata == 1
  expect_equal(
    factors[singleton_rows, , drop = FALSE],
    matrix(
      1,
      nrow = sum(singleton_rows),
      ncol = bootstrap_reps
    ),
    tolerance = 1e-12
  )
  for (stratum in setdiff(unique(data$strata), 1)) {
    rows <- data$strata == stratum
    psu_factors <- vapply(
      split(factors[rows, 1], as.character(data$psu[rows])),
      function(values) unique(values)[[1]],
      numeric(1)
    )
    expect_equal(sum(psu_factors), 4, tolerance = 1e-12)
    expect_equal(
      psu_factors / (4 / 3),
      round(psu_factors / (4 / 3)),
      tolerance = 1e-12
    )
  }
  expect_match(
    fit$diagnostics$bootstrap_singleton_control,
    "lonely_psu does not control Rao-Wu bootstrap"
  )
  expect_match(
    fit$diagnostics$bootstrap_singleton_assumption,
    "zero first-stage variance"
  )
})

test_that("bootstrap_singleton is validated independently of lonely_psu", {
  expect_error(
    f_decomp(
      data = phase9_survey_data(),
      dep_var = "y",
      group_var = "group",
      indep_vars = c("x1", "x2"),
      vce_method = "bootstrap",
      bootstrap_singleton = "average",
      quiet = TRUE
    ),
    "bootstrap_singleton must be one of"
  )
})
