test_that("wvw_decomp works for standard OLS decomposition", {
  data(lunadecomp_example)

  fit <- wvw_decomp(
    data = lunadecomp_example,
    dep_var = "health_score",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    ses_var = "ses",
    model_type = "ols",
    correction = "standard",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
  expect_s3_class(fit$results_detailed, "tbl_df")
  expect_s3_class(fit$results_grouped, "tbl_df")

  expect_true("Marginal Effects" %in% names(fit$results_detailed))
  expect_true("Elasticity" %in% names(fit$results_detailed))
  expect_true("Concentration Index" %in% names(fit$results_detailed))
  expect_true("% Contribution Explained" %in% names(fit$results_detailed))
  expect_true("% Contribution Total" %in% names(fit$results_detailed))

  # Harmonized overall inference columns (shared with oby_decomp / f_decomp)
  expect_true(all(c("Term", "Estimate", "Std_Error", "Statistic", "P_Value", "Conf_Low", "Conf_High") %in% names(fit$results_overall)))
})

test_that("wvw_decomp works for generalized concentration index", {
  data(lunadecomp_example)

  fit <- wvw_decomp(
    data = lunadecomp_example,
    dep_var = "health_score",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    ses_var = "ses",
    model_type = "ols",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
})

test_that("wvw_decomp works for Erreygers correction", {
  data(lunadecomp_example)

  fit <- wvw_decomp(
    data = lunadecomp_example,
    dep_var = "health_score",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    ses_var = "ses",
    model_type = "ols",
    correction = "erreygers",
    dep_min = 0,
    dep_max = 1,
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
})

test_that("wvw_decomp works for Wagstaff correction", {
  data(lunadecomp_example)

  fit <- wvw_decomp(
    data = lunadecomp_example,
    dep_var = "health_score",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    ses_var = "ses",
    model_type = "ols",
    correction = "wagstaff",
    dep_min = 0,
    dep_max = 1,
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
})

test_that("wvw_decomp works for logit model with binary outcome", {
  data(lunadecomp_example)

  fit <- wvw_decomp(
    data = lunadecomp_example,
    dep_var = "y_binary",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    ses_var = "ses",
    model_type = "logit",
    correction = "standard",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_detailed, "tbl_df")
})

test_that("wvw_decomp rejects invalid model_type", {
  data(lunadecomp_example)

  expect_error(
    wvw_decomp(
      data = lunadecomp_example,
      dep_var = "health_score",
      indep_vars = c("age", "rural"),
      ses_var = "ses",
      model_type = "invalid",
      quiet = TRUE
    ),
    "model_type"
  )
})