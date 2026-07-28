test_that("ke_decomp works for rank-dependent Erreygers index", {
  data(lunadecomp_example)

  fit <- ke_decomp(
    data = lunadecomp_example,
    dep_var = "health_score",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    ses_var = "ses",
    index_type = "rank",
    correction = "erreygers",
    dep_min = 0,
    dep_max = 1,
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
  expect_s3_class(fit$results_anova, "tbl_df")
  expect_s3_class(fit$results_detailed, "tbl_df")

  # results_overall harmonized with the sibling functions (shared inferential
  # spine), plus the KE-specific health bounds.
  expect_true(all(c("Term", "Estimate", "Std_Error", "Statistic", "P_Value", "Conf_Low", "Conf_High") %in% names(fit$results_overall)))
  expect_true(all(c("Min_Health", "Max_Health") %in% names(fit$results_overall)))
  expect_true("Logworth" %in% names(fit$results_anova))
  # results_detailed keeps 'Estimate' (marginal effect) plus KE-specific Prob_F/Logworth.
  expect_true(all(c("Estimate", "Statistic", "P_Value", "Prob_F", "Logworth") %in% names(fit$results_detailed)))
})

test_that("ke_decomp works for rank-dependent standard index", {
  data(lunadecomp_example)

  fit <- ke_decomp(
    data = lunadecomp_example,
    dep_var = "health_score",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    ses_var = "ses",
    index_type = "rank",
    correction = "standard",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
  expect_s3_class(fit$results_detailed, "tbl_df")
})

test_that("ke_decomp works for level-dependent generalized index", {
  data(lunadecomp_example)

  dat <- lunadecomp_example
  dat$income <- exp(dat$ses + 10)

  fit <- ke_decomp(
    data = dat,
    dep_var = "health_score",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    ses_var = "income",
    index_type = "level",
    correction = "generalized",
    vce_method = "linearized",
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
  expect_s3_class(fit$results_anova, "tbl_df")
  expect_s3_class(fit$results_detailed, "tbl_df")
})

test_that("ke_decomp rejects invalid index_type", {
  data(lunadecomp_example)

  expect_error(
    ke_decomp(
      data = lunadecomp_example,
      dep_var = "health_score",
      indep_vars = c("age", "rural"),
      ses_var = "ses",
      index_type = "invalid",
      quiet = TRUE
    ),
    "index_type"
  )
})

test_that("ke_decomp rejects missing socioeconomic variable for rank index", {
  data(lunadecomp_example)

  expect_error(
    ke_decomp(
      data = lunadecomp_example,
      dep_var = "health_score",
      indep_vars = c("age", "rural"),
      index_type = "rank",
      quiet = TRUE
    ),
    "rank"
  )
})