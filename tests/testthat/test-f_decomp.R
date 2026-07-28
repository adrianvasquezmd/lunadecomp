test_that("f_decomp works for logit Fairlie decomposition", {
  data(lunadecomp_example)

  fit <- f_decomp(
    data = lunadecomp_example,
    dep_var = "y_binary",
    group_var = "group",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    model_type = "logit",
    ref_method = "pooled",
    vce_method = "linearized",
    reps = 10,
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
  expect_s3_class(fit$results_detailed_explained, "tbl_df")
  expect_true(all(c("Group_0", "Group_1", "Difference", "Explained", "Unexplained") %in% fit$results_overall$Term))
  expect_true(all(c("Estimate", "Std_Error", "P_Value", "Conf_Low", "Conf_High") %in% names(fit$results_overall)))
})

test_that("f_decomp works for probit Fairlie decomposition", {
  data(lunadecomp_example)

  fit <- f_decomp(
    data = lunadecomp_example,
    dep_var = "y_binary",
    group_var = "group",
    indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
    model_type = "probit",
    ref_method = "pooled",
    vce_method = "linearized",
    reps = 10,
    quiet = TRUE
  )

  expect_true(is.list(fit))
  expect_s3_class(fit$results_overall, "tbl_df")
  expect_s3_class(fit$results_detailed_explained, "tbl_df")
})

test_that("f_decomp rejects non-binary dependent variables", {
  data(lunadecomp_example)

  expect_error(
    f_decomp(
      data = lunadecomp_example,
      dep_var = "y_continuous",
      group_var = "group",
      indep_vars = c("age", "rural", "insurance"),
      model_type = "logit",
      reps = 5,
      quiet = TRUE
    ),
    "dep_var.*binary|dep_var.*0 and 1|strictly binary"
  )
})

test_that("f_decomp rejects invalid model_type", {
  data(lunadecomp_example)

  expect_error(
    f_decomp(
      data = lunadecomp_example,
      dep_var = "y_binary",
      group_var = "group",
      indep_vars = c("age", "rural"),
      model_type = "ols",
      reps = 5,
      quiet = TRUE
    ),
    "model_type.*logit.*probit|model_type must be"
  )
})
