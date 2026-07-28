make_fairlie_phase8_data <- function() {
  n_per_group <- 40L
  id <- seq_len(2L * n_per_group)
  group <- rep(0:1, each = n_per_group)
  x1 <- as.numeric(scale(
    seq(-2.4, 2.4, length.out = length(id)) +
      0.5 * sin(id / 4) - 0.35 * group
  ))
  x2 <- as.integer(((id * 7L + group * 3L) %% 17L) < 8L)
  x3 <- as.integer(((id * 5L + group) %% 13L) < 6L)
  y <- as.integer(
    ((id * 11L + group * 5L + x2 * 2L) %% 23L) <
      (10L + group + x3)
  )
  data.frame(
    id = id,
    group = group,
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
}

fit_fairlie_phase8 <- function(
    data,
    method,
    model = "logit",
    seed = 88101L
) {
  f_decomp(
    data = data,
    dep_var = "y",
    group_var = "group",
    indep_vars = c("x1", "x2", "x3"),
    ref_method = "group0",
    model_type = model,
    randomize_order = FALSE,
    vce_method = method,
    reps = 5L,
    boot_reps = 12L,
    seed = seed,
    quiet = TRUE
  )
}

test_that("ordinary jackknife deletes exactly one analytic row", {
  data <- make_fairlie_phase8_data()
  fit <- fit_fairlie_phase8(data, "jackknife")
  replication <- fit$raw$replication
  weights <- replication$replicate_weights
  n <- nrow(data)

  expect_equal(dim(weights), c(n, n))
  expect_equal(colSums(weights == 0), rep(1, n))
  expect_equal(
    weights[weights > 0],
    rep(n / (n - 1), n * (n - 1))
  )
  expect_equal(replication$scale, (n - 1) / n)
  expect_true(replication$ordinary_jackknife)
  expect_identical(replication$engine, "resample::jackknife")
  expect_identical(replication$simulation, "delete-one")
  expect_identical(
    replication$statistic_type,
    "retained row indices"
  )
  expect_identical(
    replication$fit_weight_type,
    "JK1 replicate factor"
  )
  expect_equal(
    unname(replication$replicates[1, "Group_0"]),
    mean(data$y[data$group == 0 & data$id != 1]),
    tolerance = 0
  )
})

test_that("ordinary bootstrap expands integer multiplicities", {
  data <- make_fairlie_phase8_data()
  fit <- fit_fairlie_phase8(data, "bootstrap")
  replication <- fit$raw$replication
  weights <- replication$replicate_weights

  expect_equal(dim(weights), c(nrow(data), 12))
  expect_true(all(weights == round(weights)))
  expect_equal(colSums(weights), rep(nrow(data), 12))
  expect_equal(replication$scale, 1 / 11)
  expect_true(replication$ordinary_bootstrap)
  expect_identical(replication$engine, "boot::boot")
  expect_identical(replication$simulation, "ordinary")
  expect_identical(replication$statistic_type, "indices")
  expect_identical(replication$fit_weight_type, "frequency")

  first_weights <- weights[, 1]
  expected_group0 <- weighted.mean(
    data$y[data$group == 0],
    first_weights[data$group == 0]
  )
  expect_equal(
    unname(replication$replicates[1, "Group_0"]),
    expected_group0,
    tolerance = 0
  )
})

test_that("replicate VCE uses replicate means and the complete vector", {
  data <- make_fairlie_phase8_data()
  for (method in c("jackknife", "bootstrap")) {
    fit <- fit_fairlie_phase8(data, method)
    replication <- fit$raw$replication
    replicate_estimates <- replication$replicates
    manual_vcov <- replication$scale * crossprod(sweep(
      replicate_estimates,
      2,
      colMeans(replicate_estimates),
      "-"
    ))

    expect_identical(
      replication$center_method,
      "replicate mean"
    )
    expect_equal(
      unname(replication$center),
      unname(colMeans(replicate_estimates)),
      tolerance = 0
    )
    expect_equal(
      fit$raw$vcov,
      manual_vcov,
      tolerance = 1e-15
    )
    expect_false(anyNA(fit$raw$vcov))
    expect_true(fit$raw$overall_vcov_complete)
    expect_true(fit$raw$detailed_vcov_complete)
  }
})

test_that("every sampling replicate reuses all Fairlie iterations", {
  data <- make_fairlie_phase8_data()
  for (method in c("jackknife", "bootstrap")) {
    fit <- fit_fairlie_phase8(data, method)
    replication <- fit$raw$replication

    expect_equal(replication$matching_reps_per_replicate, 5)
    expect_true(replication$common_random_numbers)
    expect_equal(replication$matching_seed, 88101)
    expect_equal(replication$failed_replicates, 0)
  }
})

test_that("ordinary bootstrap is exactly reproducible from seed", {
  data <- make_fairlie_phase8_data()
  first <- fit_fairlie_phase8(data, "bootstrap", model = "probit")
  second <- fit_fairlie_phase8(data, "bootstrap", model = "probit")

  expect_equal(
    first$raw$replication$replicate_weights,
    second$raw$replication$replicate_weights,
    tolerance = 0
  )
  expect_equal(
    first$raw$replication$replicates,
    second$raw$replication$replicates,
    tolerance = 0
  )
  expect_equal(first$raw$vcov, second$raw$vcov, tolerance = 0)
})

test_that("ordinary jackknife preserves the unseeded matching RNG stream", {
  data <- make_fairlie_phase8_data()
  set.seed(88102)
  expected_matching_seed <- sample.int(.Machine$integer.max, 1)
  expected_rng_state <- .Random.seed

  set.seed(88102)
  fit <- fit_fairlie_phase8(
    data,
    "jackknife",
    seed = NULL
  )

  expect_identical(
    fit$raw$replication$matching_seed,
    expected_matching_seed
  )
  expect_identical(.Random.seed, expected_rng_state)
})
