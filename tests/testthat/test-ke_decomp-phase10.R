.phase10_ke_data <- function(n = 30L) {
  id <- seq_len(n)
  data.frame(
    id = id,
    health = 0.15 + 0.65 * stats::plogis(
      -1 + 0.45 * sin(id / 4) + 0.025 * id
    ),
    ses = rep(seq_len(n / 2L), each = 2L),
    x1 = sin(id / 5) + id / 60,
    x2 = as.integer(id %% 3L == 0L),
    weight = 0.7 + (id %% 7L) / 5,
    strata = rep(seq_len(3L), each = n / 3L),
    psu = rep(rep(seq_len(5L), each = 2L), 3L)
  )
}

.phase10_ke_fit <- function(
    data = .phase10_ke_data(),
    index_type = "rank",
    correction = "generalized",
    method = "linearized",
    survey = FALSE,
    ...
) {
  suppressWarnings(ke_decomp(
    data = data,
    dep_var = "health",
    indep_vars = c("x1", "x2"),
    ses_var = "ses",
    index_type = index_type,
    correction = correction,
    use_svy = survey,
    weight_var = if (survey) "weight" else NULL,
    strata_var = if (survey) "strata" else NULL,
    psu_var = if (survey) "psu" else NULL,
    vce_method = method,
    boot_reps = 12L,
    seed = 202610L,
    quiet = TRUE,
    ...
  ))
}

.phase10_manual_fixed_replicates <- function(fit) {
  replication <- fit$raw$replication
  weights <- replication$replicate_weights
  analytic <- fit$raw$sample$analytic_data
  matrix <- fit$raw$model$matrix
  health <- as.numeric(analytic$health)
  position <- fit$raw$target$position
  index_type <- fit$raw$settings$index_type

  estimates <- vapply(
    seq_len(ncol(weights)),
    function(column) {
      current_weight <- weights[, column]
      keep <- current_weight > 0
      mean_health <- sum(current_weight * health) /
        sum(current_weight)
      target <- if (index_type == "rank") {
        2 * position * health - mean_health
      } else {
        position * health - mean_health
      }
      coefficients <- stats::lm.wfit(
        x = matrix[keep, , drop = FALSE],
        y = target[keep],
        w = current_weight[keep]
      )$coefficients
      coefficients[is.na(coefficients)] <- 0
      coefficients
    },
    numeric(ncol(matrix))
  )
  t(estimates)
}

test_that("KE phase 10 validates scalar API arguments", {
  data <- .phase10_ke_data()
  base <- function(...) {
    ke_decomp(
      data, "health", c("x1", "x2"),
      ses_var = "ses",
      correction = "generalized",
      quiet = TRUE,
      ...
    )
  }

  expect_error(base(level = 0), "level")
  expect_error(base(level = 1), "level")
  expect_error(base(level = NA_real_), "level")
  expect_error(base(level = c(0.9, 0.95)), "level")
  expect_error(base(seed = 1.5), "seed")
  expect_error(base(seed = -1), "seed")
  expect_error(base(use_svy = NA), "use_svy")
  expect_error(base(relax = NA), "relax")
  expect_error(base(quiet = 1), "quiet")
  expect_error(
    base(vce_method = c("linearized", "bootstrap")),
    "vce_method"
  )
  expect_error(base(correction = NA_character_), "correction")
  expect_error(base(index_type = character()), "index_type")
  expect_error(base(lonely_psu = "invented"), "lonely_psu")
  expect_error(
    base(bootstrap_singleton = "average"),
    "bootstrap_singleton"
  )
  expect_error(base(boot_reps = "20", vce_method = "bootstrap"),
               "boot_reps")
  expect_error(
    ke_decomp(
      data, NULL, c("x1", "x2"),
      ses_var = "ses",
      correction = "generalized",
      quiet = TRUE
    ),
    "dep_var"
  )
  expect_error(base(weight_var = c("weight", "weight")), "weight_var")
  expect_error(
    base(precalc_level_var = "ses"),
    "precalc_level_var"
  )
  expect_error(
    base(
      index_type = "level",
      precalc_rank_var = "ses"
    ),
    "precalc_rank_var"
  )
  expect_error(
    ke_decomp(
      data, "health", "absent",
      ses_var = "ses",
      correction = "generalized",
      quiet = TRUE
    ),
    "absent from data"
  )
  expect_error(
    ke_decomp(
      numeric(), "health", "x1",
      ses_var = "ses", correction = "generalized"
    ),
    "data must be"
  )
  expect_error(
    ke_decomp(
      data.frame(), "health", "x1",
      ses_var = "ses", correction = "generalized"
    ),
    "at least one observation"
  )
})

test_that("KE phase 10 parses numeric text weights and rejects invalid weights", {
  data <- .phase10_ke_data()
  expected <- data$weight
  data$weight <- format(data$weight, scientific = FALSE, trim = TRUE)
  fit <- .phase10_ke_fit(data, survey = TRUE)

  expect_equal(
    unname(fit$raw$design$weights),
    expected[fit$raw$sample$analytic_source_rows],
    tolerance = 2e-15
  )

  for (invalid in list(0, -1, Inf, "not-a-weight")) {
    broken <- .phase10_ke_data()
    broken$weight <- as.character(broken$weight)
    broken$weight[1] <- invalid
    expect_error(
      .phase10_ke_fit(broken, survey = TRUE),
      "weight_var|sampling weights"
    )
  }
})

test_that("KE phase 10 drops missing required rows but rejects infinities", {
  missing <- .phase10_ke_data()
  missing$health[2] <- NA_real_
  missing$ses[4] <- NA_real_
  missing$x1[6] <- NA_real_
  fit <- .phase10_ke_fit(missing)

  expect_equal(fit$summary_stats$N, nrow(missing) - 3L)
  expect_identical(
    fit$raw$sample$excluded_source_rows$missing_required_values,
    c(2L, 4L, 6L)
  )

  for (variable in c("health", "ses", "x1")) {
    broken <- .phase10_ke_data()
    broken[[variable]][1] <- Inf
    expect_error(
      .phase10_ke_fit(broken),
      "finite"
    )
  }

  broken_outcome <- .phase10_ke_data()
  broken_outcome$health <- as.character(broken_outcome$health)
  broken_outcome$health[1] <- "invalid"
  expect_error(.phase10_ke_fit(broken_outcome), "finite numeric")
})

test_that("KE phase 10 handles total ties and constant outcomes explicitly", {
  tied <- .phase10_ke_data()
  tied$ses <- 1
  tied$health <- 0.4
  fit <- .phase10_ke_fit(tied)

  expect_equal(fit$results_overall$Estimate, 0, tolerance = 1e-15)
  expect_equal(
    unname(fit$raw$target$fractional_rank),
    rep(0.5, nrow(tied)),
    tolerance = 1e-15
  )
  expect_equal(fit$raw$target$target, rep(0, nrow(tied)),
               tolerance = 1e-15)
  expect_true(is.na(fit$model_metrics$R_Squared))

  relative <- tied
  relative$health <- 0
  expect_error(
    .phase10_ke_fit(relative, index_type = "rank",
                    correction = "standard"),
    "strictly positive mean"
  )
  relative$health <- 0.4
  relative$ses <- -1
  expect_error(
    .phase10_ke_fit(relative, index_type = "level"),
    "strictly positive mean"
  )
})

test_that("KE phase 10 reports degenerate factors and preserves collinearity", {
  factor_data <- .phase10_ke_data()
  factor_data$one_level <- factor("only")
  expect_error(
    suppressWarnings(ke_decomp(
      factor_data, "health", "one_level",
      ses_var = "ses", correction = "generalized",
      quiet = TRUE
    )),
    "fewer than two observed levels"
  )

  collinear <- .phase10_ke_data()
  collinear$xcopy <- 2 * collinear$x1
  fit <- suppressWarnings(ke_decomp(
    collinear, "health", c("x1", "xcopy", "x2"),
    ses_var = "ses", correction = "generalized",
    quiet = TRUE
  ))
  expect_true(any(fit$raw$model$aliased))
  expect_true(all(
    fit$raw$model$coefficients[fit$raw$model$aliased] == 0
  ))
  expect_true(all(
    fit$results_detailed$Omitted[fit$raw$model$aliased]
  ))
})

test_that("KE phase 10 precalculated rank reproduces the internal point", {
  data <- .phase10_ke_data()
  internal <- .phase10_ke_fit(data)
  rank_by_row <- numeric(nrow(data))
  rank_by_row[internal$raw$target$source_row] <-
    internal$raw$target$fractional_rank
  data$rank_pre <- rank_by_row
  data$ignored_ses <- NA_real_

  fixed <- suppressWarnings(ke_decomp(
    data, "health", c("x1", "x2"),
    ses_var = "ignored_ses",
    precalc_rank_var = "rank_pre",
    index_type = "rank",
    correction = "generalized",
    quiet = TRUE
  ))

  expect_equal(
    fixed$raw$estimates$overall,
    internal$raw$estimates$overall,
    tolerance = 2e-14
  )
  expect_equal(
    fixed$raw$estimates$detailed,
    internal$raw$estimates$detailed,
    tolerance = 1e-14
  )
  expect_equal(fixed$summary_stats$N, nrow(data))
  expect_true(fixed$raw$position$precalculated)
  expect_match(fixed$raw$position$replicate_semantics, "remain fixed")
})

test_that("KE phase 10 precalculated level reproduces the internal point", {
  data <- .phase10_ke_data()
  internal <- .phase10_ke_fit(data, index_type = "level")
  level_by_row <- numeric(nrow(data))
  level_by_row[internal$raw$target$source_row] <-
    internal$raw$target$relative_level
  data$level_pre <- level_by_row

  fixed <- suppressWarnings(ke_decomp(
    data, "health", c("x1", "x2"),
    precalc_level_var = "level_pre",
    index_type = "level",
    correction = "generalized",
    quiet = TRUE
  ))

  expect_equal(
    fixed$raw$estimates$overall,
    internal$raw$estimates$overall,
    tolerance = 1e-15
  )
  expect_equal(
    fixed$raw$estimates$detailed,
    internal$raw$estimates$detailed,
    tolerance = 1e-14
  )
  expect_true(fixed$raw$position$precalculated)

  invalid <- data
  invalid$level_pre[1] <- -0.01
  expect_error(
    suppressWarnings(ke_decomp(
      invalid, "health", c("x1", "x2"),
      precalc_level_var = "level_pre",
      index_type = "level",
      correction = "generalized",
      quiet = TRUE
    )),
    "nonnegative"
  )
})

test_that("KE phase 10 keeps supplied positions fixed in every replicate", {
  for (index_type in c("rank", "level")) {
    data <- .phase10_ke_data()
    point <- .phase10_ke_fit(data, index_type = index_type)
    position <- point$raw$target$position
    source <- point$raw$target$source_row
    position_by_row <- numeric(nrow(data))
    position_by_row[source] <- position

    if (index_type == "rank") {
      data$position_pre <- position_by_row
    } else {
      data$position_pre <- position_by_row
    }

    for (survey in c(FALSE, TRUE)) {
      for (method in c("jackknife", "bootstrap")) {
        arguments <- list(
          data = data,
          dep_var = "health",
          indep_vars = c("x1", "x2"),
          index_type = index_type,
          correction = "generalized",
          use_svy = survey,
          weight_var = if (survey) "weight" else NULL,
          strata_var = if (survey) "strata" else NULL,
          psu_var = if (survey) "psu" else NULL,
          vce_method = method,
          boot_reps = 12L,
          seed = 202610L,
          quiet = TRUE
        )
        if (index_type == "rank") {
          arguments$precalc_rank_var <- "position_pre"
        } else {
          arguments$precalc_level_var <- "position_pre"
        }
        fit <- suppressWarnings(do.call(ke_decomp, arguments))
        expected <- .phase10_manual_fixed_replicates(fit)

        expect_equal(
          fit$raw$replication$estimates,
          expected,
          tolerance = 2e-14,
          info = paste(index_type, survey, method)
        )
        expect_match(
          fit$raw$position$replicate_semantics,
          "remain fixed"
        )
      }
    }
  }
})

test_that("KE phase 10 validates precalculated-position domains", {
  data <- .phase10_ke_data()
  data$rank_pre <- 0.5
  data$rank_pre[1] <- 1.01
  expect_error(
    suppressWarnings(ke_decomp(
      data, "health", c("x1", "x2"),
      precalc_rank_var = "rank_pre",
      correction = "generalized",
      quiet = TRUE
    )),
    "outside \\[0, 1\\]"
  )
  expect_error(
    suppressWarnings(ke_decomp(
      data, "health", c("x1", "x2"),
      precalc_rank_var = "rank_pre",
      correction = "generalized",
      relax = TRUE,
      quiet = TRUE
    )),
    "outside \\[0, 1\\]"
  )

  data$level_pre <- 0
  expect_error(
    suppressWarnings(ke_decomp(
      data, "health", c("x1", "x2"),
      precalc_level_var = "level_pre",
      index_type = "level",
      correction = "generalized",
      quiet = TRUE
    )),
    "at least one positive"
  )
})
