#' Decomposition Methods for Health Inequality
#'
#' @description
#' Four regression-based methods for analyzing two-group outcome gaps or
#' socioeconomic health inequality. Each function returns an invisible named
#' list of estimates, model summaries, diagnostics, and unrounded audit data.
#'
#' @section Choosing a method:
#' * [oby_decomp()] partitions a favored-minus-disadvantaged outcome gap into
#'   explained and unexplained components using OLS, LPM, logit, or probit.
#' * [f_decomp()] applies Fairlie's ranked matching to binary outcomes and
#'   reports detailed explained effects. Matching precision is available in
#'   `raw$matching$precision`.
#' * [wvw_decomp()] reports additive determinant contributions to a
#'   concentration-type index, including a residual and, for nonlinear models,
#'   an approximation error.
#' * [ke_decomp()] regresses a rank- or level-dependent composite target on
#'   determinants. Its coefficients are marginal effects, not additive index
#'   contributions. The aggregate index is reported without inference.
#'
#' @section Survey designs and domains:
#' Declare survey weights, strata, and PSUs with `use_svy = TRUE` and the
#' corresponding variable arguments. Unit weights are used when no weight
#' variable is supplied. Without a PSU variable, each design-source row is
#' treated as a PSU. The interface represents one PSU stage under an
#' ultimate-cluster approximation.
#'
#' To analyze a domain or common complete-case sample, pass the full parent
#' dataset and set `analysis_subset` to a logical vector of length `nrow(data)`
#' with no missing values and at least one `TRUE`. This requires survey mode.
#' Valid positive-weight design rows in participating strata are retained,
#' including PSUs with no analytical observations. Model estimates and
#' empirical socioeconomic positions use only analytical records.
#'
#' All four functions accept `vce_method = "linearized"`, `"jackknife"`, or
#' `"bootstrap"`. Fairlie requires jackknife or bootstrap when strata or PSUs
#' are declared. Its linearized detailed covariance is incomplete. KE
#' linearized coefficient inference is conditional on the constructed target.
#' Consult each function's variance section for design assumptions and limits.
#'
#' @section Reproducibility and interpretation:
#' An explicit `seed` makes stochastic calculations reproducible and restores
#' the caller's random-number state. Fairlie's `reps` controls matching
#' iterations; `boot_reps` controls sampling bootstrap replicates. Check
#' matching precision and stability as these counts increase.
#'
#' Specify comparison groups, socioeconomic ordering, outcome bounds, survey
#' design, and variance method before interpreting results. Decomposition
#' components describe statistical associations and are not causal effects.
#' Implementation credits and methodological references appear in each
#' function's help. Third-party notices are installed in `COPYRIGHTS`.
#'
#' @seealso [lunadecomp_example], [oby_decomp()], [f_decomp()],
#'   [wvw_decomp()], [ke_decomp()]
#' @keywords package
"_PACKAGE"
