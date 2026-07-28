#' Fairlie Nonlinear Decomposition for Binary Outcome Gaps
#'
#' @description
#' Decomposes differences in binary outcomes between two groups using the
#' nonlinear decomposition approach proposed by Fairlie for logit and probit
#' models.
#'
#' The function estimates group-specific and reference binary-response models,
#' then partitions the observed outcome gap into an explained component
#' attributable to differences in observed covariates and an unexplained
#' residual component. Detailed explained contributions are calculated by
#' sequentially replacing the covariate distribution of one group with that of
#' the other group. Multi-category predictors are swapped jointly as blocks to
#' preserve valid counterfactual factor profiles.
#'
#' @details
#' `f_decomp()` is designed for binary dependent variables coded as 0 and 1.
#' It supports logit and probit links. The method is useful when the standard
#' linear Oaxaca-Blinder decomposition is not appropriate because the outcome is
#' binary and the analyst wants decomposition results on the response
#' probability scale.
#'
#' Let `group_var` define two comparison groups coded as 0 and 1. The observed
#' gap is computed as:
#'
#' \deqn{
#'   \Delta = Pr(Y = 1 \mid G = 0) - Pr(Y = 1 \mid G = 1).
#' }
#'
#' The explained component is obtained from counterfactual predicted
#' probabilities generated from the selected reference coefficient vector. The
#' unexplained component is computed as the difference between the observed gap
#' and the explained component.
#'
#' The detailed decomposition is path-dependent in nonlinear models because the
#' contribution of each covariate can depend on the order in which covariate
#' distributions are swapped. When `randomize_order = TRUE`, the function
#' randomizes the order of blocks across Monte Carlo iterations and reports the
#' average contribution. This follows Fairlie's recommendation to assess or
#' average over ordering sensitivity.
#'
#' When the two groups have unequal sizes and no sampling weights are supplied,
#' every observation in the smaller group is retained and a simple random
#' sample without replacement is drawn from the larger group in each
#' replication. Both samples are then ordered by their reference probabilities
#' and matched by rank. `raw$matching$detailed_monte_carlo_se` reports the
#' simulation error of each detailed contribution across these replications.
#' When `randomize_order = TRUE`, it also includes the finite-replication
#' variation caused by randomizing the block order. It quantifies algorithmic
#' noise and is distinct from the sampling standard errors reported in the
#' public results table.
#'
#' @section Group coding and the direction of the gap:
#' As in [oby_decomp()], the favored (socially advantaged) group is coded 0 and
#' the disadvantaged group is coded 1, and the reported gap is always
#'
#' \deqn{\Delta = Pr(Y = 1 \mid \text{favored}) - Pr(Y = 1 \mid \text{disadvantaged}).}
#'
#' Which group is "favored" is a substantive choice: whatever level is coded 0
#' becomes the base whose covariate distribution is revalued in the explained
#' component. "Favored" refers to social advantage, not to the level of the
#' outcome, so for an adverse outcome (e.g. mortality) the favored group may
#' yield \eqn{\Delta < 0}. Either code `group_var` as 0/1 (group 0 favored by
#' convention), leave its two natural labels and declare `favored_group`, or use
#' `group_levels = c(favored, disadvantaged)` to select and order two levels
#' from a grouping variable with two or more levels. The recoding and resulting
#' mapping are reported in the console and returned in `summary_stats` and
#' `raw$group_mapping`. Because the detailed Fairlie contributions are
#' path-dependent, prefer `randomize_order = TRUE` and keep the group coding
#' fixed across comparisons.
#'
#' @section Reference coefficient structure:
#' The argument `ref_method` defines the coefficient vector used to construct
#' counterfactual predictions:
#'
#' \itemize{
#'   \item `"pooled"`: coefficients from a pooled model including a translated
#'   group indicator. The point at which that indicator is removed from the
#'   decomposition reference is selected by `pooled_anchor`.
#'   \item `"neumark"`: coefficients from a pooled model excluding the group
#'   indicator.
#'   \item `"cotton"`: weighted average of group-specific coefficients using
#'   group sample weights.
#'   \item `"reimers"`: equally weighted average of group-specific coefficients.
#'   \item `"group0"`: coefficients from the group coded as 0.
#'   \item `"group1"`: coefficients from the group coded as 1.
#' }
#'
#' @section Survey designs:
#' If `use_svy = TRUE`, the function uses the supplied weights, strata, and
#' primary sampling units to estimate survey-adjusted logit or probit models.
#' Supplying any of `weight_var`, `psu_var`, or `strata_var` activates survey
#' mode automatically.
#'
#' A supplied `weight_var` has probability-weight (`pweight`) semantics. Missing
#' weights follow complete-case exclusion; zero-weight observations are removed
#' from the analytic sample; and negative, non-finite, or nonnumeric weights are
#' rejected. Multiplying every weight by the same positive constant does not
#' change the point estimator. Frequency-weight and importance-weight semantics
#' are not inferred from this single argument. Thus the survey VCE is
#' conditional on the positive-weight analytic sample; this intentionally
#' follows `fairlie` and is not the same as retaining structural zero-weight
#' records in a generic survey design.
#'
#' The current survey interface is an ultimate-cluster design: one final
#' sampling weight, one optional first-stage PSU identifier, and one optional
#' stratum identifier. It does not accept lower-stage identifiers, finite
#' population corrections, calibration replicate weights, or a first-stage PPS
#' specification. For such designs, the built-in replication methods should
#' not be treated as a substitute for agency-supplied replicate weights.
#'
#' Weighted detailed matching follows Ben Jann's `fairlie`: each replication
#' draws
#' \eqn{\lfloor(N_0 + N_1)/2\rfloor} observations independently from each
#' probability-ranked group, with replacement and selection probability
#' proportional to the weight. Consequently, weighted detailed contributions
#' contain Monte Carlo error even when the two groups have the same size.
#' Increase `reps` and inspect
#' `raw$matching$detailed_monte_carlo_se` when precise detailed estimates are
#' required.
#'
#' @section Variance estimation:
#' Standard errors can be estimated with analytic linearization,
#' jackknife replicate variance estimation, or bootstrap replicate variance
#' estimation. The linearized option is fast and useful during development;
#' replicate-based approaches are preferable for final inference when feasible.
#'
#' Outside a complex survey, `"linearized"` always uses the Huber--White
#' maximum-likelihood sandwich for the Logit or Probit reference coefficients.
#' This is the documented equivalent of
#' `fairlie, vce(robust)`: each fitted equation uses Stata's
#' \eqn{N/(N-1)} finite-sample correction, and Probit uses the
#' observed-information Hessian rather than the Fisher-information bread. There
#' is no separate `robust` switch. Ordinary linearized tests and confidence
#' intervals use the standard normal distribution.
#'
#' For each detailed block and matching replication, the function constructs
#' the Fairlie delta-method gradient and evaluates
#' \eqn{d_x^\top V_\beta d_x}. The reported variance is the average of these
#' per-replication variances, exactly as in Ben Jann's implementation.
#' Covariances between detailed blocks are set to zero for compatibility and
#' `raw$detailed_vcov_complete` is therefore `FALSE`; joint tests of detailed
#' contributions should not be constructed from this matrix.
#'
#' A Taylor-linearized VCE for the complete Fairlie estimator is not offered
#' when strata or PSUs are declared. The detailed estimator contains ranked
#' matching, possible PPS subsampling, and path averaging; propagating only the
#' coefficient covariance would omit the design variation of the empirical
#' covariate distributions and their covariance with the fitted model.
#' Consequently, complex designs must use `"jackknife"` or `"bootstrap"`.
#' This restriction prevents a coefficient-only approximation from being
#' reported as a complete survey VCE.
#'
#' In ordinary replicate modes, [resample::jackknife()] generates one
#' leave-one-out sample per analytic observation and the covariance uses scale
#' \eqn{(N-1)/N}. [boot::boot()] generates ordinary samples of \eqn{N}
#' observations with replacement and the covariance uses scale
#' \eqn{1/(B-1)}. Both are centered on the replicate mean, matching Stata's
#' default non-MSE prefix convention. Bootstrap multiplicities are expanded
#' into physical rows for the Fairlie matching, and the jackknife-deleted row
#' is physically removed rather than retained with zero weight.
#'
#' Every sampling replicate re-runs all `reps` Fairlie matching iterations.
#' A common matching random-number stream is reused across sampling replicates
#' to keep Monte Carlo matching noise from inflating the sampling VCE. The
#' replicate estimates, weights, centering, scale, and complete covariance
#' matrix are returned in `raw$replication`. If any replicates fail, the
#' default is to stop for a survey design because silently dropping a PSU
#' replicate can invalidate the design variance. With `relax = TRUE`, failed
#' replicates are dropped and the result is explicitly marked exploratory.
#' For `"reimers"` and `"cotton"`, the
#' reference-coefficient covariance combines the group-specific covariances
#' assuming independence between the two group models (their samples are
#' disjoint). Reimers fixes the group-0 coefficient weight at 0.5. Cotton fixes
#' it at the full-sample weighted share of group 0. These reference weights are
#' decomposition settings and remain fixed across jackknife and bootstrap
#' replicates; they are not re-estimated from the composition of each replicate.
#' Consequently, linearized and replicate inference are conditional on the
#' observed Reimers/Cotton reference weight.
#'
#' Survey jackknife uses JK1 without strata and JKn with strata. JKn variance
#' is centered on the mean of the replicates belonging to the affected
#' stratum. Survey bootstrap uses the Rao--Wu rescaled construction: within
#' each stratum it samples \eqn{n_h-1} PSUs with replacement and multiplies
#' their sampling weights by \eqn{n_h/(n_h-1)}. Its final non-MSE variance uses
#' Stata's \eqn{1/B} convention; the generator's native finite-replicate scale
#' is retained in `raw$replication$generator_scale`.
#'
#' `lonely_psu` does not control the construction of Rao--Wu bootstrap
#' replicates. Use `bootstrap_singleton` for that purpose. With
#' `bootstrap_singleton = "fail"` (the default), a singleton-PSU stratum stops
#' survey bootstrap estimation. With `bootstrap_singleton = "certainty"`,
#' [svrep::as_bootstrap_design()] generates Rao--Wu--Yue--Beaumont replicate
#' weights and holds the singleton PSU fixed with replicate factor 1. This
#' assigns zero first-stage variance to that stratum and should therefore be
#' selected only when treating the singleton as a certainty PSU is
#' substantively defensible. Regular strata continue to use Rao--Wu
#' \eqn{n_h-1} resampling. `lonely_psu` can still affect Taylor covariance
#' reported for fitted survey-model coefficients, but it does not determine the
#' Fairlie bootstrap VCE.
#'
#' This Rao--Wu construction assumes simple random sampling of PSUs within
#' strata (with the usual with-replacement/ultimate-cluster variance
#' approximation). It is not a generic bootstrap for first-stage PPS designs.
#'
#' Detailed contributions are not rescaled to the full-sample explained total.
#' This matches Ben Jann's implementation: with unequal groups or weighted
#' matching, their sum targets the mean matched-sample path and can differ from
#' the full-sample explained component because of Monte Carlo sampling.
#'
#' @section Interactions:
#' Formula interactions such as `"x1:x2"` or `"x1 * x2"` are deliberately
#' rejected. In a sequential nonlinear decomposition, replacing a main effect
#' without coherently recomputing every interaction that contains it can create
#' counterfactual rows that do not correspond to any valid covariate profile.
#' There is no unique detailed allocation rule for such a path.
#'
#' A precomputed interaction column cannot be distinguished programmatically
#' from an ordinary numeric predictor. It should therefore be used only when
#' every main effect and precomputed interaction belonging to the interacting
#' system is placed in one `groupings` block. The resulting contribution is for
#' the joint system; the function does not support separate hierarchical
#' attribution within that block.
#'
#' @section Predicted-outcome overlap diagnostic:
#' The function reports the overlap of the two groups' reference-model
#' predicted outcome probabilities. This is a descriptive diagnostic for the
#' amount of extrapolation implicit in rank matching; it is not a propensity
#' score, a causal common-support test, or a trimming rule. Observations are
#' never removed by this diagnostic. The ranges, their intersection, and
#' weighted and unweighted off-overlap proportions are returned in
#' `raw$predicted_outcome_overlap`.
#'
#'
#' @section Percentage contribution columns:
#' The detailed Fairlie table reports explained contributions by variable or
#' factor block. Because the Fairlie procedure decomposes the explained
#' component in detail, each row is reported with two percentage denominators:
#'
#' \itemize{
#'   \item `% Contribution Explained`: the variable contribution divided by the
#'   total explained component.
#'   \item `% Contribution Total`: the variable contribution divided by the
#'   observed total gap.
#' }
#'
#' The summary table also reports the percentage of the total gap attributable
#' to the total explained and total unexplained components.
#'
#' @param data A data frame containing the outcome, group, predictors, and
#' optional survey design variables.
#' @param dep_var Character string. Name of the dependent binary variable. The
#' variable must be coded as 0 and 1.
#' @param group_var Character string. Name of the two-group comparison variable.
#' Coded as 0 and 1 (0 = favored/base/reference, 1 = disadvantaged/comparison),
#' unless `favored_group` is supplied for a two-level variable or
#' `group_levels` explicitly selects and orders two observed levels.
#' @param favored_group Optional. The level of `group_var` that identifies the
#' favored (socially advantaged) group. When supplied, `group_var` may use any
#' two labels: the favored level is recoded to 0 and the other to 1 before
#' estimation, so that the reported gap is
#' `Difference = E(Y | favored) - E(Y | disadvantaged)`. When `NULL` (default),
#' `group_var` must already be coded 0/1 unless `group_levels` is supplied.
#' Retained for backward compatibility; new analyses should prefer
#' `group_levels`.
#' @param group_levels Optional ordered vector of length two:
#' `c(favored, disadvantaged)`. The selected levels are recoded internally as
#' 0 and 1 and all other observed group levels are filtered before estimation.
#' Cannot be combined with `favored_group`.
#' @param indep_vars Character vector. Names of independent variables included
#' in the decomposition model. Categorical predictors should be coded as factors
#' or character variables.
#' @param groupings Optional named list. Each element should contain one or more
#' predictor names or prefixes to aggregate detailed effects into conceptual
#' domains. If empty, variables are grouped automatically by model term, and
#' factor levels are swapped jointly. Every specification must match a modeled
#' term, column, or prefix; a model column cannot belong to more than one block,
#' and columns generated by the same multi-level factor cannot be split across
#' blocks. Unspecified columns are collected in an explicit residual block.
#' @param ref_method Character string. Counterfactual reference structure.
#' Options are `"pooled"`, `"neumark"`, `"cotton"`, `"reimers"`, `"group0"`,
#' and `"group1"`.
#' @param pooled_anchor Character string. Point of the group indicator used for
#' the pooled reference. Options are `"favored"` (0), `"disadvantaged"` (1),
#' and `"centered"` (the weighted mean of the disadvantaged-group indicator).
#' Used only when `ref_method = "pooled"`; the default is `"favored"`, matching
#' `fairlie, pooled(groupvar)` after the favored group is coded 0.
#' @param model_type Character string. Binary-response model. Options are
#' `"logit"` and `"probit"`.
#' @param randomize_order Logical. If `TRUE`, randomizes the order in which
#' variable blocks are swapped across Monte Carlo iterations to reduce
#' path-dependence in detailed contributions. The resulting finite-repetition
#' error is returned in `raw$matching$detailed_monte_carlo_se`.
#' @param use_svy Logical. If `TRUE`, applies complex survey design settings.
#' @param weight_var Optional character string. Name of a nonnegative sampling
#' weight variable, interpreted with probability-weight (`pweight`) semantics.
#' Missing weights are omitted with other incomplete rows, zero weights are
#' excluded, and negative, non-finite, or nonnumeric values are rejected.
#' @param psu_var Optional character string. Name of the primary sampling unit
#' variable used for an ultimate-cluster variance approximation. Lower sampling
#' stages and stage-specific weights are not currently represented.
#' @param strata_var Optional character string. Name of the stratification
#' variable.
#' @param vce_method Character string. Variance estimation method. Options are
#' `"linearized"`, `"jackknife"`, and `"bootstrap"`. Outside a complex survey,
#' `"linearized"` uses the Huber--White maximum-likelihood sandwich. With
#' declared strata or PSUs, use `"jackknife"` or `"bootstrap"`; a complete
#' Taylor linearization of the ranked Fairlie estimator is not available.
#' @param reps Integer. Number of Monte Carlo iterations used for the Fairlie
#' matching and ordering procedure. For unequal unweighted groups,
#' `raw$matching$detailed_monte_carlo_se` can be used to assess whether `reps`
#' is sufficiently large.
#' @param boot_reps Integer. Number of bootstrap replications when
#' `vce_method = "bootstrap"`. Ordinary observation-level samples are generated
#' by [boot::boot()]; survey bootstrap continues to use replicate weights from
#' the declared survey design.
#' @param lonely_psu Character string. Handling of strata with one PSU. Options
#' are `"fail"`, `"remove"`, `"certainty"`, `"adjust"`, and `"average"`.
#' This controls Taylor-linearized survey covariance and supported jackknife
#' conversions; it does not control Rao--Wu bootstrap replicate construction.
#' See `bootstrap_singleton`.
#' @param bootstrap_singleton Character string. Handling of a stratum with one
#' PSU when `vce_method = "bootstrap"` is used with a survey design. Options
#' are `"fail"` (the default) and `"certainty"`. `"fail"` requires at least
#' two PSUs per stratum and stops rather than making an assumption. The latter uses
#' [svrep::as_bootstrap_design()] and fixes the singleton PSU at replicate
#' factor 1, thereby assigning zero first-stage variance to that stratum. It
#' does not implement the Taylor `"adjust"` or `"average"` rules.
#' @param level Numeric. Confidence level for confidence intervals. Default is
#' `0.95`.
#' @param seed Optional integer. Random seed for reproducibility of the Fairlie
#' matching/ordering and replicate resampling. Defaults to `NULL`
#' (non-deterministic); set an explicit value to reproduce results. The prior
#' state of the global random number generator is restored on exit.
#' @param relax Logical. If `TRUE`, selected structural validation problems are
#' stored as diagnostics rather than stopping execution.
#' @param quiet Logical. If `TRUE`, suppresses console output and progress bars.
#'
#' @return
#' A named list with the following components (harmonized with [oby_decomp()]):
#' \describe{
#'   \item{summary_stats}{Sample size, population size, number of strata, number
#'   of PSUs, design degrees of freedom, group-specific sample sizes, and the
#'   labels mapped to the favored (code 0) and disadvantaged (code 1) groups.}
#'   \item{model_metrics}{Model fit statistics for group-specific and reference
#'   models.}
#'   \item{models_coefficients}{Unrounded coefficient tables for the models
#'   actually required by the selected reference. The reference table is always
#'   included; an unused group model is not estimated or returned.}
#'   \item{diagnostics}{Warnings and diagnostic messages generated during
#'   validation, estimation, or variance estimation.}
#'   \item{results_overall}{Overall decomposition table with the same columns as
#'   `oby_decomp()`: group predictions (`Group_0`, `Group_1`), total
#'   `Difference`, `Explained`, and `Unexplained` components, each with standard
#'   error, test statistic, p-value, and confidence interval.}
#'   \item{results_detailed_explained}{Unrounded detailed explained
#'   contributions by variable or grouped factor block, reported without
#'   rescaling, including `% Contribution Explained`, `% Contribution Total`,
#'   standard errors, and confidence intervals.}
#'   \item{results_detailed_unexplained}{`NULL`. The Fairlie method decomposes
#'   the explained component in detail but does not partition the unexplained
#'   component by predictor; this slot is kept for structural parity with
#'   `oby_decomp()`.}
#'   \item{raw}{Unrounded overall and detailed estimates, standard errors and
#'   available covariance matrices; group and reference coefficients; the full
#'   reference-model coefficient vector and covariance; reference predictions;
#'   analytic-sample rows and a row-level exclusion audit in `raw$sample_flow`;
#'   decomposition blocks; fitted-model inventory; matching diagnostics,
#'   including the separate Monte Carlo error of the detailed matching
#'   estimates; replicate estimates, weights, centering, scale, and engine
#'   metadata in `raw$replication`; group mapping; weight metadata in
#'   `raw$weighting`;
#'   reference classification and fixed Reimers/Cotton weights in
#'   `raw$reference_structure`; the interaction contract in
#'   `raw$interaction_handling`; the descriptive predicted-outcome overlap
#'   audit in `raw$predicted_outcome_overlap`;
#'   linearized-VCE conventions in `raw$linearized_vce`; and estimation
#'   settings. Covariance matrices that are incomplete under the current
#'   linearized method retain unavailable elements as `NA` rather than treating
#'   them as zero.}
#' }
#'
#' @examples
#' data(lunadecomp_example)
#'
#' fit_fairlie <- f_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "y_binary",
#'   group_var = "group",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   model_type = "logit",
#'   ref_method = "pooled",
#'   vce_method = "linearized",
#'   reps = 25,
#'   quiet = TRUE
#' )
#'
#' fit_fairlie$results_overall
#' fit_fairlie$results_detailed_explained
#'
#' \donttest{
#' fit_probit <- f_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "y_binary",
#'   group_var = "group",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   model_type = "probit",
#'   ref_method = "pooled",
#'   vce_method = "linearized",
#'   reps = 50,
#'   quiet = TRUE
#' )
#'
#' fit_probit$results_overall
#' }
#'
#' @author
#' Adrian Vasquez-Mejia, MD, MSc \cr
#' Oscar J. Mujica, MD, MPH, PHE, FACE \cr
#' Antonio Sanhueza, MPH, MSc, PhD
#'
#' @references
#' Fairlie, R. W. (2005). An extension of the Blinder-Oaxaca decomposition
#' technique to logit and probit models. \emph{Journal of Economic and Social
#' Measurement}, 30(4), 305--316. \doi{10.3233/JEM-2005-0259}
#'
#' Blinder, A. S. (1973). Wage discrimination: Reduced form and structural
#' estimates. \emph{The Journal of Human Resources}, 8(4), 436.
#' \doi{10.2307/144855}
#'
#' Oaxaca, R. (1973). Male-female wage differentials in urban labor markets.
#' \emph{International Economic Review}, 14(3), 693.
#' \doi{10.2307/2525981}
#'
#' Oaxaca, R. L., & Ransom, M. R. (1994). On discrimination and the
#' decomposition of wage differentials. \emph{Journal of Econometrics}, 61(1),
#' 5--21. \doi{10.1016/0304-4076(94)90074-4}
#'
#' @importFrom dplyr arrange desc
#' @importFrom tibble tibble add_row
#' @importFrom rlang .data
#' @importFrom survey svydesign svyglm as.svrepdesign withReplicates
#' @importFrom stats as.formula coef complete.cases cov.wt glm median model.matrix pt qt quasibinomial runif terms var vcov weighted.mean
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export

f_decomp <- function(
    data, dep_var, group_var, indep_vars, groupings = list(),
    favored_group = NULL,
    ref_method = "pooled", model_type = "logit", randomize_order = TRUE,
    use_svy = FALSE, weight_var = NULL, psu_var = NULL, strata_var = NULL,
    vce_method = "linearized", reps = 1000, boot_reps = 500,
    lonely_psu = "adjust", level = 0.95, seed = NULL, relax = FALSE,
    quiet = FALSE, group_levels = NULL, pooled_anchor = "favored",
    bootstrap_singleton = "fail"
) {

  # C1/F3: restore the caller's RNG state on exit instead of leaving it altered.
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  if (
    !use_svy &&
      (!is.null(weight_var) || !is.null(psu_var) || !is.null(strata_var))
  ) {
    use_svy <- TRUE
  }
  complex_survey_design <- use_svy &&
    (!is.null(psu_var) || !is.null(strata_var))

  # === 1. STRUCTURAL VALIDATION & DATA CONFIGURATION ===
  valid_vce <- c("linearized", "jackknife", "bootstrap")
  valid_ref <- c("pooled", "neumark", "cotton", "reimers", "group0", "group1")
  valid_model <- c("logit", "probit")
  valid_lonely <- c("fail", "remove", "certainty", "adjust", "average")
  valid_bootstrap_singleton <- c("fail", "certainty")
  valid_pooled_anchor <- c("favored", "disadvantaged", "centered")

  if (!(vce_method %in% valid_vce)) {
    stop(
      "Error: vce_method must be one of: ",
      paste(valid_vce, collapse = ", "),
      call. = FALSE
    )
  }

  if (complex_survey_design && identical(vce_method, "linearized")) {
    stop(
      "Error: a complete Taylor-linearized VCE is not available for the ",
      "ranked Fairlie estimator with declared strata or PSUs. Use ",
      "vce_method = \"jackknife\" or \"bootstrap\". A coefficient-only ",
      "delta approximation would omit design variation in the empirical ",
      "covariate distributions.",
      call. = FALSE
    )
  }

  if (!(ref_method %in% valid_ref)) {
    stop(
      "Error: ref_method must be one of: ",
      paste(valid_ref, collapse = ", "),
      call. = FALSE
    )
  }

  if (!(model_type %in% valid_model)) {
    stop(
      "Error: model_type must be one of: ",
      paste(valid_model, collapse = ", "),
      call. = FALSE
    )
  }

  if (
    !is.character(pooled_anchor) ||
      length(pooled_anchor) != 1 ||
      is.na(pooled_anchor) ||
      !(pooled_anchor %in% valid_pooled_anchor)
  ) {
    stop(
      "Error: pooled_anchor must be one of: ",
      paste(valid_pooled_anchor, collapse = ", "),
      call. = FALSE
    )
  }

  if (!(lonely_psu %in% valid_lonely)) {
    stop(
      "Error: lonely_psu must be one of: ",
      paste(valid_lonely, collapse = ", "),
      call. = FALSE
    )
  }

  if (
    !is.character(bootstrap_singleton) ||
      length(bootstrap_singleton) != 1L ||
      is.na(bootstrap_singleton) ||
      !(bootstrap_singleton %in% valid_bootstrap_singleton)
  ) {
    stop(
      "Error: bootstrap_singleton must be one of: ",
      paste(valid_bootstrap_singleton, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.data.frame(data)) {
    stop("Error: data must be a data frame.", call. = FALSE)
  }

  if (!is.character(dep_var) || length(dep_var) != 1) {
    stop("Error: dep_var must be a single character string.", call. = FALSE)
  }

  if (!is.character(group_var) || length(group_var) != 1) {
    stop("Error: group_var must be a single character string.", call. = FALSE)
  }

  if (!is.character(indep_vars) || length(indep_vars) < 1) {
    stop("Error: indep_vars must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.list(groupings)) {
    stop("Error: groupings must be a named list.", call. = FALSE)
  }

  if (length(groupings) > 0) {
    grouping_names <- names(groupings)
    if (
      is.null(grouping_names) ||
        any(!nzchar(grouping_names)) ||
        anyDuplicated(grouping_names)
    ) {
      stop(
        "Error: every element of groupings must have a unique, non-empty name.",
        call. = FALSE
      )
    }
    invalid_grouping_values <- vapply(
      groupings,
      function(value) {
        !is.character(value) ||
          length(value) < 1 ||
          anyNA(value) ||
          any(!nzchar(value))
      },
      logical(1)
    )
    if (any(invalid_grouping_values)) {
      stop(
        "Error: each groupings element must be a non-empty character vector ",
        "of model terms, columns, or prefixes.",
        call. = FALSE
      )
    }
  }

  if (!is.null(group_levels) && !is.null(favored_group)) {
    stop(
      "Error: use either group_levels or favored_group, not both.",
      call. = FALSE
    )
  }

  if (
    !is.null(group_levels) &&
      (!is.atomic(group_levels) || length(group_levels) != 2 ||
        anyNA(group_levels))
  ) {
    stop(
      "Error: group_levels must be an ordered vector of length two: ",
      "c(favored, disadvantaged).",
      call. = FALSE
    )
  }

  if (
    !is.null(group_levels) &&
      length(unique(as.character(group_levels))) != 2
  ) {
    stop(
      "Error: group_levels must contain two distinct values.",
      call. = FALSE
    )
  }

  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("Error: level must be a numeric value between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(reps) || length(reps) != 1 || reps < 1) {
    stop("Error: reps must be a numeric value greater than or equal to 1.", call. = FALSE)
  }

  if (!is.numeric(boot_reps) || length(boot_reps) != 1 || boot_reps < 2) {
    stop("Error: boot_reps must be a numeric value greater than or equal to 2.", call. = FALSE)
  }

  old_lonely_psu <- getOption("survey.lonely.psu")
  options(survey.lonely.psu = lonely_psu)
  on.exit(options(survey.lonely.psu = old_lonely_psu), add = TRUE)

  requested_formula <- stats::as.formula(
    paste("~", paste(indep_vars, collapse = " + "))
  )
  if (any(attr(stats::terms(requested_formula), "order") > 1)) {
    stop(
      paste0(
        "CRITICAL: Interaction terms detected. Nonlinear path-dependent ",
        "decomposition with interactions is theoretically ambiguous and ",
        "disabled. If an interaction is substantively required, precompute it ",
        "and exchange the complete interacting system in one groupings block."
      ),
      call. = FALSE
    )
  }

  req_vars <- unique(c(dep_var, group_var, indep_vars, weight_var, psu_var, strata_var))
  req_vars <- req_vars[!is.na(req_vars) & nzchar(req_vars)]

  missing_vars <- setdiff(req_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "Error: the following variables are missing from data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  df <- data[, req_vars, drop = FALSE]
  df$.f_source_row <- seq_len(nrow(data))
  N_initial <- nrow(df)
  complete_rows <- stats::complete.cases(df)
  missing_source_rows <- df$.f_source_row[!complete_rows]
  df <- df[complete_rows, , drop = FALSE]
  n_missing_dropped <- N_initial - nrow(df)

  for (v in indep_vars) {
    if (is.character(df[[v]])) df[[v]] <- as.factor(df[[v]])
  }
  df <- droplevels(df)

  df[[dep_var]] <- as.numeric(as.character(df[[dep_var]]))
  if (any(is.na(df[[dep_var]]))) {
    stop("CRITICAL DATA ERROR: Numeric conversion introduced NA values in: ", dep_var, ". Check variable coding.")
  }

  n_zero_weight_dropped <- 0L
  if (!is.null(weight_var)) {
    numeric_weight <- suppressWarnings(
      as.numeric(as.character(df[[weight_var]]))
    )
    if (anyNA(numeric_weight)) {
      stop(
        "Error: weight_var must contain numeric weights. Numeric conversion ",
        "introduced missing values.",
        call. = FALSE
      )
    }
    if (any(!is.finite(numeric_weight))) {
      stop(
        "Error: weight_var must contain only finite weights.",
        call. = FALSE
      )
    }
    if (any(numeric_weight < 0)) {
      stop(
        "Error: weight_var cannot contain negative weights.",
        call. = FALSE
      )
    }
    df[[weight_var]] <- numeric_weight
  }

  # -----------------------------------------------------------------------------
  # GROUP CODING & FAVORED / DISADVANTAGED ASSIGNMENT
  # The reported gap is always Difference = E(Y | favored) - E(Y | disadvantaged),
  # with favored coded 0 (base / reference) and disadvantaged coded 1.
  # -----------------------------------------------------------------------------
  grp_levels_raw <- unique(as.character(df[[group_var]]))
  grp_levels_raw <- grp_levels_raw[!is.na(grp_levels_raw)]
  available_group_levels <- grp_levels_raw
  excluded_group_levels <- character(0)
  n_group_filtered <- 0L
  group_filtered_source_rows <- integer(0)

  if (!is.null(group_levels)) {
    selected_group_levels <- as.character(group_levels)
    missing_group_levels <- setdiff(
      selected_group_levels,
      available_group_levels
    )
    if (length(missing_group_levels) > 0) {
      stop(
        "Error: group_levels contains value(s) not observed in group_var: ",
        paste(missing_group_levels, collapse = ", "),
        ". Available levels are: ",
        paste(available_group_levels, collapse = ", "),
        call. = FALSE
      )
    }
    favored_chr <- selected_group_levels[[1]]
    disadvantaged_chr <- selected_group_levels[[2]]
    keep_group <- as.character(df[[group_var]]) %in%
      selected_group_levels
    n_group_filtered <- sum(!keep_group)
    group_filtered_source_rows <- df$.f_source_row[!keep_group]
    excluded_group_levels <- setdiff(
      available_group_levels,
      selected_group_levels
    )
    df <- droplevels(df[keep_group, , drop = FALSE])
    df[[group_var]] <- ifelse(
      as.character(df[[group_var]]) == favored_chr,
      0,
      1
    )
    favored_label <- favored_chr
    disadvantaged_label <- disadvantaged_chr
    group_selection_mode <- "group_levels"
    if (!quiet) {
      message(sprintf(
        paste0(
          "Note: group_levels selected '%s' as favored (coded 0) and ",
          "'%s' as disadvantaged (coded 1); %d observation(s) from ",
          "other group levels were filtered. Gap = E(Y|favored) - ",
          "E(Y|disadvantaged)."
        ),
        favored_label,
        disadvantaged_label,
        n_group_filtered
      ))
    }
  } else if (!is.null(favored_group)) {
    if (length(grp_levels_raw) != 2) {
      stop(
        "Error: favored_group can only infer the comparison when group_var ",
        "contains exactly two groups. Use group_levels = c(favored, ",
        "disadvantaged) to select two levels explicitly.",
        call. = FALSE
      )
    }
    favored_chr <- as.character(favored_group)
    if (length(favored_chr) != 1 || !(favored_chr %in% grp_levels_raw)) {
      stop(
        "Error: favored_group must be one of the two levels of group_var: ",
        paste(grp_levels_raw, collapse = ", "), call. = FALSE
      )
    }
    disadvantaged_chr <- setdiff(grp_levels_raw, favored_chr)
    df[[group_var]] <- ifelse(as.character(df[[group_var]]) == favored_chr, 0, 1)
    favored_label <- favored_chr
    disadvantaged_label <- disadvantaged_chr
    group_selection_mode <- "favored_group"
    if (!quiet) {
      message(sprintf(
        "Note: favored_group = '%s' coded as 0 (reference); '%s' coded as 1. Gap = E(Y|favored) - E(Y|disadvantaged).",
        favored_label, disadvantaged_label
      ))
    }
  } else {
    if (length(grp_levels_raw) != 2) {
      stop(
        "Error: group_var must contain exactly two distinct groups unless ",
        "group_levels = c(favored, disadvantaged) selects the comparison.",
        call. = FALSE
      )
    }
    df[[group_var]] <- suppressWarnings(as.numeric(as.character(df[[group_var]])))
    if (any(is.na(df[[group_var]]))) {
      stop(
        "Error: group_var is not coded 0/1. Either recode it as 0 and 1, or ",
        "supply group_levels = c(favored, disadvantaged).",
        call. = FALSE
      )
    }
    if (!setequal(unique(df[[group_var]]), c(0, 1))) {
      stop(
        "Error: group_var must be coded as 0 and 1 (or use group_levels).",
        call. = FALSE
      )
    }
    favored_label <- "0"
    disadvantaged_label <- "1"
    group_selection_mode <- "implicit_0_1"
  }

  if (!is.null(weight_var)) {
    positive_weight <- df[[weight_var]] > 0
    n_zero_weight_dropped <- sum(!positive_weight)
    zero_weight_source_rows <- df$.f_source_row[!positive_weight]
    df <- droplevels(df[positive_weight, , drop = FALSE])
  } else {
    zero_weight_source_rows <- integer(0)
  }

  N_obs <- nrow(df)
  if (!quiet && n_missing_dropped > 0) {
    cat(sprintf(
      "\n[!] Data Prep: %d observations with missing values were omitted.\n",
      n_missing_dropped
    ))
  }
  if (!quiet && n_group_filtered > 0) {
    cat(sprintf(
      "[!] Group selection: %d observations outside group_levels were filtered. Active N = %d\n",
      n_group_filtered,
      N_obs
    ))
  }
  if (!quiet && n_zero_weight_dropped > 0) {
    cat(sprintf(
      paste0(
        "[!] Weight selection: %d observation(s) with zero weight were ",
        "excluded. Active N = %d\n"
      ),
      n_zero_weight_dropped,
      N_obs
    ))
  }
  if (N_obs == 0) {
    stop(
      "CRITICAL: Zero observations remaining after missing-data removal and group selection."
    )
  }

  unique_y <- unique(df[[dep_var]])
  unique_g <- unique(df[[group_var]])
  if (!all(unique_y %in% c(0, 1))) {
    stop("Error: dep_var must be strictly binary and coded as 0 and 1.", call. = FALSE)
  }
  if (!all(unique_g %in% c(0, 1))) {
    stop("Error: group_var must be strictly binary and coded as 0 and 1.", call. = FALSE)
  }
  if (!all(c(0, 1) %in% unique_y)) {
    stop("Error: dep_var must contain both outcome categories 0 and 1.", call. = FALSE)
  }
  if (!all(c(0, 1) %in% unique_g)) {
    if (!is.null(weight_var) && n_zero_weight_dropped > 0) {
      stop(
        "Error: both selected groups must contain at least one observation ",
        "with positive weight.",
        call. = FALSE
      )
    }
    stop("Error: group_var must contain observations in both groups 0 and 1.", call. = FALSE)
  }

  df$svy_weight <- if (use_svy && !is.null(weight_var)) as.numeric(df[[weight_var]]) else rep(1, N_obs)
  df$svy_strata <- if (use_svy && !is.null(strata_var)) as.character(df[[strata_var]]) else "1"
  df$svy_psu    <- if (use_svy && !is.null(psu_var)) as.character(df[[psu_var]]) else as.character(seq_len(N_obs))

  idx0 <- which(df[[group_var]] == 0); idx1 <- which(df[[group_var]] == 1)
  N_0 <- length(idx0); N_1 <- length(idx1)
  w0 <- df$svy_weight[idx0]; w1 <- df$svy_weight[idx1]
  reference_group0_weight <- switch(
    ref_method,
    reimers = 0.5,
    cotton = sum(w0) / sum(c(w0, w1)),
    NA_real_
  )

  # === 2. DESIGN MATRIX & DICCIONARIO INTELIGENTE ===
  form_base <- requested_formula

  X_full <- stats::model.matrix(form_base, data = df)
  x_cols <- colnames(X_full); nvars <- length(x_cols)
  X0 <- X_full[idx0, , drop = FALSE]; X1 <- X_full[idx1, , drop = FALSE]

  # The pooled reference uses a translated version of the internal 0/1 group
  # indicator. Keeping it in a dedicated column lets the coefficient extractor
  # omit the group-control coefficient while retaining the common intercept and
  # covariate coefficients used by the Fairlie counterfactual.
  pooled_group_column <- ".f_pooled_group"
  while (pooled_group_column %in% names(df)) {
    pooled_group_column <- paste0(pooled_group_column, "_")
  }

  assign_attr <- attr(X_full, "assign"); term_labels <- attr(stats::terms(form_base), "term.labels")
  display_dict <- list(); auto_groups <- list()

  for (i in seq_along(term_labels)) {
    term <- term_labels[i]; cols <- which(assign_attr == i)
    if(length(cols) > 0) {
      auto_groups[[term]] <- cols
      if (is.factor(df[[term]])) {
        ref_level <- levels(df[[term]])[1]
        for (cn in x_cols[cols]) display_dict[[cn]] <- sprintf("%s (%s vs. Ref: %s)", term, substring(cn, nchar(term) + 1), ref_level)
      } else {
        for (cn in x_cols[cols]) display_dict[[cn]] <- sprintf("%s (Continuous)", term)
      }
    }
  }
  display_dict[["(Intercept)"]] <- "Baseline Reference (Intercept)"
  diagnostics <- list()
  if (use_svy && vce_method == "bootstrap") {
    diagnostics$bootstrap_singleton_control <- paste0(
      "Methodology Note: lonely_psu does not control Rao-Wu bootstrap ",
      "replicate construction. Survey-bootstrap singleton strata are ",
      "controlled by bootstrap_singleton = '",
      bootstrap_singleton,
      "'."
    )
    if (!quiet) message(diagnostics$bootstrap_singleton_control)
  }

  if (length(groupings) == 0) {
    g_idx <- auto_groups
    group_names <- character(length(g_idx))
    for(k in seq_along(g_idx)) {
      term <- names(g_idx)[k]
      if(length(g_idx[[k]]) == 1) {
        group_names[k] <- display_dict[[ x_cols[g_idx[[k]]] ]]
      } else {
        ref_level <- levels(df[[term]])[1]
        group_names[k] <- sprintf("%s (Factor block: %d categories vs. Ref: %s)", term, length(g_idx[[k]]), ref_level)
      }
    }
  } else {
    g_idx <- list()
    # Match priority: exact model term > exact column name > prefix (literal).
    # Prefix matching can over-capture (e.g. "edu" also grabs "education..."),
    # so flag it when a prefix spans columns from more than one model term.
    ambiguous_prefixes <- character(0)
    unmatched_grouping_specs <- character(0)
    match_group_cols <- function(v) {
      if (v %in% names(auto_groups)) return(auto_groups[[v]])
      if (v %in% x_cols) return(which(x_cols == v))
      hits <- which(startsWith(x_cols, v))
      if (length(hits) > 0) {
        terms_hit <- unique(assign_attr[hits])
        if (length(terms_hit) > 1) ambiguous_prefixes <<- c(ambiguous_prefixes, v)
      }
      hits
    }
    for (g in names(groupings)) {
      matched_by_spec <- lapply(groupings[[g]], match_group_cols)
      unmatched <- lengths(matched_by_spec) == 0
      if (any(unmatched)) {
        unmatched_grouping_specs <- c(
          unmatched_grouping_specs,
          paste0(g, "=", groupings[[g]][unmatched])
        )
      }
      matched_cols <- unlist(matched_by_spec)
      matched_cols <- setdiff(matched_cols, which(x_cols == "(Intercept)"))
      if (length(matched_cols) > 0) g_idx[[g]] <- unique(matched_cols)
    }
    if (length(unmatched_grouping_specs) > 0) {
      stop(
        "Grouping Error: the following grouping specifications did not match ",
        "a model term, column, or prefix: ",
        paste(unmatched_grouping_specs, collapse = ", "),
        call. = FALSE
      )
    }
    if (length(ambiguous_prefixes) > 0) {
      diagnostics$groupings <- sprintf("Grouping Warning: prefix(es) %s matched columns from more than one model term. Verify the intended variables; use exact term names to disambiguate.", paste(sQuote(unique(ambiguous_prefixes)), collapse = ", "))
      if (!quiet) warning(diagnostics$groupings, call. = FALSE)
    }
    # Coverage guard: columns not covered by 'groupings' would otherwise be
    # silently absorbed by whichever block happens to be swapped last (an
    # order-dependent misattribution under randomize_order). Assign them to an
    # explicit residual block instead.
    covered <- unique(unlist(g_idx))
    uncovered <- setdiff(which(x_cols != "(Intercept)"), covered)
    if (length(uncovered) > 0) {
      g_idx[["Other (ungrouped predictors)"]] <- uncovered
      diagnostics$groupings_coverage <- sprintf(
        "Grouping Note: %d model column(s) not covered by 'groupings' were assigned to an 'Other (ungrouped predictors)' block: %s.",
        length(uncovered), paste(x_cols[uncovered], collapse = ", ")
      )
      if (!quiet) warning(diagnostics$groupings_coverage, call. = FALSE)
    }

    column_owners <- lapply(seq_along(x_cols), function(column) {
      names(g_idx)[vapply(
        g_idx,
        function(columns) column %in% columns,
        logical(1)
      )]
    })
    overlapping_columns <- which(lengths(column_owners) > 1)
    if (length(overlapping_columns) > 0) {
      overlap_description <- vapply(
        overlapping_columns,
        function(column) {
          sprintf(
            "%s {%s}",
            x_cols[[column]],
            paste(column_owners[[column]], collapse = ", ")
          )
        },
        character(1)
      )
      stop(
        "Grouping Error: model columns cannot belong to more than one block: ",
        paste(overlap_description, collapse = "; "),
        call. = FALSE
      )
    }

    split_terms <- character(0)
    for (term in names(auto_groups)) {
      term_columns <- auto_groups[[term]]
      if (length(term_columns) > 1) {
        term_owners <- unique(unlist(column_owners[term_columns]))
        if (length(term_owners) > 1) {
          split_terms <- c(split_terms, term)
        }
      }
    }
    if (length(split_terms) > 0) {
      stop(
        "Grouping Error: all columns generated by a multi-column factor term ",
        "must remain in one decomposition block. Split term(s): ",
        paste(split_terms, collapse = ", "),
        call. = FALSE
      )
    }
    group_names <- names(g_idx)
  }
  block_ids <- names(g_idx)
  k_groups <- length(g_idx)

  # === 3. MODEL VALIDATION SUITE ===
  # F8: SEs from a covariance diagonal without masking non-positive variances via
  # abs(); such terms (non-PSD covariance) are flagged and returned as NA.
  safe_se <- function(v) {
    neg <- which(v < -1e-8)
    if (length(neg) > 0) {
      diagnostics$neg_var <<- sprintf(
        "Numerical Warning: %d non-positive variance term(s) detected (non-PSD covariance); their SE set to NA.",
        length(neg)
      )
      if (!relax && !quiet) warning(diagnostics$neg_var)
    }
    v[!is.na(v) & v < 0] <- NA_real_
    sqrt(v)
  }

  ordinary_linearized_vce <- vce_method == "linearized" &&
    !complex_survey_design
  if (vce_method == "linearized") {
    diagnostics$vce <- if (ordinary_linearized_vce) {
      paste0(
        "Methodology Note: 'linearized' uses the Huber-White robust ",
        "sandwich with the maximum-likelihood N/(N-1) correction for ",
        "Logit/Probit coefficients."
      )
    } else {
      paste0(
        "Methodology Note: 'linearized' uses the declared survey design ",
        "for Taylor-linearized coefficient covariance."
      )
    }
    if (!quiet) message(diagnostics$vce)
  }

  if (ncol(X_full) >= min(N_0, N_1)) {
    diagnostics$df_obs <- sprintf("CRITICAL MODEL VALIDATION: Predictors (%d) exceeds observations in one of the groups (Min N = %d). Unidentified model.", ncol(X_full), min(N_0, N_1))
    if (!relax) stop(diagnostics$df_obs) else if (!quiet) warning(diagnostics$df_obs)
  }

  # F2: count PSUs within strata (nested), not globally.
  n_psu <- if (use_svy) {
    nrow(unique(df[, c("svy_strata", "svy_psu")]))
  } else {
    NA_integer_
  }
  psu_per_stratum <- if (use_svy) {
    table(unique(df[, c("svy_strata", "svy_psu")])$svy_strata)
  } else {
    integer(0)
  }
  bootstrap_singleton_strata <- if (
    use_svy && vce_method == "bootstrap"
  ) {
    names(psu_per_stratum)[psu_per_stratum < 2]
  } else {
    character(0)
  }
  bootstrap_has_singletons <- length(bootstrap_singleton_strata) > 0L
  if (
    bootstrap_has_singletons &&
      bootstrap_singleton == "fail"
  ) {
    stop(
      "Error: Rao-Wu survey bootstrap requires at least two PSUs in ",
      "every stratum. Set bootstrap_singleton = \"certainty\" only when ",
      "assigning zero first-stage variance to singleton strata is defensible; ",
      "otherwise collapse or redesign singleton strata, or use appropriate ",
      "externally supplied replicate weights.",
      call. = FALSE
    )
  }
  if (
    bootstrap_has_singletons &&
      bootstrap_singleton == "certainty"
  ) {
    diagnostics$bootstrap_singleton_assumption <- sprintf(
      paste0(
        "Bootstrap Assumption: %d singleton-PSU stratum/strata are held ",
        "fixed with replicate factor 1 and contribute zero first-stage ",
        "variance (bootstrap_singleton = 'certainty')."
      ),
      length(bootstrap_singleton_strata)
    )
    if (!quiet) {
      warning(
        diagnostics$bootstrap_singleton_assumption,
        call. = FALSE
      )
    }
  }
  if (use_svy && vce_method %in% c("jackknife", "bootstrap")) {
    diagnostics$survey_replication_scope <- paste0(
      "Survey Replication Scope: ultimate-cluster variance using one final ",
      "weight, optional strata, and one PSU stage; no FPC, lower-stage ",
      "identifiers, calibration replicate weights, or first-stage PPS ",
      "resampling is represented."
    )
    if (!quiet) message(diagnostics$survey_replication_scope)
  }
  if (complex_survey_design && ncol(X_full) >= n_psu) {
    diagnostics$df_psu <- sprintf("Model Validation Warning: Predictors (%d) exceeds or equals PSUs (%d). VCE matrix is not full rank.", ncol(X_full), n_psu)
    if (!relax && !quiet) warning(diagnostics$df_psu)
  }

  var_check0 <- apply(X0[, -1, drop=FALSE], 2, function(x) stats::var(x, na.rm=TRUE))
  var_check1 <- apply(X1[, -1, drop=FALSE], 2, function(x) stats::var(x, na.rm=TRUE))
  if (any(var_check0 == 0) || any(var_check1 == 0)) {
    bad_v <- unique(c(names(var_check0)[var_check0 == 0], names(var_check1)[var_check1 == 0]))
    diagnostics$zero_variance <- paste(
      "Model Note: zero within-group variance detected for:",
      paste(bad_v, collapse = ", "),
      ". The affected coefficient is allowed to be omitted and represented as 0, matching Stata."
    )
    if (!quiet) message(diagnostics$zero_variance)
  }

  exact_collinearity <- qr(
    X_full,
    tol = 1e-7,
    LAPACK = FALSE
  )$rank < ncol(X_full)
  if (exact_collinearity) {
    diagnostics$vif <- paste0(
      "Model Note: exact collinearity detected. Aliased coefficients are ",
      "allowed to be omitted and represented as 0, matching Stata."
    )
    if (!quiet) message(diagnostics$vif)
  }
  cor_X <- suppressWarnings(stats::cov.wt(X_full[, -1, drop = FALSE], wt = df$svy_weight, cor = TRUE)$cor)
  if (!exact_collinearity && !is.null(cor_X) && !any(is.na(cor_X))) {
    inv_cor <- tryCatch(solve(cor_X), error = function(e) NULL)
    if (is.null(inv_cor)) {
      diagnostics$vif <- paste0(
        "Model Note: exact collinearity detected. Aliased coefficients are ",
        "allowed to be omitted and represented as 0, matching Stata."
      )
      if (!quiet) message(diagnostics$vif)
    } else if (any(diag(inv_cor) > 15)) {
      diagnostics$vif <- "Model Validation Warning: Severe Multicollinearity detected (VIF > 15). Decomposed components may be unstable."
      if (!relax && !quiet) warning(diagnostics$vif)
    }
  }

  if (use_svy && !is.null(weight_var)) {
    med_w <- stats::median(df$svy_weight, na.rm = TRUE)
    if (med_w > 0 && (max(df$svy_weight, na.rm = TRUE) / med_w) > 50) {
      diagnostics$weights <- sprintf("Sample Design Warning: Highly skewed survey weights (Max/Median = %.1f). Resampling VCE strongly advised.", max(df$svy_weight, na.rm=TRUE)/med_w)
      if (!quiet) warning(diagnostics$weights)
    }
  }

  # === 4. BASE ESTIMATION ===
  fam <- if(model_type == "logit") stats::quasibinomial(link="logit") else stats::quasibinomial(link="probit")
  F_link <- fam$linkinv; f_dens <- fam$mu.eta
  ids_f <- if (is.null(psu_var)) ~1 else stats::as.formula(paste0("~", psu_var))
  strata_f <- if (is.null(strata_var)) NULL else stats::as.formula(paste0("~", strata_var))
  form_reg <- stats::as.formula(paste(dep_var, "~", paste(indep_vars, collapse = " + ")))
  fit_control <- stats::glm.control(
    epsilon = if (exact_collinearity) 1e-8 else 1e-16,
    maxit = 100
  )

  mcfadden_r2 <- function(m) {
      res <- tryCatch(1 - (m$deviance / m$null.deviance), error = function(e) NA)
      return(if (is.null(res)) NA else res)
  }

  estimate_system <- function(
      dataset, w_col = "svy_weight", compute_covariance = TRUE
  ) {
    i0 <- which(dataset[[group_var]] == 0 & dataset[[w_col]] > 0)
    i1 <- which(dataset[[group_var]] == 1 & dataset[[w_col]] > 0)
    need_g0 <- ref_method %in% c("group0", "reimers", "cotton")
    need_g1 <- ref_method %in% c("group1", "reimers", "cotton")
    need_pooled <- ref_method %in% c("pooled", "neumark")
    comparison_rows <- c(i0, i1)
    require_outcome_variation <- function(rows, label) {
      observed <- unique(dataset[[dep_var]][
        rows[dataset[[w_col]][rows] > 0]
      ])
      if (!all(c(0, 1) %in% observed)) {
        stop(
          "Error: the ", label, " model cannot be estimated because ",
          "dep_var does not contain both outcome categories 0 and 1 in ",
          "its estimation sample.",
          call. = FALSE
        )
      }
    }
    if (need_g0) require_outcome_variation(i0, "group 0")
    if (need_g1) require_outcome_variation(i1, "group 1")
    if (need_pooled) {
      require_outcome_variation(comparison_rows, "pooled reference")
    }
    pooled_group_mean <- sum(
      dataset[[w_col]][comparison_rows] *
        dataset[[group_var]][comparison_rows]
    ) / sum(dataset[[w_col]][comparison_rows])
    pooled_anchor_value <- switch(
      pooled_anchor,
      favored = 0,
      disadvantaged = 1,
      centered = pooled_group_mean
    )
    if (ref_method == "pooled") {
      dataset[[pooled_group_column]] <-
        dataset[[group_var]] - pooled_anchor_value
      form_pool <- stats::as.formula(paste(
        dep_var,
        "~",
        pooled_group_column,
        "+",
        paste(indep_vars, collapse = " + ")
      ))
    }

    # Robust extraction: glm/svyglm drop aliased (within-group collinear)
    # coefficients, so indexing coef()/vcov() directly by x_cols fails with an
    # obscure "subscript out of bounds". Pad to the full parameter set with 0s
    # (coefficient and variance) and flag it.
    aliased_flag <- FALSE
    robust_binary_vcov <- function(
        m,
        model_formula,
        model_data,
        model_weights
    ) {
      cf <- stats::coef(m)
      model_matrix <- stats::model.matrix(
        model_formula,
        data = model_data
      )
      keep <- names(cf)[!is.na(cf)]
      model_matrix <- model_matrix[, keep, drop = FALSE]
      response <- as.numeric(model_data[[dep_var]])
      beta <- cf[keep]
      eta <- as.vector(model_matrix %*% beta)
      mu <- if (model_type == "logit") {
        stats::plogis(eta)
      } else {
        stats::pnorm(eta)
      }
      mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
      variance <- mu * (1 - mu)
      derivative <- if (model_type == "logit") {
        variance
      } else {
        stats::dnorm(eta)
      }
      residual <- response - mu
      score_values <- model_weights * residual *
        derivative / variance
      hessian_weights <- model_weights *
        derivative^2 / variance
      if (model_type == "probit") {
        derivative_prime <- -eta * derivative
        hessian_weights <- model_weights * (
          derivative^2 / variance -
            residual * derivative_prime / variance +
            residual * derivative^2 * (1 - 2 * mu) /
              variance^2
        )
      }
      bread <- MASS::ginv(crossprod(
        model_matrix,
        model_matrix * hessian_weights
      ))
      covariance_kept <- bread %*%
        crossprod(model_matrix * score_values) %*%
        bread
      observations <- nrow(model_matrix)
      covariance_kept <- (
        observations / (observations - 1)
      ) * covariance_kept
      full_covariance <- matrix(
        0,
        nrow = length(cf),
        ncol = length(cf),
        dimnames = list(names(cf), names(cf))
      )
      full_covariance[keep, keep] <- covariance_kept
      full_covariance
    }

    survey_binary_vcov <- function(
        m,
        model_formula,
        model_data,
        model_weights,
        survey_design
    ) {
      cf <- stats::coef(m)
      model_matrix <- stats::model.matrix(
        model_formula,
        data = model_data
      )
      keep <- names(cf)[!is.na(cf)]
      model_matrix <- model_matrix[, keep, drop = FALSE]
      response <- as.numeric(model_data[[dep_var]])
      beta <- cf[keep]
      eta <- as.vector(model_matrix %*% beta)
      mu <- if (model_type == "logit") {
        stats::plogis(eta)
      } else {
        stats::pnorm(eta)
      }
      mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
      variance <- mu * (1 - mu)
      derivative <- if (model_type == "logit") {
        variance
      } else {
        stats::dnorm(eta)
      }
      residual <- response - mu
      score_values <- model_weights * residual *
        derivative / variance
      hessian_weights <- model_weights *
        derivative^2 / variance
      if (model_type == "probit") {
        derivative_prime <- -eta * derivative
        hessian_weights <- model_weights * (
          derivative^2 / variance -
            residual * derivative_prime / variance +
            residual * derivative^2 * (1 - 2 * mu) /
              variance^2
        )
      }
      bread <- MASS::ginv(crossprod(
        model_matrix,
        model_matrix * hessian_weights
      ))
      linearized_values <- (model_matrix * score_values) %*% bread
      sampling_weights <- as.numeric(stats::weights(
        survey_design,
        type = "sampling"
      ))
      total_variable_names <- paste0(
        ".f_score_if_",
        seq_len(ncol(linearized_values))
      )
      design_for_totals <- survey_design
      for (column in seq_along(total_variable_names)) {
        design_for_totals$variables[[total_variable_names[[column]]]] <-
          linearized_values[, column] / sampling_weights
      }
      total_formula <- stats::reformulate(total_variable_names)
      covariance_kept <- stats::vcov(survey::svytotal(
        total_formula,
        design_for_totals
      ))
      dimnames(covariance_kept) <- list(keep, keep)
      full_covariance <- matrix(
        0,
        nrow = length(cf),
        ncol = length(cf),
        dimnames = list(names(cf), names(cf))
      )
      full_covariance[keep, keep] <- covariance_kept
      full_covariance
    }

    extract_cv <- function(
        m,
        model_formula,
        model_data,
        model_weights,
        survey_design = NULL
    ) {
      cf <- stats::coef(m)
      V <- if (!compute_covariance) {
        matrix(
          0,
          nrow = length(cf),
          ncol = length(cf),
          dimnames = list(names(cf), names(cf))
        )
      } else if (ordinary_linearized_vce) {
        robust_binary_vcov(
          m,
          model_formula,
          model_data,
          model_weights
        )
      } else if (use_svy && !is.null(survey_design)) {
        survey_binary_vcov(
          m,
          model_formula,
          model_data,
          model_weights,
          survey_design
        )
      } else {
        stats::vcov(m)
      }
      b_out <- rep(0, nvars); names(b_out) <- x_cols
      V_out <- matrix(0, nvars, nvars, dimnames = list(x_cols, x_cols))
      keep <- intersect(names(cf)[!is.na(cf)], x_cols)
      b_out[keep] <- cf[keep]
      keepV <- intersect(rownames(V), keep)
      V_out[keepV, keepV] <- V[keepV, keepV, drop = FALSE]
      if (length(keep) < length(x_cols)) aliased_flag <<- TRUE
      list(b = b_out, V = V_out, full_V = V)
    }

    m0 <- NULL
    m1 <- NULL
    m_ref <- NULL
    b0 <- NULL
    b1 <- NULL
    V0 <- NULL
    V1 <- NULL
    V0_full <- NULL
    V1_full <- NULL
    V_ref_full <- NULL

    if (use_svy) {
      design_formula <- stats::as.formula(paste0("~", w_col))
      d_all <- survey::svydesign(
        ids = ids_f,
        strata = strata_f,
        weights = design_formula,
        data = dataset,
        nest = TRUE
      )
      if (need_g0) {
        # Build the design before selecting the domain. This retains the
        # original first-stage PSU counts used by Taylor linearization instead
        # of pretending that PSUs with no group-0 observations were never
        # sampled.
        d_g0 <- d_all[dataset[[group_var]] == 0 &
          dataset[[w_col]] > 0, ]
        m0 <- suppressWarnings(survey::svyglm(
          form_reg, design = d_g0, family = fam,
          control = stats::glm.control(epsilon = 1e-16, maxit = 100)
        ))
        e0 <- extract_cv(
          m0,
          form_reg,
          dataset[i0, , drop = FALSE],
          dataset[[w_col]][i0],
          d_g0
        )
        b0 <- e0$b
        V0 <- e0$V
        V0_full <- e0$full_V
      }
      if (need_g1) {
        d_g1 <- d_all[dataset[[group_var]] == 1 &
          dataset[[w_col]] > 0, ]
        m1 <- suppressWarnings(survey::svyglm(
          form_reg, design = d_g1, family = fam,
          control = stats::glm.control(epsilon = 1e-16, maxit = 100)
        ))
        e1 <- extract_cv(
          m1,
          form_reg,
          dataset[i1, , drop = FALSE],
          dataset[[w_col]][i1],
          d_g1
        )
        b1 <- e1$b
        V1 <- e1$V
        V1_full <- e1$full_V
      }
      if (ref_method == "neumark") {
        m_ref <- suppressWarnings(survey::svyglm(
          form_reg, design = d_all, family = fam,
          control = stats::glm.control(epsilon = 1e-16, maxit = 100)
        ))
        er <- extract_cv(
          m_ref,
          form_reg,
          dataset,
          dataset[[w_col]],
          d_all
        )
        b <- er$b; bV <- er$V; V_ref_full <- er$full_V
      } else if (ref_method == "pooled") {
        m_ref <- suppressWarnings(survey::svyglm(
          form_pool, design = d_all, family = fam,
          control = stats::glm.control(epsilon = 1e-16, maxit = 100)
        ))
        er <- extract_cv(
          m_ref,
          form_pool,
          dataset,
          dataset[[w_col]],
          d_all
        )
        b <- er$b; bV <- er$V; V_ref_full <- er$full_V
      }
    } else {
      d0 <- dataset[i0, , drop = FALSE]
      d1 <- dataset[i1, , drop = FALSE]
      d_all <- dataset
      d0$.lunadecomp_weight <- d0[[w_col]]
      d1$.lunadecomp_weight <- d1[[w_col]]
      d_all$.lunadecomp_weight <- d_all[[w_col]]

      if (need_g0) {
        m0 <- suppressWarnings(stats::glm(
          form_reg, data = d0, weights = .lunadecomp_weight, family = fam,
          control = fit_control
        ))
        e0 <- extract_cv(
          m0,
          form_reg,
          d0,
          d0$.lunadecomp_weight
        )
        b0 <- e0$b
        V0 <- e0$V
        V0_full <- e0$full_V
      }
      if (need_g1) {
        m1 <- suppressWarnings(stats::glm(
          form_reg, data = d1, weights = .lunadecomp_weight, family = fam,
          control = fit_control
        ))
        e1 <- extract_cv(
          m1,
          form_reg,
          d1,
          d1$.lunadecomp_weight
        )
        b1 <- e1$b
        V1 <- e1$V
        V1_full <- e1$full_V
      }
      if (ref_method == "neumark") {
        m_ref <- suppressWarnings(stats::glm(
          form_reg, data = d_all, weights = .lunadecomp_weight,
          family = fam, control = fit_control
        ))
        er <- extract_cv(
          m_ref,
          form_reg,
          d_all,
          d_all$.lunadecomp_weight
        )
        b <- er$b; bV <- er$V; V_ref_full <- er$full_V
      } else if (ref_method == "pooled") {
        m_ref <- suppressWarnings(stats::glm(
          form_pool, data = d_all, weights = .lunadecomp_weight,
          family = fam, control = fit_control
        ))
        er <- extract_cv(
          m_ref,
          form_pool,
          d_all,
          d_all$.lunadecomp_weight
        )
        b <- er$b; bV <- er$V; V_ref_full <- er$full_V
      }
    }

    if (ref_method == "group0") {
      b <- b0
      bV <- V0
      m_ref <- m0
      V_ref_full <- V0_full
    }
    if (ref_method == "group1") {
      b <- b1
      bV <- V1
      m_ref <- m1
      V_ref_full <- V1_full
    }
    if (ref_method %in% c("reimers", "cotton")) {
      # The Reimers/Cotton coefficient weight is a full-sample decomposition
      # setting. In particular, Cotton's observed weighted group share remains
      # fixed across sampling replicates instead of changing when a replicate
      # happens to contain more observations from one group.
      p0_ref <- reference_group0_weight
      b <- p0_ref * b0 + (1 - p0_ref) * b1
      bV <- (p0_ref^2) * V0 + ((1 - p0_ref)^2) * V1
    }

    required_models <- Filter(
      Negate(is.null),
      list(group0 = m0, group1 = m1, pooled = if (need_pooled) m_ref else NULL)
    )
    perfect_classification <- function(model) {
      coefficients <- stats::coef(model)
      probabilities <- stats::fitted(model)
      response <- model$y
      if (
        is.null(response) ||
          anyNA(probabilities) ||
          !any(abs(coefficients[!is.na(coefficients)]) > 15)
      ) {
        return(FALSE)
      }
      all((probabilities >= 0.5) == (response == 1))
    }
    separated_models <- names(required_models)[vapply(
      required_models,
      perfect_classification,
      logical(1)
    )]
    cvg <- all(vapply(
      required_models,
      function(model) isTRUE(model$converged),
      logical(1)
    ))
    out_aliased <- aliased_flag

    b[is.na(b)] <- 0
    bV[is.na(bV)] <- 0
    ref_full_b <- if (!is.null(m_ref)) stats::coef(m_ref) else b
    ref_full_V <- if (!is.null(V_ref_full)) V_ref_full else bV
    fitted_model_names <- names(required_models)

    out <- list(
      b = b,
      bV = bV,
      b0 = b0,
      b1 = b1,
      V0 = V0,
      V1 = V1,
      reference_model_b = ref_full_b,
      reference_model_V = ref_full_V,
      conv = cvg,
      separated_models = separated_models,
      aliased = out_aliased,
      fitted_models = fitted_model_names,
      pooled_anchor = list(
        requested = pooled_anchor,
        applicable = identical(ref_method, "pooled"),
        effective = if (identical(ref_method, "pooled")) {
          pooled_anchor
        } else {
          NA_character_
        },
        reference_value = if (identical(ref_method, "pooled")) {
          pooled_anchor_value
        } else {
          NA_real_
        },
        disadvantaged_share = pooled_group_mean,
        model_term = if (identical(ref_method, "pooled")) {
          pooled_group_column
        } else {
          NA_character_
        }
      )
    )
    out$r2 <- c(
      Group0 = if (!is.null(m0)) mcfadden_r2(m0) else NA_real_,
      Group1 = if (!is.null(m1)) mcfadden_r2(m1) else NA_real_,
      Reference = if (!is.null(m_ref)) mcfadden_r2(m_ref) else NA_real_
    )
    out$dev <- c(
      Group0 = if (!is.null(m0)) m0$deviance else NA_real_,
      Group1 = if (!is.null(m1)) m1$deviance else NA_real_,
      Reference = if (!is.null(m_ref)) m_ref$deviance else NA_real_
    )
    return(out)
  }

  sys_eval <- estimate_system(df)
  if (length(sys_eval$separated_models) > 0) {
    diagnostics$separation <- paste0(
      "CRITICAL MODEL VALIDATION: perfect separation detected in: ",
      paste(sys_eval$separated_models, collapse = ", "),
      ". Stata fairlie also rejects this estimation sample."
    )
    stop(diagnostics$separation, call. = FALSE)
  }
  if (!sys_eval$conv) {
    diagnostics$convergence <- "CRITICAL MODEL VALIDATION: GLM algorithm did not converge for at least one base group or reference."
    if (!relax) stop(diagnostics$convergence) else if (!quiet) warning(diagnostics$convergence)
  }
  b_star <- sys_eval$b; V_star <- sys_eval$bV
  if (isTRUE(sys_eval$aliased)) {
    diagnostics$aliased <- "Note: aliased (collinear within a group) coefficient(s) detected; they were set to 0 with zero variance. Check within-group collinearity."
    if (!quiet) message(diagnostics$aliased)
  }

  max_beta <- max(abs(c(sys_eval$b0[-1], sys_eval$b1[-1], b_star[-1])), na.rm = TRUE)
  if (!is.infinite(max_beta) && max_beta > 15) {
    diagnostics$separation <- sprintf("Separation Warning: Extreme coefficients detected (|beta| Max = %.1f). Risk of perfect separation (Hauck-Donner effect).", max_beta)
    if (!relax && !quiet) warning(diagnostics$separation)
  }

  # === 5. MACRO METRICS & COMMON SUPPORT AUDIT ===
  xb0_full <- as.vector(X0 %*% b_star); xb1_full <- as.vector(X1 %*% b_star)
  P0_full <- F_link(xb0_full); P1_full <- F_link(xb1_full)

  group0_prediction_range <- range(P0_full)
  group1_prediction_range <- range(P1_full)
  min_p_global <- max(
    group0_prediction_range[[1]],
    group1_prediction_range[[1]]
  )
  max_p_global <- min(
    group0_prediction_range[[2]],
    group1_prediction_range[[2]]
  )
  overlap_exists <- min_p_global <= max_p_global
  off_g0_flag <- if (overlap_exists) {
    P0_full < min_p_global | P0_full > max_p_global
  } else {
    rep(TRUE, length(P0_full))
  }
  off_g1_flag <- if (overlap_exists) {
    P1_full < min_p_global | P1_full > max_p_global
  } else {
    rep(TRUE, length(P1_full))
  }
  off_g0_unweighted <- mean(off_g0_flag)
  off_g1_unweighted <- mean(off_g1_flag)
  off_g0_weighted <- stats::weighted.mean(off_g0_flag, w0)
  off_g1_weighted <- stats::weighted.mean(off_g1_flag, w1)
  off_g0 <- if (!is.null(weight_var)) off_g0_weighted else off_g0_unweighted
  off_g1 <- if (!is.null(weight_var)) off_g1_weighted else off_g1_unweighted
  overlap_diagnostic <- list(
    quantity = "reference-model predicted outcome probability",
    causal_propensity_score = FALSE,
    trimming_applied = FALSE,
    group_0_range = stats::setNames(
      group0_prediction_range,
      c("min", "max")
    ),
    group_1_range = stats::setNames(
      group1_prediction_range,
      c("min", "max")
    ),
    intersection = stats::setNames(
      c(min_p_global, max_p_global),
      c("min", "max")
    ),
    overlap_exists = overlap_exists,
    off_overlap_unweighted = c(
      group_0 = off_g0_unweighted,
      group_1 = off_g1_unweighted
    ),
    off_overlap_weighted = c(
      group_0 = off_g0_weighted,
      group_1 = off_g1_weighted
    ),
    alert_threshold = 0.05
  )

  if (off_g0 > 0.05 || off_g1 > 0.05) {
    overlap_scale <- if (!is.null(weight_var)) "weighted" else "unweighted"
    diagnostics$overlap <- sprintf(
      paste0(
        "Predicted-Outcome Overlap Alert: limited overlap in reference-model ",
        "predicted outcome probabilities. Group 0 outside the shared range: ",
        "%.1f%%; Group 1: %.1f%% (%s proportions). This is descriptive only; ",
        "no observations were trimmed."
      ),
      off_g0 * 100,
      off_g1 * 100,
      overlap_scale
    )
    if (!relax && !quiet) warning(diagnostics$overlap)
  }

  p0_rank <- order(P0_full); p1_rank <- order(P1_full)
  mean_g0 <- stats::weighted.mean(df[[dep_var]][idx0], w0); mean_g1 <- stats::weighted.mean(df[[dep_var]][idx1], w1)
  diff_total <- mean_g0 - mean_g1
  tot_explained <- stats::weighted.mean(P0_full, w0) - stats::weighted.mean(P1_full, w1)
  tot_unexplained <- diff_total - tot_explained

  # F5: design-consistent SEs for the group means and the observed gap. The gap
  # SE comes from a design-based linear model of the outcome on the group
  # indicator, which correctly accounts for weighting, stratification and any
  # shared clustering across groups (so it is not assumed to be sqrt(se0^2+se1^2)).
  des_means <- survey::svydesign(ids = ids_f, strata = strata_f, weights = ~svy_weight, data = df, nest = TRUE)
  sm <- survey::svyby(stats::as.formula(paste0("~", dep_var)),
                      stats::as.formula(paste0("~", group_var)), des_means, survey::svymean)
  se_mean0 <- as.numeric(sm$se[sm[[group_var]] == 0])
  se_mean1 <- as.numeric(sm$se[sm[[group_var]] == 1])
  m_gap <- suppressWarnings(survey::svyglm(stats::as.formula(paste(dep_var, "~", group_var)), design = des_means))
  se_gap <- as.numeric(sqrt(stats::vcov(m_gap)[2, 2]))  # coef on group = mean1 - mean0

  # === 6. ITERATIVE PERMUTATION ENGINE (FAIRLIE) ===
  n_samp <- if (use_svy) {
    trunc((nrow(X0) + nrow(X1)) / 2)
  } else {
    min(nrow(X0), nrow(X1))
  }

  sample_ecdf <- function(n, weights) {
    if (length(weights) <= 1) return(rep(1, n))
    ub <- cumsum(weights) / sum(weights)
    u <- sort(stats::runif(n)); res <- integer(n); j <- 1
    for (i in seq_len(n)) { while (u[i] > ub[j] && j < length(ub)) j <- j + 1; res[i] <- j }
    return(res)
  }

  sample_rank_positions <- function(population_size, sample_size) {
    if (sample_size >= population_size) {
      return(seq_len(population_size))
    }
    sort(sample.int(population_size, sample_size, replace = FALSE))
  }

core_simulation <- function(runs, beta, vcov_beta, w_0, w_1, ord_0, ord_1, m_0, m_1, n_s, calc_var = TRUE) {
    c_res <- numeric(k_groups); v_res <- numeric(k_groups)
    c_m2 <- numeric(k_groups)
    dens_left <- NULL; dens_right <- NULL
    x_left <- NULL; x_right <- NULL
    prob_left <- NULL; prob_right <- NULL

    for (i in seq_len(runs)) {
      if (use_svy) {
        s0 <- ord_0[sample_ecdf(n_s, w_0[ord_0])]; s1 <- ord_1[sample_ecdf(n_s, w_1[ord_1])]
      } else {
        s0 <- ord_0[sample_rank_positions(length(ord_0), n_s)]
        s1 <- ord_1[sample_rank_positions(length(ord_1), n_s)]
      }

      k_order <- if (randomize_order) sample(seq_len(k_groups)) else seq_len(k_groups)
      v_sel <- rep(0, nvars)

      for (j in seq_len(k_groups)) {
        jj <- k_order[j]
        if (j == 1) {
          x_left <- m_0[s0, , drop = FALSE]; eta_left <- as.vector(x_left %*% beta); prob_left <- F_link(eta_left)
          if (calc_var) dens_left <- f_dens(eta_left)
        } else {
          x_left <- x_right; prob_left <- prob_right; if (calc_var) dens_left <- dens_right
        }

        if (j == k_groups) {
          x_right <- m_1[s1, , drop = FALSE]
        } else {
          v_sel[g_idx[[jj]]] <- 1; idx_act <- which(v_sel == 1); idx_pas <- which(v_sel == 0)
          x_right <- matrix(0, nrow = n_s, ncol = nvars)
          if (length(idx_act) > 0) x_right[, idx_act] <- m_1[s1, idx_act, drop = FALSE]
          if (length(idx_pas) > 0) x_right[, idx_pas] <- m_0[s0, idx_pas, drop = FALSE]
        }

        eta_right <- as.vector(x_right %*% beta); prob_right <- F_link(eta_right)
        if (calc_var) dens_right <- f_dens(eta_right)

        delta_p <- mean(prob_left - prob_right)
        contribution_delta <- delta_p - c_res[jj]
        c_res[jj] <- c_res[jj] + contribution_delta / i
        c_m2[jj] <- c_m2[jj] +
          contribution_delta * (delta_p - c_res[jj])

        if (calc_var) {
          jacob_mat <- colMeans(as.vector(dens_left) * x_left - as.vector(dens_right) * x_right)
          var_update <- as.numeric(t(jacob_mat) %*% vcov_beta %*% jacob_mat)
          v_res[jj] <- v_res[jj] + (var_update - v_res[jj]) / i
        }
      }
    }
    monte_carlo_se <- if (runs > 1) {
      sqrt((c_m2 / (runs - 1)) / runs)
    } else {
      rep(NA_real_, k_groups)
    }
    return(list(c = c_res, v = v_res, mc_se = monte_carlo_se))
  }

  if (!quiet) {
    cat(sprintf("\n[*] Initializing Nonlinear Decomposition Framework (%d iterations)...\n", reps))
    if (!randomize_order) cat("    -> Variable sequence locked (Path Dependence Mitigation Disabled).\n")
  }

  # Design-based SEs for the group means and observed gap are shared by both
  # engines (they describe observed statistics, not the decomposition components).
  se_mean0_ov <- se_mean0; se_mean1_ov <- se_mean1; se_gap_ov <- se_gap
  overall_vcov_raw <- NULL
  overall_vcov_complete <- FALSE
  detailed_vcov_complete <- FALSE
  replication_details <- NULL

  if (vce_method == "linearized") {
    if (!quiet) { cat("    -> VCE Engine: Analytic Delta Method\n"); pb <- utils::txtProgressBar(min = 0, max = 1, style = 3) }
    res_lin <- core_simulation(reps, b_star, V_star, w0, w1, p0_rank, p1_rank, X0, X1, n_samp, TRUE)
    if (!quiet) { utils::setTxtProgressBar(pb, 1); close(pb) }
    mean_contribs <- res_lin$c; se_contribs <- safe_se(res_lin$v)
    matching_mc_se <- res_lin$mc_se
    detailed_vcov_raw <- diag(
      as.numeric(res_lin$v),
      nrow = k_groups,
      ncol = k_groups
    )
    dimnames(detailed_vcov_raw) <- list(block_ids, block_ids)
    # Exact delta-method SE for the full-sample Total Explained. Summing per-block
    # variances is wrong on two counts: it ignores cross-block covariances (all
    # blocks share V_star) and it targets the matched-subsample quantity rather
    # than the full-sample Explained reported in the overall table.
    dens0_full <- f_dens(xb0_full); dens1_full <- f_dens(xb1_full)
    J_tot <- as.numeric(colSums(X0 * (dens0_full * w0)) / sum(w0) - colSums(X1 * (dens1_full * w1)) / sum(w1))
    se_tot_exp <- safe_se(as.numeric(t(J_tot) %*% V_star %*% J_tot))
    # F4/F5: in linearized mode the unexplained SE combines the (design) gap SE and
    # the (delta) explained SE assuming independence; this is a conservative
    # approximation. Replicate engines return an exact, covariance-aware SE.
    se_unexp_ov <- sqrt(se_gap_ov^2 + se_tot_exp^2)
    diagnostics$linearized_unexp <- "Methodology Note: In 'linearized' mode the Total Unexplained SE is a conservative delta approximation (independence of gap and explained assumed). Use 'jackknife'/'bootstrap' for exact covariance-aware inference."
  } else {
    ordinary_resampling <- !use_svy
    survey_generator_scale <- NA_real_
    survey_replicate_type <- NA_character_
    survey_replicate_engine <- NA_character_
    replicate_center_groups <- NULL
    if (ordinary_resampling && vce_method == "bootstrap") {
      # Generate the ordinary observation-level samples with boot::boot().
      # Fairlie matching still runs in the existing replicate pipeline below:
      # the integer multiplicities must be expanded into physical rows before
      # refitting and matching. simple = TRUE preserves the established
      # seed-to-replicate mapping without retaining an N x R index array.
      bootstrap_samples <- boot::boot(
        data = seq_len(N_obs),
        statistic = function(data, indices) {
          tabulate(indices, nbins = N_obs)
        },
        R = boot_reps,
        sim = "ordinary",
        stype = "i",
        simple = TRUE
      )
      replicate_weights <- t(bootstrap_samples$t)
      replicate_scale <- 1 / (boot_reps - 1)
      replicate_rscales <- rep(1, boot_reps)
    } else if (ordinary_resampling) {
      # resample::jackknife() supplies the deterministic leave-one-out
      # samples. Isolate its internal RNG initialization so the subsequent
      # Fairlie matching seed remains identical to the established pipeline.
      preserve_rng_state <- function(expression) {
        had_rng_state <- exists(".Random.seed", envir = .GlobalEnv)
        rng_state <- if (had_rng_state) {
          get(".Random.seed", envir = .GlobalEnv)
        } else {
          NULL
        }
        on.exit({
          if (had_rng_state) {
            assign(".Random.seed", rng_state, envir = .GlobalEnv)
          } else if (exists(".Random.seed", envir = .GlobalEnv)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        }, add = TRUE)
        force(expression)
      }
      row_ids <- seq_len(N_obs)
      jackknife_samples <- preserve_rng_state(
        resample::jackknife(
          data = row_ids,
          statistic = function(retained_rows) {
            if (length(retained_rows) == N_obs) {
              return(rep(1, N_obs))
            }
            weights <- rep(N_obs / (N_obs - 1), N_obs)
            weights[
              setdiff(row_ids, retained_rows)
            ] <- 0
            weights
          },
          statisticNames = paste0("weight_", row_ids),
          trace = FALSE
        )
      )
      replicate_weights <- unname(
        t(jackknife_samples$replicates)
      )
      replicate_scale <- (N_obs - 1) / N_obs
      replicate_rscales <- rep(1, N_obs)
    } else {
      des_base <- survey::svydesign(
        ids = ids_f,
        strata = strata_f,
        weights = ~svy_weight,
        data = df,
        nest = TRUE
      )
      des_rep_type <- if (
        use_svy && !is.null(strata_var)
      ) {
        "JKn"
      } else {
        "JK1"
      }
      des_rep <- if (vce_method == "bootstrap") {
        if (
          bootstrap_has_singletons &&
            bootstrap_singleton == "certainty"
        ) {
          # RWYB extends Rao--Wu to permit a final-stage singleton. The
          # singleton factor is fixed at 1, while regular strata retain n_h - 1
          # resampling. SRSWR states the same negligible-FPC/ultimate-cluster
          # approximation used by the established subbootstrap route.
          survey_replicate_engine <- "svrep::as_bootstrap_design"
          svrep::as_bootstrap_design(
            des_base,
            type = "Rao-Wu-Yue-Beaumont",
            replicates = boot_reps,
            mse = FALSE,
            samp_method_by_stage = "SRSWR"
          )
        } else {
          # Preserve the previously validated Stata-compatible generator
          # exactly whenever every stratum contains at least two PSUs.
          survey_replicate_engine <- "survey::as.svrepdesign"
          survey::as.svrepdesign(
            des_base,
            type = "subbootstrap",
            replicates = boot_reps
          )
        }
      } else {
        survey_replicate_engine <- "survey::as.svrepdesign"
        survey::as.svrepdesign(des_base, type = des_rep_type)
      }
      replicate_factors <- as.matrix(des_rep$repweights)
      replicate_weights <- replicate_factors * df$svy_weight
      survey_generator_scale <- des_rep$scale
      survey_replicate_type <- if (vce_method == "bootstrap") {
        if (
          bootstrap_has_singletons &&
            bootstrap_singleton == "certainty"
        ) {
          paste0(
            "Rao-Wu-Yue-Beaumont bootstrap with singleton certainty"
          )
        } else {
          "Rao-Wu rescaled bootstrap"
        }
      } else {
        des_rep_type
      }
      replicate_scale <- if (vce_method == "bootstrap") {
        1 / ncol(replicate_factors)
      } else {
        survey_generator_scale
      }
      replicate_rscales <- des_rep$rscales
      if (vce_method == "jackknife" && !is.null(strata_var)) {
        candidate_groups <- vapply(
          seq_len(ncol(replicate_factors)),
          function(rep_index) {
            changed <- abs(
              replicate_factors[, rep_index] - 1
            ) > 1e-12
            changed_strata <- unique(df$svy_strata[changed])
            if (length(changed_strata) == 1L) {
              changed_strata
            } else {
              NA_character_
            }
          },
          character(1)
        )
        if (!anyNA(candidate_groups)) {
          replicate_center_groups <- candidate_groups
        }
      }
    }
    n_replicates <- ncol(replicate_weights)
    matching_seed <- if (!is.null(seed)) {
      as.integer(seed)
    } else {
      sample.int(.Machine$integer.max, 1)
    }
    with_matching_seed <- function(expression) {
      had_seed <- exists(".Random.seed", envir = .GlobalEnv)
      previous_seed <- if (had_seed) {
        get(".Random.seed", envir = .GlobalEnv)
      } else {
        NULL
      }
      on.exit({
        if (had_seed) {
          assign(".Random.seed", previous_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(matching_seed)
      force(expression)
    }

    if (!quiet) {
      cat(sprintf(
        "    -> VCE Engine: Empirical %s (%d replications)\n",
        toupper(vce_method),
        n_replicates
      ))
      pb <- utils::txtProgressBar(
        min = 0,
        max = n_replicates,
        style = 3
      )
      env_pb <- new.env()
      env_pb$count <- 0
    }
    base_matching <- with_matching_seed(core_simulation(
      reps,
      b_star,
      V_star,
      w0,
      w1,
      p0_rank,
      p1_rank,
      X0,
      X1,
      n_samp,
      FALSE
    ))
    mean_contribs <- base_matching$c
    matching_mc_se <- base_matching$mc_se

    # Replicates return per-block contributions plus the group means and explained
    # component, so all overall SEs (means, gap, explained, unexplained) are
    # obtained consistently from one covariance.
    n_ret <- k_groups + 3
    rep_wrapper <- function(wt_rep, dataset) {
      if (!quiet) {
        env_pb$count <- env_pb$count + 1
        utils::setTxtProgressBar(
          pb,
          min(env_pb$count, n_replicates)
        )
      }
      tryCatch({
        if (ordinary_resampling && vce_method == "bootstrap") {
          multiplicity <- as.integer(round(wt_rep))
          expanded_rows <- rep(
            seq_len(nrow(dataset)),
            times = multiplicity
          )
          if (length(expanded_rows) == 0) {
            return(rep(NA_real_, n_ret))
          }
          replicate_data <- dataset[expanded_rows, , drop = FALSE]
          replicate_data$svy_weight <- 1
        } else {
          positive_rows <- which(wt_rep > 0)
          if (length(positive_rows) == 0) {
            return(rep(NA_real_, n_ret))
          }
          replicate_data <- dataset[positive_rows, , drop = FALSE]
          replicate_data$svy_weight <- if (ordinary_resampling) {
            1
          } else {
            wt_rep[positive_rows]
          }
        }
        sys_rep <- estimate_system(
          replicate_data,
          "svy_weight",
          compute_covariance = FALSE
        )
        if (!sys_rep$conv) return(rep(NA_real_, n_ret))

        i0_r <- which(replicate_data[[group_var]] == 0)
        i1_r <- which(replicate_data[[group_var]] == 1)
        if (length(i0_r) == 0 || length(i1_r) == 0) {
          return(rep(NA_real_, n_ret))
        }
        w0_r <- replicate_data$svy_weight[i0_r]
        w1_r <- replicate_data$svy_weight[i1_r]
        X_rep <- stats::model.matrix(
          form_base,
          data = replicate_data
        )
        if (!identical(colnames(X_rep), x_cols)) {
          X_aligned <- matrix(
            0,
            nrow = nrow(X_rep),
            ncol = length(x_cols),
            dimnames = list(NULL, x_cols)
          )
          common_columns <- intersect(colnames(X_rep), x_cols)
          X_aligned[, common_columns] <-
            X_rep[, common_columns, drop = FALSE]
          X_rep <- X_aligned
        }
        X0_r <- X_rep[i0_r, , drop = FALSE]
        X1_r <- X_rep[i1_r, , drop = FALSE]

        pr_0 <- F_link(as.vector(X0_r %*% sys_rep$b))
        pr_1 <- F_link(as.vector(X1_r %*% sys_rep$b))
        n_s_r <- if (!ordinary_resampling && use_svy) {
          trunc((nrow(X0_r) + nrow(X1_r)) / 2)
        } else {
          min(nrow(X0_r), nrow(X1_r))
        }

        sim_r <- with_matching_seed(core_simulation(
          reps,
          sys_rep$b,
          sys_rep$bV,
          w0_r,
          w1_r,
          order(pr_0),
          order(pr_1),
          X0_r,
          X1_r,
          n_s_r,
          FALSE
        ))
        mean0_r <- stats::weighted.mean(
          replicate_data[[dep_var]][i0_r],
          w0_r
        )
        mean1_r <- stats::weighted.mean(
          replicate_data[[dep_var]][i1_r],
          w1_r
        )
        exp_r <- stats::weighted.mean(pr_0, w0_r) - stats::weighted.mean(pr_1, w1_r)
        return(c(sim_r$c, mean0_r, mean1_r, exp_r))
      }, error = function(e) return(rep(NA_real_, n_ret)))
    }

    replicate_estimates <- vapply(
      seq_len(n_replicates),
      function(rep_index) {
        rep_wrapper(
          replicate_weights[, rep_index],
          df
        )
      },
      numeric(n_ret)
    )
    res_mat <- list(
      theta = c(mean_contribs, mean_g0, mean_g1, tot_explained),
      replicates = t(replicate_estimates)
    )
    if (!quiet) close(pb)

    valid_idx <- stats::complete.cases(res_mat$replicates)
    fail_rate <- sum(!valid_idx) / length(valid_idx)
    if (!ordinary_resampling && any(!valid_idx) && !relax) {
      stop(
        "CRITICAL: ",
        sum(!valid_idx),
        " survey replicate(s) failed. Dropping a PSU replicate can ",
        "invalidate the design variance. Inspect the sparse/separated ",
        "replicates or rerun with relax = TRUE for an explicitly ",
        "exploratory approximation.",
        call. = FALSE
      )
    }
    if (fail_rate > 0.05) {
      diagnostics$resampling <- sprintf("Statistical Alert: %.1f%% of resampling iterations failed to converge. SEs might be less reliable.", fail_rate * 100)
      if (!relax && fail_rate > 0.2) stop(diagnostics$resampling) else if (!quiet) warning(diagnostics$resampling)
    }
    if (sum(valid_idx) < 2) stop("CRITICAL: Resampling variance estimation failed due to sparsity.")

    valid_replicates <- res_mat$replicates[
      valid_idx,
      ,
      drop = FALSE
    ]
    colnames(valid_replicates) <- c(
      block_ids,
      "Group_0",
      "Group_1",
      "Explained"
    )
    valid_center_groups <- if (!is.null(replicate_center_groups)) {
      replicate_center_groups[valid_idx]
    } else {
      NULL
    }
    if (is.null(valid_center_groups)) {
      center_vec <- colMeans(valid_replicates)
      diffs <- sweep(
        valid_replicates,
        2,
        center_vec,
        "-"
      )
      replicate_center_method <- "replicate mean"
    } else {
      group_centers <- lapply(
        unique(valid_center_groups),
        function(group) {
          colMeans(
            valid_replicates[
              valid_center_groups == group,
              ,
              drop = FALSE
            ]
          )
        }
      )
      names(group_centers) <- unique(valid_center_groups)
      center_vec <- do.call(
        rbind,
        lapply(valid_center_groups, function(group) {
          group_centers[[group]]
        })
      )
      colnames(center_vec) <- colnames(valid_replicates)
      diffs <- valid_replicates - center_vec
      replicate_center_method <- "stratum-specific replicate mean"
    }
    rscales <- if (
      length(replicate_rscales) == length(valid_idx)
    ) {
      replicate_rscales[valid_idx]
    } else {
      replicate_rscales
    }
    # With relax = TRUE, this adjustment is only an exploratory safeguard; the
    # default survey path stops above rather than silently dropping a PSU
    # replicate.
    rep_adjust <- length(valid_idx) / sum(valid_idx)
    if (rep_adjust > 1) {
      diagnostics$rep_adjust <- sprintf("Variance Note: %d failed replicate(s) dropped; variance rescaled by %.3f to compensate.", sum(!valid_idx), rep_adjust)
    }
    effective_scale <- replicate_scale * rep_adjust
    V_CC <- effective_scale * crossprod(diffs * sqrt(rscales))
    se_contribs <- safe_se(diag(V_CC)[1:k_groups])
    detailed_vcov_raw <- V_CC[
      seq_len(k_groups),
      seq_len(k_groups),
      drop = FALSE
    ]
    dimnames(detailed_vcov_raw) <- list(block_ids, block_ids)
    detailed_vcov_complete <- TRUE
    im0 <- k_groups + 1; im1 <- k_groups + 2; iex <- k_groups + 3
    se_mean0_ov <- safe_se(V_CC[im0, im0]); se_mean1_ov <- safe_se(V_CC[im1, im1])
    se_tot_exp <- safe_se(V_CC[iex, iex])
    L_gap <- c(1, -1, 0); L_unx <- c(1, -1, -1)          # gap = m0-m1 ; unexp = m0-m1-exp
    Vsub <- V_CC[c(im0, im1, iex), c(im0, im1, iex)]
    se_gap_ov <- safe_se(as.numeric(t(L_gap) %*% Vsub %*% L_gap))
    se_unexp_ov <- safe_se(as.numeric(t(L_unx) %*% Vsub %*% L_unx))
    overall_transform <- rbind(
      Group_0 = c(1, 0, 0),
      Group_1 = c(0, 1, 0),
      Difference = c(1, -1, 0),
      Explained = c(0, 0, 1),
      Unexplained = c(1, -1, -1)
    )
    overall_vcov_raw <- overall_transform %*% Vsub %*%
      t(overall_transform)
    overall_vcov_complete <- TRUE
    detailed_overall_vcov <- V_CC[
      seq_len(k_groups),
      c(im0, im1, iex),
      drop = FALSE
    ] %*% t(overall_transform)
    replicate_output_terms <- c(
      rownames(overall_transform),
      paste0("Exp_", block_ids)
    )
    replicate_full_vcov <- rbind(
      cbind(overall_vcov_raw, t(detailed_overall_vcov)),
      cbind(detailed_overall_vcov, detailed_vcov_raw)
    )
    dimnames(replicate_full_vcov) <- list(
      replicate_output_terms,
      replicate_output_terms
    )
    replicate_overall_estimates <- cbind(
      Group_0 = valid_replicates[, "Group_0"],
      Group_1 = valid_replicates[, "Group_1"],
      Difference = valid_replicates[, "Group_0"] -
        valid_replicates[, "Group_1"],
      Explained = valid_replicates[, "Explained"],
      Unexplained = valid_replicates[, "Group_0"] -
        valid_replicates[, "Group_1"] -
        valid_replicates[, "Explained"]
    )
    replicate_output_estimates <- cbind(
      replicate_overall_estimates,
      valid_replicates[, block_ids, drop = FALSE]
    )
    colnames(replicate_output_estimates)[
      seq.int(ncol(replicate_overall_estimates) + 1, ncol(replicate_output_estimates))
    ] <- paste0("Exp_", block_ids)
    if (is.matrix(center_vec)) {
      replicate_output_center <- cbind(
        Group_0 = center_vec[, "Group_0"],
        Group_1 = center_vec[, "Group_1"],
        Difference = center_vec[, "Group_0"] -
          center_vec[, "Group_1"],
        Explained = center_vec[, "Explained"],
        Unexplained = center_vec[, "Group_0"] -
          center_vec[, "Group_1"] -
          center_vec[, "Explained"],
        center_vec[, block_ids, drop = FALSE]
      )
      colnames(replicate_output_center)[
        seq.int(6, ncol(replicate_output_center))
      ] <- paste0("Exp_", block_ids)
    } else {
      replicate_output_center <- colMeans(replicate_output_estimates)
    }
    replication_details <- list(
      center_method = replicate_center_method,
      center = replicate_output_center,
      center_groups = valid_center_groups,
      replicates = replicate_output_estimates,
      vcov = replicate_full_vcov,
      replicate_weights = replicate_weights[, valid_idx, drop = FALSE],
      requested_replicates = length(valid_idx),
      valid_replicates = sum(valid_idx),
      failed_replicates = sum(!valid_idx),
      scale = effective_scale,
      generator_scale = if (is.na(survey_generator_scale)) {
        replicate_scale
      } else {
        survey_generator_scale
      },
      rscales = rscales,
      ordinary_resampling = ordinary_resampling,
      ordinary_bootstrap = ordinary_resampling &&
        vce_method == "bootstrap",
      ordinary_jackknife = ordinary_resampling &&
        vce_method == "jackknife",
      engine = if (
        ordinary_resampling && vce_method == "bootstrap"
      ) {
        "boot::boot"
      } else if (ordinary_resampling) {
        "resample::jackknife"
      } else {
        survey_replicate_engine
      },
      simulation = if (
        ordinary_resampling && vce_method == "bootstrap"
      ) {
        "ordinary"
      } else if (ordinary_resampling) {
        "delete-one"
      } else {
        NULL
      },
      statistic_type = if (
        ordinary_resampling && vce_method == "bootstrap"
      ) {
        "indices"
      } else if (ordinary_resampling) {
        "retained row indices"
      } else {
        NULL
      },
      fit_weight_type = if (
        ordinary_resampling && vce_method == "bootstrap"
      ) {
        "frequency"
      } else if (ordinary_resampling) {
        "JK1 replicate factor"
      } else {
        "survey replicate weight"
      },
      survey_replicate_type = survey_replicate_type,
      bootstrap_singleton = list(
        requested = bootstrap_singleton,
        applied = !ordinary_resampling &&
          vce_method == "bootstrap" &&
          bootstrap_has_singletons &&
          bootstrap_singleton == "certainty",
        n_singleton_strata = length(bootstrap_singleton_strata),
        singleton_strata = bootstrap_singleton_strata,
        first_stage_variance = if (
          !ordinary_resampling &&
            vce_method == "bootstrap" &&
            bootstrap_has_singletons &&
            bootstrap_singleton == "certainty"
        ) {
          "zero for singleton strata"
        } else {
          "not modified"
        }
      ),
      matching_reps_per_replicate = reps,
      common_random_numbers = TRUE,
      matching_seed = matching_seed
    )
  }

  raw_mean_contribs <- stats::setNames(
    as.numeric(mean_contribs),
    block_ids
  )
  raw_se_contribs <- stats::setNames(
    as.numeric(se_contribs),
    block_ids
  )
  raw_matching_mc_se <- stats::setNames(
    as.numeric(matching_mc_se),
    block_ids
  )
  sum_blocks <- sum(raw_mean_contribs)
  matching_difference <- sum_blocks - tot_explained
  matching_tolerance <- 1e-12 * max(1, abs(tot_explained))
  if (abs(matching_difference) > matching_tolerance) {
    diagnostics$matching_target <- sprintf(
      paste0(
        "Matching Note: detailed contributions sum to %.12g, whereas the ",
        "full-sample Total Explained is %.12g (difference %.12g). ",
        "Contributions are reported without rescaling, matching Ben Jann's ",
        "estimator."
      ),
      sum_blocks,
      tot_explained,
      matching_difference
    )
  }
  names(mean_contribs) <- group_names
  names(se_contribs) <- group_names

  # === 7. PROFESSIONAL EXPORT & OUTPUT GENERATION ===
  q_val <- 1 - ((1 - level) / 2)
  # Replicate survey inference uses PSU-minus-strata degrees of freedom,
  # including a weights-only design where observations are the first-stage
  # units. Ordinary robust linearization and ordinary bootstrap are
  # normal-asymptotic; ordinary delete-one jackknife retains N-1.
  design_df <- if (use_svy && vce_method != "linearized") {
    n_strata_design <- length(unique(df$svy_strata))
    max(1, n_psu - n_strata_design)
  } else if (vce_method == "jackknife") {
    max(1, nrow(df) - 1)
  } else {
    Inf
  }

  # --- results_overall: 5 rows with full inference (harmonized with oby_decomp) ---
  overall_terms <- c(
    "Group_0", "Group_1", "Difference", "Explained", "Unexplained"
  )
  ov_est <- c(mean_g0, mean_g1, diff_total, tot_explained, tot_unexplained)
  names(ov_est) <- overall_terms
  ov_se  <- as.numeric(c(se_mean0_ov, se_mean1_ov, se_gap_ov, se_tot_exp, se_unexp_ov))
  names(ov_se) <- overall_terms
  ov_stat <- ov_est / ifelse(is.na(ov_se) | ov_se == 0, NA_real_, ov_se)
  ov_p   <- 2 * stats::pt(
    abs(ov_stat),
    df = design_df,
    lower.tail = FALSE
  )
  ov_lo  <- ov_est - stats::qt(q_val, df = design_df) * ov_se
  ov_hi  <- ov_est + stats::qt(q_val, df = design_df) * ov_se
  if (is.null(overall_vcov_raw)) {
    overall_vcov_raw <- matrix(
      NA_real_,
      nrow = length(overall_terms),
      ncol = length(overall_terms),
      dimnames = list(overall_terms, overall_terms)
    )
    diag(overall_vcov_raw) <- ov_se^2
  }
  tb_overall <- tibble::tibble(
    Term = overall_terms,
    Estimate = unname(ov_est),
    Std_Error = unname(ov_se),
    Statistic = unname(ov_stat),
    P_Value = unname(ov_p),
    Conf_Low = unname(ov_lo),
    Conf_High = unname(ov_hi)
  )

  # --- results_detailed_explained: per variable/factor block (Fairlie) ---
  pct_contribution_explained <- if (abs(tot_explained) > 1e-12) {
    (mean_contribs / tot_explained) * 100
  } else {
    rep(NA_real_, length(mean_contribs))
  }
  pct_contribution_total <- if (abs(diff_total) > 1e-12) {
    (mean_contribs / diff_total) * 100
  } else {
    rep(NA_real_, length(mean_contribs))
  }

  stat_vec <- mean_contribs / ifelse(se_contribs == 0 | is.na(se_contribs), NA_real_, se_contribs)
  tb_detailed <- tibble::tibble(
    Term = names(mean_contribs),
    Estimate = unname(mean_contribs),
    `% Contribution Explained` = unname(pct_contribution_explained),
    `% Contribution Total` = unname(pct_contribution_total),
    Std_Error = unname(se_contribs),
    Statistic = unname(stat_vec),
    P_Value = unname(2 * stats::pt(
      abs(stat_vec),
      df = design_df,
      lower.tail = FALSE
    )),
    Conf_Low = unname(
      mean_contribs - stats::qt(q_val, df = design_df) * se_contribs
    ),
    Conf_High = unname(
      mean_contribs + stats::qt(q_val, df = design_df) * se_contribs
    )
  ) %>% dplyr::arrange(dplyr::desc(abs(.data$Estimate)))

  build_coef_table <- function(beta, V) {
    if (is.null(beta) || is.null(V)) return(NULL)
    se <- safe_se(diag(V))
    stat <- beta / ifelse(is.na(se) | se == 0, NA_real_, se)
    p_v <- 2 * stats::pt(
      abs(stat),
      df = design_df,
      lower.tail = FALSE
    )
    terms_disp <- unname(sapply(names(beta), function(x) ifelse(!is.null(display_dict[[x]]), display_dict[[x]], x)))
    tibble::tibble(
      Term = terms_disp,
      Estimate = unname(beta),
      Std_Error = unname(se),
      Statistic = unname(stat),
      P_Value = unname(p_v)
    )
  }

  tb_coefs <- Filter(
    Negate(is.null),
    list(
      Group_0 = build_coef_table(sys_eval$b0, sys_eval$V0),
      Group_1 = build_coef_table(sys_eval$b1, sys_eval$V1),
      Reference = build_coef_table(b_star, V_star)
    )
  )
  fit_metrics <- tibble::tibble(
    Model = c("Group 0", "Group 1", "Reference"),
    Estimated = c(
      !is.null(sys_eval$b0),
      !is.null(sys_eval$b1),
      TRUE
    ),
    Pseudo_R2 = unname(sys_eval$r2),
    Deviance = unname(sys_eval$dev)
  )

  if (!quiet) {
    cat("\n", rep("=", 95), "\n NONLINEAR DECOMPOSITION ANALYSIS (FAIRLIE FRAMEWORK)\n", rep("=", 95), "\n", sep="")
    cat(sprintf(" Number of Strata = %-15s N (Group 0) = %d\n", ifelse(use_svy, length(unique(df$svy_strata)), "N/A"), N_0))
    cat(sprintf(" Number of PSUs   = %-15s N (Group 1) = %d\n", ifelse(use_svy, n_psu, "N/A"), N_1))
    cat(sprintf(" Ref. Method      = %-15s Pr(Y=1|G=0) = %.6f\n", toupper(ref_method), mean_g0))
    cat(sprintf(" VCE Engine       = %-15s Pr(Y=1|G=1) = %.6f\n", toupper(vce_method), mean_g1))
    cat(sprintf(" Permutations     = %-15d Difference  = %.6f\n", reps, diff_total))
    cat(sprintf(" Random Ordering  = %-15s Design df   = %s\n", as.character(randomize_order), design_df))
    cat(sprintf(" Group 0 (favored)      : %s = %s\n", group_var, favored_label))
    cat(sprintf(" Group 1 (disadvantaged): %s = %s\n", group_var, disadvantaged_label))
    cat(" Gap orientation : Difference = E(Y|favored) - E(Y|disadvantaged)\n")
    if (identical(ref_method, "pooled")) {
      cat(sprintf(
        " Pooled anchor   : %s (group-indicator reference value = %.6f)\n",
        pooled_anchor,
        sys_eval$pooled_anchor$reference_value
      ))
    }
    cat(rep("-", 95), "\n", sep="")

    if (length(diagnostics) > 0) {
        cat("\n [MODEL VALIDATION & DIAGNOSTICS SUITE]\n")
        for (diag_name in names(diagnostics)) cat(" *", diagnostics[[diag_name]], "\n\n")
        cat(rep("-", 95), "\n", sep="")
    } else {
        cat("\n [DIAGNOSTICS: ALL ECONOMETRIC ASSUMPTIONS PASSED]\n"); cat(rep("-", 95), "\n", sep="")
    }

    cat("\n [DECOMPOSITION SUMMARY]\n"); print(as.data.frame(tb_overall), row.names = FALSE)
    cat("\n [DETAILED EXPLAINED COMPONENTS]\n"); print(as.data.frame(tb_detailed[, c("Term", "Estimate", "% Contribution Explained", "% Contribution Total", "Std_Error", "Statistic", "P_Value")]), row.names = FALSE)
    cat(rep("-", 95), "\n\n", sep="")
  }

  raw_detail_names <- paste0("Exp_", block_ids)
  names(raw_mean_contribs) <- raw_detail_names
  names(raw_se_contribs) <- raw_detail_names
  names(raw_matching_mc_se) <- raw_detail_names
  dimnames(detailed_vcov_raw) <- list(
    raw_detail_names,
    raw_detail_names
  )
  raw_estimates <- c(ov_est, raw_mean_contribs)
  raw_standard_errors <- c(ov_se, raw_se_contribs)
  raw_vcov <- matrix(
    NA_real_,
    nrow = length(raw_estimates),
    ncol = length(raw_estimates),
    dimnames = list(names(raw_estimates), names(raw_estimates))
  )
  raw_vcov[overall_terms, overall_terms] <- overall_vcov_raw
  raw_vcov[raw_detail_names, raw_detail_names] <- detailed_vcov_raw
  if (!is.null(replication_details)) {
    raw_vcov <- replication_details$vcov[
      names(raw_estimates),
      names(raw_estimates),
      drop = FALSE
    ]
  }

  return(invisible(list(
    summary_stats = list(
      N = N_obs, Pop = sum(df$svy_weight),
      Strata = if (use_svy) length(unique(df$svy_strata)) else NA,
      PSUs = n_psu, DF = design_df, N0 = N_0, N1 = N_1,
      favored = favored_label, disadvantaged = disadvantaged_label,
      group_selection_mode = group_selection_mode,
      available_group_levels = available_group_levels,
      selected_group_levels = c(favored_label, disadvantaged_label),
      excluded_group_levels = excluded_group_levels,
      n_missing_dropped = n_missing_dropped,
      n_group_filtered = n_group_filtered,
      n_zero_weight_dropped = n_zero_weight_dropped
    ),
    model_metrics = fit_metrics,
    models_coefficients = tb_coefs,
    diagnostics = diagnostics,
    results_overall = tb_overall,
    results_detailed_explained = tb_detailed,
    results_detailed_unexplained = NULL,
    raw = list(
      estimates = raw_estimates,
      standard_errors = raw_standard_errors,
      vcov = raw_vcov,
      overall_estimates = ov_est,
      overall_standard_errors = ov_se,
      overall_vcov = overall_vcov_raw,
      overall_vcov_complete = overall_vcov_complete,
      detailed_estimates = raw_mean_contribs,
      detailed_standard_errors = raw_se_contribs,
      detailed_vcov = detailed_vcov_raw,
      detailed_vcov_complete = detailed_vcov_complete,
      group_coefficients = list(
        group_0 = sys_eval$b0,
        group_1 = sys_eval$b1
      ),
      group_vcov = list(
        group_0 = sys_eval$V0,
        group_1 = sys_eval$V1
      ),
      reference_coefficients = b_star,
      reference_vcov = V_star,
      reference_model_coefficients = sys_eval$reference_model_b,
      reference_model_vcov = sys_eval$reference_model_V,
      reference_linear_predictor = list(
        group_0 = stats::setNames(
          xb0_full,
          df$.f_source_row[idx0]
        ),
        group_1 = stats::setNames(
          xb1_full,
          df$.f_source_row[idx1]
        )
      ),
      reference_predictions = list(
        group_0 = stats::setNames(
          P0_full,
          df$.f_source_row[idx0]
        ),
        group_1 = stats::setNames(
          P1_full,
          df$.f_source_row[idx1]
        )
      ),
      analytic_sample = df$.f_source_row,
      sample_flow = list(
        initial_rows = seq_len(N_initial),
        missing_rows = as.integer(missing_source_rows),
        group_filtered_rows = as.integer(group_filtered_source_rows),
        zero_weight_rows = as.integer(zero_weight_source_rows),
        analytic_rows = as.integer(df$.f_source_row)
      ),
      model_terms = x_cols,
      decomposition_blocks = stats::setNames(
        lapply(g_idx, function(columns) x_cols[columns]),
        block_ids
      ),
      fitted_models = sys_eval$fitted_models,
      matching = list(
        sample_size = n_samp,
        reps = reps,
        randomize_order = randomize_order,
        sampling_method = if (use_svy) {
          "ranked_pps_with_replacement"
        } else if (N_0 == N_1) {
          "complete_rank_match"
        } else {
          "ranked_srs_without_replacement"
        },
        sampled_group = if (use_svy) {
          "both"
        } else if (N_0 > N_1) {
          "group_0"
        } else if (N_1 > N_0) {
          "group_1"
        } else {
          "none"
        },
        group_sizes = c(group_0 = N_0, group_1 = N_1),
        stochastic_subsampling = if (
          use_svy
        ) {
          TRUE
        } else {
          N_0 != N_1
        },
        detailed_monte_carlo_se = raw_matching_mc_se,
        detailed_sum = sum_blocks,
        full_sample_explained = tot_explained,
        difference = matching_difference,
        rescaled = FALSE
      ),
      replication = replication_details,
      group_mapping = list(
        favored = favored_label,
        disadvantaged = disadvantaged_label,
        selected = c(favored_label, disadvantaged_label),
        available = available_group_levels,
        excluded = excluded_group_levels,
        selection_mode = group_selection_mode,
        n_filtered = n_group_filtered
      ),
      pooled_anchor = sys_eval$pooled_anchor,
      reference_structure = list(
        method = ref_method,
        fairlie_native_comparator = ref_method %in%
          c("group0", "group1", "neumark", "pooled"),
        validation_status = if (
          ref_method %in% c("reimers", "cotton")
        ) {
          "package extension validated by internal identities"
        } else {
          "validated against Ben Jann's fairlie 1.0.7"
        },
        group_0_coefficient_weight = reference_group0_weight,
        group_1_coefficient_weight = if (
          ref_method %in% c("reimers", "cotton")
        ) {
          1 - reference_group0_weight
        } else {
          NA_real_
        },
        weight_fixed_across_sampling_replicates =
          ref_method %in% c("reimers", "cotton"),
        reference_weight_inference = if (
          ref_method %in% c("reimers", "cotton")
        ) {
          "conditional on the full-sample reference weight"
        } else {
          "not applicable"
        }
      ),
      interaction_handling = list(
        formula_interactions = "rejected",
        precomputed_interactions = paste0(
          "not detectable; valid only when all constituent main effects and ",
          "interaction columns are exchanged in one joint groupings block"
        ),
        separate_hierarchical_attribution = FALSE
      ),
      predicted_outcome_overlap = overlap_diagnostic,
      weighting = list(
        active = !is.null(weight_var),
        variable = if (!is.null(weight_var)) {
          weight_var
        } else {
          NA_character_
        },
        interpretation = if (!is.null(weight_var)) {
          "pweight"
        } else {
          "none"
        },
        zero_weight_rows_excluded = n_zero_weight_dropped,
        total_weight = sum(df$svy_weight),
        group_total_weight = c(
          group_0 = sum(w0),
          group_1 = sum(w1)
        ),
        scale_invariant_point_estimator = TRUE
      ),
      model_type = model_type,
      model_family = "quasibinomial",
      ref_method = ref_method,
      vce_method = vce_method,
      linearized_vce = list(
        coefficient_covariance = if (ordinary_linearized_vce) {
          "Huber-White maximum-likelihood sandwich"
        } else if (vce_method == "linearized") {
          "Taylor survey linearization"
        } else {
          "not applicable"
        },
        finite_sample_correction = if (ordinary_linearized_vce) {
          "N/(N-1) by fitted Logit/Probit equation"
        } else if (vce_method == "linearized") {
          "survey design convention"
        } else {
          "not applicable"
        },
        probit_bread = if (
          ordinary_linearized_vce && model_type == "probit"
        ) {
          "observed information"
        } else if (
          ordinary_linearized_vce && model_type == "logit"
        ) {
          "observed information (equal to Fisher information for Logit)"
        } else {
          "not applicable"
        },
        inference_distribution = if (is.infinite(design_df)) {
          "normal"
        } else {
          paste0("t(", design_df, ")")
        }
      ),
      survey_mode = use_svy,
      complex_survey_design = complex_survey_design,
      survey_design_scope = list(
        variance_structure = if (use_svy) {
          "ultimate cluster"
        } else {
          "not applicable"
        },
        final_weight = !is.null(weight_var),
        strata = !is.null(strata_var),
        psu = !is.null(psu_var),
        psu_per_stratum = psu_per_stratum,
        finite_population_correction = FALSE,
        lower_sampling_stages = FALSE,
        first_stage_pps_resampling = FALSE,
        external_replicate_weights = FALSE,
        zero_weight_records_retained_for_variance = FALSE,
        lonely_psu = lonely_psu,
        bootstrap_singleton = bootstrap_singleton,
        bootstrap_singleton_strata = bootstrap_singleton_strata
      ),
      survey_linearization = if (complex_survey_design) {
        "not available for the complete ranked Fairlie estimator; use jackknife or bootstrap"
      } else {
        "not applicable"
      }
    )
  )))
}
