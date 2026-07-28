#' Oaxaca-Blinder and Yun Decomposition for Two-Group Outcome Gaps
#'
#' @description
#' Decomposes mean outcome differences between two groups using
#' Oaxaca-Blinder-type counterfactual decompositions for ordinary least-squares
#' and linear probability models, and Yun-type decompositions for nonlinear
#' binary-response models.
#'
#' The function estimates group-specific models and partitions the observed
#' outcome gap into explained and unexplained components. For continuous
#' outcomes, use `model_type = "ols"`. For binary outcomes estimated with a
#' linear probability model, use `model_type = "lpm"`. For nonlinear
#' binary-response decompositions on the response-probability scale, use
#' `model_type = "logit"` or `model_type = "probit"`.
#'
#' @details
#' Let `group_var` define two comparison groups coded as 0 and 1. The function
#' estimates separate outcome models for each group and constructs a
#' counterfactual reference structure according to `ref_method`.
#'
#' For `model_type = "ols"`, the decomposition is performed on the outcome mean
#' scale using a Gaussian identity-link model. For `model_type = "lpm"`, the
#' function applies the same linear decomposition to a binary outcome, yielding
#' components on the probability scale. For `model_type = "logit"` or
#' `model_type = "probit"`, the decomposition is performed on the response
#' probability scale using the inverse-link transformation.
#'
#' The total observed gap is reported as:
#'
#' \deqn{
#'   \Delta = E(Y \mid G = 0) - E(Y \mid G = 1)
#' }
#'
#' and is partitioned into an explained component and an unexplained component.
#' The explained component represents the part of the gap associated with
#' differences in observed covariates. The unexplained component represents the
#' residual part associated with differences in coefficients and unobserved
#' factors. In applied inequality research, this component should not be
#' interpreted mechanically as discrimination or inequity unless the substantive
#' assumptions are justified.
#'
#' @section Group coding and the direction of the gap:
#' The decomposition always compares two selected levels of `group_var`.
#' Internally the favored (socially advantaged) level is coded 0 and the
#' disadvantaged level is coded 1, and the reported gap is always
#'
#' \deqn{\Delta = E(Y \mid \text{favored}) - E(Y \mid \text{disadvantaged}).}
#'
#' Which group is "favored" is a substantive choice, not a property of the data:
#' it follows directly from the formula above. Whatever level is assigned to code
#' 0 becomes the minuend and the reference structure whose endowments are
#' revalued in the explained component. The sign of \eqn{\Delta} therefore
#' depends jointly on the coding and on the polarity of the outcome: for a
#' desirable outcome (e.g. insurance coverage) the favored group typically yields
#' \eqn{\Delta > 0}, whereas for an adverse outcome (e.g. mortality, unmet need)
#' the favored group may yield \eqn{\Delta < 0}. "Favored" refers to social
#' advantage, not to the numeric level of the outcome.
#'
#' There are three ways to declare the comparison:
#' \itemize{
#'   \item Recommended: pass an ordered two-element vector through
#'   `group_levels = c(favored, disadvantaged)`. The two selected levels are
#'   retained, all other levels of `group_var` are filtered out, and the
#'   selected levels are recoded internally to 0 and 1.
#'   \item For a variable that already contains exactly two groups, pass
#'   `favored_group` with the advantaged level. The other observed level is
#'   inferred as disadvantaged. This interface is retained for backward
#'   compatibility.
#'   \item Code `group_var` yourself as 0/1. If both `group_levels` and
#'   `favored_group` are `NULL`, group 0 is treated as the favored/reference
#'   group by convention.
#' }
#' The selected mapping is reported in the console, in
#' `summary_stats$favored` / `summary_stats$disadvantaged`, and in
#' `raw$group_mapping`.
#'
#' In nonlinear models (`"logit"`, `"probit"`) the split between explained and
#' unexplained depends on the point of the inverse-link curve at which the
#' counterfactual is evaluated. When `ref_method = "pooled"`, `pooled_anchor`
#' makes that choice explicit. The default `"favored"` uses the favored group's
#' baseline index, `"disadvantaged"` uses the disadvantaged group's baseline
#' index, and `"centered"` uses the weighted mean of the two group positions.
#' The latter is symmetric with respect to the selected groups but represents a
#' synthetic position on the link scale. In linear models (`"ols"`, `"lpm"`)
#' the reported decomposition is invariant to this choice because the reference
#' intercept cancels.
#'
#' `ref_method = "neumark"` omits the group indicator from the reference model.
#' This may transfer residual group differences into the pooled slopes when
#' group membership is associated with the predictors, thereby increasing the
#' explained component (Jann, 2008). The choice of anchor should therefore
#' follow the counterfactual question rather than the sign or desirability of
#' the outcome.
#'
#' @section Reference coefficient structure:
#' The argument `ref_method` defines the counterfactual coefficient vector used
#' for the decomposition:
#'
#' \itemize{
#'   \item `"pooled"`: coefficients from a pooled model including the group
#'   indicator. The reference intercept is selected by `pooled_anchor`.
#'   \item `"neumark"`: coefficients from a pooled model excluding the group
#'   indicator.
#'   \item `"cotton"`: weighted average of group-specific coefficients using
#'   group sample weights.
#'   \item `"reimers"`: equally weighted average of group-specific coefficients.
#'   \item `"group1"`: coefficients from the group coded as 0.
#'   \item `"group2"`: coefficients from the group coded as 1.
#' }
#'
#' @section Categorical predictors:
#' Factor and character predictors are expanded into treatment-contrast
#' indicators using `model.matrix()`. A factor with \eqn{K} observed levels
#' contributes \eqn{K-1} indicators, and the first factor level is the omitted
#' reference category. Each indicator receives its own detailed contribution
#' and standard error.
#'
#' Supplying the original factor term in `groupings`, for example
#' `groupings = list(Education = "education")`, aggregates all of its generated
#' indicators. The grouped standard error is calculated from the full
#' covariance submatrix, including covariance between indicator contributions.
#'
#' Without normalization, detailed categorical contributions, especially the
#' unexplained allocation between indicators and the intercept, can depend on
#' the omitted reference category. Use `normalize` to express the effects of
#' selected factors as deviation contrasts from their unweighted grand mean
#' (Yun, 2005). Supply factor names individually, for example
#' `normalize = c("education", "insurance")`, or use `normalize = "all"` to
#' normalize every factor included in `indep_vars`.
#'
#' A normalized factor contributes one detailed row for every observed level,
#' including the originally omitted level. The transformation preserves fitted
#' values and the overall decomposition, makes normalized detailed results
#' invariant to the choice of omitted category, and propagates the complete
#' coefficient-and-mean covariance matrix. The grand mean is an equal-weight
#' mean over categories, not a prevalence-weighted mean.
#'
#' If a normalized factor interacts with a continuous numeric variable, the
#' interaction is normalized automatically and independently. Its
#' equal-category mean is transferred to the continuous main-effect
#' coefficient, just as the mean of the factor's main effects is transferred to
#' the intercept. For example, with
#' `indep_vars = c("education * experience")` and
#' `normalize = "education"`, both `education` and
#' `education:experience` are normalized. The model must be hierarchical: the
#' factor main effect, continuous main effect, and interaction must all be
#' present. Normalization of factor-by-factor and higher-order interactions is
#' not currently supported.
#'
#' @section Survey designs:
#' If `use_svy = TRUE`, the function uses the supplied survey weights, strata,
#' and primary sampling units to account for complex survey designs. If
#' `weight_var` is supplied and `use_svy = FALSE`, survey mode is activated
#' automatically. `weight_var` is interpreted as a sampling weight (the
#' counterpart of Stata's `pweight`), not as a frequency, analytic, or
#' importance weight. Multiplying every sampling weight by the same positive
#' constant leaves estimates and their design-based covariance unchanged.
#' Missing weights are excluded with the other incomplete analysis variables;
#' zero weights are retained in the analytic observation count but contribute
#' no population mass. Negative and infinite weights are rejected.
#'
#' PSU identifiers are nested within strata, so the same PSU labels may be
#' reused in different strata. For a stratum containing one PSU,
#' `lonely_psu = "adjust"` matches Stata's `singleunit(centered)`;
#' `"certainty"` and `"remove"` assign a zero between-PSU contribution and
#' match `singleunit(certainty)`; `"fail"` stops; and `"average"` imputes the
#' average variance contribution of the non-lonely strata, matching Stata's
#' `singleunit(scaled)`.
#'
#' @section Variance estimation:
#' The variance engine is selected with `vce_method`:
#'
#' \itemize{
#'   \item `"linearized"`: delta-method variance outside survey mode and
#'   Taylor-linearized design variance in survey mode. Outside survey mode it
#'   always uses the Huber--White sandwich covariance, matching
#'   `oaxaca, vce(robust)` in comparable settings.
#'   \item `"jackknife"`: replicate-weight jackknife variance estimation.
#'   \item `"bootstrap"`: bootstrap replicate variance estimation.
#' }
#'
#' In survey mode the Taylor-linearized covariance generalizes the
#' model-robust sandwich covariance by estimating variation from the declared
#' weights, strata, and PSUs. Consequently, there is no separate robust switch:
#' `linearized` selects the appropriate robust covariance for the sampling
#' context.
#'
#' Outside survey mode, `"linearized"` and `"bootstrap"` use standard-normal
#' inference; ordinary delete-one `"jackknife"` uses a Student t distribution
#' with \eqn{N-1} degrees of freedom. In survey mode all three methods use a
#' Student t distribution with design degrees of freedom equal to the number
#' of nested PSUs minus the number of strata. `level` controls the confidence
#' level. Two-sided p-values are evaluated from the upper tail directly so that
#' very small, nonzero probabilities are not rounded spuriously to zero.
#'
#' Outside survey mode the Huber--White covariance is robust to
#' heteroskedasticity under independent observations. It does not correct an
#' incorrectly specified outcome model, dependence between observations,
#' outliers, or an incorrectly declared survey design.
#'
#' To reproduce the finite-sample conventions of `oaxaca, vce(robust)`, the
#' group-specific OLS/LPM sandwiches use the `regress` correction
#' \eqn{N_g/(N_g-k)}, whereas group-specific Logit/Probit sandwiches use the
#' maximum-likelihood correction \eqn{N_g/(N_g-1)}. Pooled-reference joint
#' covariances use the corresponding `suest` correction. The Probit bread is
#' constructed from the observed-information Hessian used by Stata, rather than
#' the Fisher-information matrix used internally by standard GLM Fisher
#' scoring.
#'
#' For final inferential analyses with complex survey data, replicate-based
#' variance estimation is generally preferable when computationally feasible.
#'
#' Two standard approximations apply to the linearized method. First, for
#' nonlinear models (`"logit"`, `"probit"`) the variance of the group-level
#' predictions propagates the covariate means through the \emph{average}
#' derivative of the inverse link rather than through the individual-level
#' influence of each prediction, so the standard errors of `Group_0`, `Group_1`,
#' and `Difference` are first-order approximations; replicate methods capture
#' this variability exactly. Second, the covariance between estimated
#' coefficients and covariate means is set to zero, the standard assumption
#' under correct model specification (Jann, 2008).
#'
#' Ordinary bootstrap resamples observations with replacement. Survey
#' bootstrap uses the Rao--Wu rescaled construction: within each stratum it
#' samples \eqn{n_h-1} PSUs with replacement and multiplies their sampling
#' weights by \eqn{n_h/(n_h-1)}. Its final variance uses Stata's
#' `svy bootstrap` convention for one bootstrap sample per replicate-weight
#' variable: deviations from the replicate mean with scale \eqn{1/B}. The
#' native `survey` generator scale \eqn{1/(B-1)} is retained in
#' `raw$replication$generator_scale` for auditability. Rao--Wu requires at
#' least two PSUs in every stratum and assumes simple or stratified random
#' sampling of PSUs under the ultimate-cluster/with-replacement approximation;
#' it is not a first-stage PPS bootstrap.
#'
#' Ordinary and unstratified jackknife variances use the replicate-mean
#' convention. Stratified JKn uses stratum-specific replicate means, which is
#' algebraically equivalent to Stata's default pseudovalue formula. Bootstrap
#' also uses the replicate mean. These are the default, non-MSE conventions in
#' Stata. Cotton's full-sample group share is held fixed across replications,
#' matching `oaxaca, weight(#)`. If any replicates fail, the variance is
#' rescaled by the ratio of total to valid replicates and a diagnostic is
#' recorded.
#'
#' Ordinary delete-one jackknife samples are generated by
#' [resample::jackknife()]. Each retained observation receives the JK1 factor
#' \eqn{N/(N-1)}, which leaves the replicate point estimate unchanged but
#' preserves the exact Stata/`survey` JK1 representation. The covariance uses
#' deviations from the replicate mean and scale \eqn{(N-1)/N}.
#'
#' @param data A data frame containing the outcome, group, predictors, and
#' optional survey design variables.
#' @param dep_var Character string. Name of the dependent variable.
#' @param group_var Character string. Name of the comparison variable. It may
#' contain two or more observed levels when `group_levels` explicitly selects
#' the two levels to compare. Otherwise it must contain exactly two groups.
#' @param favored_group Optional. The level of `group_var` that identifies the
#' favored (socially advantaged) group when `group_var` already contains exactly
#' two observed groups. The other level is inferred as disadvantaged. Cannot
#' be combined with `group_levels`. Retained for backward compatibility; for
#' new analyses, prefer the explicit ordered `group_levels` interface.
#' @param group_levels Optional ordered vector of length two:
#' `c(favored, disadvantaged)`. The first value defines group 0 and the second
#' defines group 1, so the reported gap is
#' `Difference = E(Y | favored) - E(Y | disadvantaged)`. Both values must occur
#' in `group_var` and must be distinct. If `group_var` contains additional
#' levels, their observations are filtered out before estimation.
#' @param indep_vars Character vector. Predictor names or model terms to include
#' in the decomposition model. Interactions may be supplied using ordinary R
#' formula syntax, for example `"education * experience"` or the explicitly
#' expanded terms `c("education", "experience",
#' "education:experience")`. Categorical predictors should be coded as factors
#' or character variables.
#' @param groupings Optional named list. Each element should contain one or more
#' predictor names or prefixes to aggregate detailed contributions into broader
#' conceptual domains.
#' @param normalize Optional character vector naming factor predictors whose
#' categorical effects should be normalized. Use `"all"` by itself to normalize
#' every factor in `indep_vars`; use `NULL` (default) for treatment coding with
#' the first level omitted. `"all"` cannot be combined with individual names.
#' Factor-by-continuous interactions involving a selected factor are normalized
#' automatically when the corresponding main effects are present.
#' @param ref_method Character string. Counterfactual reference structure.
#' Options are `"pooled"`, `"neumark"`, `"cotton"`, `"reimers"`, `"group1"`,
#' and `"group2"`.
#' @param pooled_anchor Character string. Baseline position of the group
#' indicator used to construct the pooled reference coefficients. Options are
#' `"favored"` (default), `"disadvantaged"`, and `"centered"`. With
#' `"centered"`, the internal indicator is centered at its analysis-weighted
#' mean (or survey-weighted population mean when weights are supplied).
#' This argument is used only when `ref_method = "pooled"`. It can materially
#' change Logit/Probit decompositions; OLS/LPM decompositions are invariant.
#' @param model_type Character string. Model family. Options are `"ols"`,
#' `"lpm"`, `"logit"`, and `"probit"`. Use `"ols"` for continuous outcomes,
#' `"lpm"` for binary outcomes estimated as a linear probability model, and
#' `"logit"` or `"probit"` for nonlinear binary-response models.
#' @param use_svy Logical. If `TRUE`, applies complex survey design settings.
#' @param weight_var Optional character string. Name of the sampling weight
#' variable. Non-missing weights must be numeric, finite, and non-negative.
#' Zero weights are retained in the analytic sample, matching Stata survey
#' conventions, but each selected group must have positive total weight.
#' @param strata_var Optional character string. Name of the stratification
#' variable.
#' @param psu_var Optional character string. Name of the primary sampling unit
#' variable.
#' @param vce_method Character string. Variance estimation method. Options are
#' `"linearized"`, `"jackknife"`, and `"bootstrap"`.
#' @param boot_reps Integer. Number of bootstrap replications when
#' `vce_method = "bootstrap"`. Without a survey design, ordinary
#' observation-level resampling is performed by [boot::boot()]. Bootstrap
#' observation indices are converted internally to frequency weights before
#' fitting each replicate. Survey bootstrap continues to use replicate weights
#' generated from the declared survey design.
#' @param lonely_psu Character string. Handling of strata with one PSU. Options
#' are `"fail"`, `"remove"`, `"certainty"`, `"adjust"`, and `"average"`.
#' @param level Numeric. Confidence level for confidence intervals. Default is
#' `0.95`.
#' @param seed Optional integer. Random seed for reproducible replicate variance
#' estimation. Defaults to `NULL` (non-deterministic); set an explicit value for
#' reproducibility. The prior state of the global random number generator is
#' restored on exit.
#' @param relax Logical. If `TRUE`, selected structural validation problems are
#' stored as diagnostics rather than stopping execution.
#' @param quiet Logical. If `TRUE`, suppresses console output.
#'
#' @section Reported contribution percentages:
#' Detailed and grouped contribution tables include two percentage measures.
#' For explained components, `% Contribution Explained` reports each term's
#' contribution as a percentage of the explained gap. For unexplained
#' components, `% Contribution Unexplained` reports each term's contribution as
#' a percentage of the unexplained gap. `% Contribution Total` reports each
#' term's contribution as a percentage of the total observed gap. These
#' percentages are descriptive summaries derived from the estimated
#' contributions and do not change the underlying decomposition.
#'
#' @return
#' A named list with the following components:
#' \describe{
#'   \item{summary_stats}{Sample size, population size, design degrees of
#'   freedom, number of strata, number of PSUs, group-specific sample sizes, and
#'   the labels mapped to the favored (code 0) and disadvantaged (code 1)
#'   groups, plus the group-selection mode and counts of observations filtered
#'   from unselected levels.}
#'   \item{model_metrics}{Model fit statistics for group-specific and reference
#'   models.}
#'   \item{diagnostics}{Warnings and diagnostic messages generated during
#'   validation, estimation, or variance estimation.}
#'   \item{results_overall}{Overall decomposition table with group predictions,
#'   total difference, explained component, and unexplained component.}
#'   \item{results_detailed_explained}{Detailed explained contributions by
#'   predictor, including percentage contribution to the explained gap and to
#'   the total observed gap.}
#'   \item{results_detailed_unexplained}{Detailed unexplained contributions by
#'   predictor, including percentage contribution to the unexplained gap and to
#'   the total observed gap.}
#'   \item{results_grouped_explained}{Optional grouped explained contributions,
#'   including percentage contribution to the explained gap and to the total
#'   observed gap.}
#'   \item{results_grouped_unexplained}{Optional grouped unexplained
#'   contributions, including percentage contribution to the unexplained gap and
#'   to the total observed gap.}
#'   \item{models}{Fitted group-specific and reference model objects.}
#'   \item{raw}{Unrounded estimates, standard errors, covariance matrix,
#'   group-specific and reference coefficients, covariate means, analytic-sample
#'   row indices, fitted and decomposition model terms, normalization metadata,
#'   pooled-anchor metadata, replicate estimates/scales/centering diagnostics,
#'   group mapping (selected, available, and excluded levels), VCE method, and
#'   survey mode.}
#' }
#'
#' @examples
#' data(lunadecomp_example)
#'
#' fit_ols <- oby_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "y_continuous",
#'   group_var = "group",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   model_type = "ols",
#'   ref_method = "pooled",
#'   vce_method = "linearized",
#'   quiet = TRUE
#' )
#'
#' fit_ols$results_overall
#'
#' # Explicit order: first favored, then disadvantaged. If the original
#' # variable had additional levels, they would be filtered before estimation.
#' fit_ordered <- oby_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "y_continuous",
#'   group_var = "group",
#'   indep_vars = c("age", "rural"),
#'   group_levels = c(0, 1),
#'   quiet = TRUE
#' )
#'
#' fit_lpm <- oby_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "y_binary",
#'   group_var = "group",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   model_type = "lpm",
#'   ref_method = "pooled",
#'   vce_method = "linearized",
#'   quiet = TRUE
#' )
#'
#' fit_lpm$results_overall
#'
#' \donttest{
#' fit_logit <- oby_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "y_binary",
#'   group_var = "group",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   model_type = "logit",
#'   ref_method = "pooled",
#'   vce_method = "linearized",
#'   quiet = TRUE
#' )
#'
#' fit_logit$results_overall
#'
#' # For an adverse outcome, changing the anchor changes the counterfactual
#' # question, not the favored-minus-disadvantaged direction of the gap.
#' fit_logit_disadvantaged_anchor <- oby_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "y_binary",
#'   group_var = "group",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   model_type = "logit",
#'   ref_method = "pooled",
#'   pooled_anchor = "disadvantaged",
#'   vce_method = "linearized",
#'   quiet = TRUE
#' )
#' }
#'
#' @author
#' Adrian Vasquez-Mejia, MD, MSc \cr
#' Oscar J. Mujica, MD, MPH, PHE, FACE \cr
#' Antonio Sanhueza, MPH, MSc, PhD
#'
#' @references
#' Oaxaca, R. (1973). Male-female wage differentials in urban labor markets.
#' \emph{International Economic Review}, 14(3), 693.
#' \doi{10.2307/2525981}
#'
#' Blinder, A. S. (1973). Wage discrimination: Reduced form and structural
#' estimates. \emph{The Journal of Human Resources}, 8(4), 436.
#' \doi{10.2307/144855}
#'
#' Neumark, D. (1988). Employers' discriminatory behavior and the estimation
#' of wage discrimination. \emph{The Journal of Human Resources}, 23(3), 279.
#' \doi{10.2307/145830}
#'
#' Oaxaca, R. L., & Ransom, M. R. (1994). On discrimination and the
#' decomposition of wage differentials. \emph{Journal of Econometrics}, 61(1),
#' 5--21. \doi{10.1016/0304-4076(94)90074-4}
#'
#' Yun, M.-S. (2004). Decomposing differences in the first moment.
#' \emph{Economics Letters}, 82(2), 275--280.
#' \doi{10.1016/j.econlet.2003.09.008}
#'
#' Jann, B. (2008). The Blinder-Oaxaca decomposition for linear regression
#' models. \emph{The Stata Journal}, 8(4), 453--479.
#' \doi{10.1177/1536867X0800800401}
#'
#' Yun, M.-S. (2005). A simple solution to the identification problem in
#' detailed wage decompositions. \emph{Economic Inquiry}, 43(4), 766--772.
#' \doi{10.1093/ei/cbi053}
#'
#' @importFrom dplyr %>% select all_of mutate row_number group_by summarise n_distinct filter
#' @importFrom tidyr drop_na
#' @importFrom tibble tibble add_row
#' @importFrom rlang .data
#' @importFrom MASS ginv
#' @importFrom Matrix bdiag
#' @importFrom survey svydesign as.svrepdesign withReplicates
#' @importFrom stats as.formula complete.cases dlogis dnorm gaussian glm.fit median pnorm pt qt quasibinomial terms weighted.mean plogis setNames
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export


oby_decomp <- function(
    data, dep_var, group_var, indep_vars, groupings = list(),
    favored_group = NULL,
    ref_method = "pooled", model_type = "ols",
    use_svy = FALSE, weight_var = NULL, strata_var = NULL, psu_var = NULL,
    vce_method = "linearized", boot_reps = 500, lonely_psu = "adjust",
    level = 0.95, seed = NULL, relax = FALSE, quiet = FALSE,
    normalize = NULL, group_levels = NULL, pooled_anchor = "favored"
) {

  # ==============================================================================
  # 1. VALIDATION, DATA PREPARATION & SAFETY CHECKS
  # ==============================================================================
  valid_vce <- c("linearized", "jackknife", "bootstrap")
  valid_ref <- c("pooled", "neumark", "cotton", "reimers", "group1", "group2")
  valid_pooled_anchor <- c("favored", "disadvantaged", "centered")
  valid_model <- c("ols", "lpm", "logit", "probit")
  valid_lonely <- c("fail", "remove", "certainty", "adjust", "average")

  if (!(vce_method %in% valid_vce)) {
    stop(
      "Error: vce_method must be one of: ",
      paste(valid_vce, collapse = ", "),
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

  if (!(model_type %in% valid_model)) {
    stop(
      "Error: model_type must be one of: ",
      paste(valid_model, collapse = ", "),
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

  if (!is.data.frame(data)) {
    stop("Error: data must be a data frame.", call. = FALSE)
  }

  if (!is.character(dep_var) || length(dep_var) != 1) {
    stop("Error: dep_var must be a single character string.", call. = FALSE)
  }

  if (!is.character(group_var) || length(group_var) != 1) {
    stop("Error: group_var must be a single character string.", call. = FALSE)
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

  if (!is.character(indep_vars) || length(indep_vars) < 1) {
    stop("Error: indep_vars must be a non-empty character vector.", call. = FALSE)
  }

  form_base <- tryCatch(
    stats::as.formula(paste("~", paste(indep_vars, collapse = " + "))),
    error = function(e) {
      stop(
        "Error: indep_vars could not be parsed as R model terms: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  model_variables <- all.vars(form_base)
  if (length(model_variables) < 1) {
    stop(
      "Error: indep_vars must contain at least one predictor variable.",
      call. = FALSE
    )
  }

  if (
    !is.null(normalize) &&
      (!is.character(normalize) || anyNA(normalize) ||
        length(normalize) < 1 || any(!nzchar(normalize)))
  ) {
    stop(
      "Error: normalize must be NULL, \"all\", or a non-empty character vector of factor names.",
      call. = FALSE
    )
  }

  if (!is.null(normalize) && "all" %in% normalize && length(normalize) != 1) {
    stop(
      "Error: normalize = \"all\" must be used by itself, not combined with individual factor names.",
      call. = FALSE
    )
  }

  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("Error: level must be a numeric value between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(boot_reps) || length(boot_reps) != 1 || boot_reps < 2) {
    stop("Error: boot_reps must be a numeric value greater than or equal to 2.", call. = FALSE)
  }

  design_var_args <- list(
    weight_var = weight_var,
    strata_var = strata_var,
    psu_var = psu_var
  )
  invalid_design_args <- names(design_var_args)[vapply(
    design_var_args,
    function(x) {
      !is.null(x) &&
        (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x))
    },
    logical(1)
  )]
  if (length(invalid_design_args) > 0) {
    stop(
      "Error: ",
      paste(invalid_design_args, collapse = ", "),
      " must be NULL or a single non-empty character string.",
      call. = FALSE
    )
  }

  req_vars <- unique(c(
    dep_var,
    group_var,
    model_variables,
    weight_var,
    strata_var,
    psu_var
  ))
  req_vars <- req_vars[!is.na(req_vars) & nzchar(req_vars)]

  missing_vars <- setdiff(req_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "Error: the following variables are missing from data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(weight_var)) {
    weight_values <- data[[weight_var]]
    if (!is.numeric(weight_values)) {
      stop(
        "Error: weight_var must identify a numeric sampling-weight variable.",
        call. = FALSE
      )
    }
    if (any(is.infinite(weight_values))) {
      stop(
        "Error: sampling weights must be finite; infinite weights are not allowed.",
        call. = FALSE
      )
    }
    if (any(weight_values < 0, na.rm = TRUE)) {
      stop(
        "Error: sampling weights must be non-negative.",
        call. = FALSE
      )
    }
  }

  df_clean <- data %>%
    dplyr::mutate(.oby_source_row = dplyr::row_number()) %>%
    dplyr::select(dplyr::all_of(c(".oby_source_row", req_vars))) %>%
    tidyr::drop_na() %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    droplevels()
  n_missing_dropped <- nrow(data) - nrow(df_clean)

  factor_predictors <- model_variables[vapply(
    df_clean[model_variables],
    is.factor,
    logical(1)
  )]
  normalized_factors <- if (is.null(normalize)) {
    character(0)
  } else if (identical(normalize, "all")) {
    factor_predictors
  } else {
    unique(normalize)
  }

  invalid_normalize <- setdiff(normalized_factors, model_variables)
  if (length(invalid_normalize) > 0) {
    stop(
      "Error: normalize contains variable(s) not included in indep_vars: ",
      paste(invalid_normalize, collapse = ", "),
      call. = FALSE
    )
  }

  nonfactor_normalize <- setdiff(normalized_factors, factor_predictors)
  if (length(nonfactor_normalize) > 0) {
    stop(
      "Error: normalize can only include factor or character predictors: ",
      paste(nonfactor_normalize, collapse = ", "),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------------------------
  # GROUP CODING & FAVORED / DISADVANTAGED ASSIGNMENT
  # ----------------------------------------------------------------------------
  # The reported gap is always Difference = E(Y | favored) - E(Y | disadvantaged),
  # with favored coded 0 (minuend / reference structure) and disadvantaged coded 1.
  grp_levels_raw <- unique(as.character(df_clean[[group_var]]))
  grp_levels_raw <- grp_levels_raw[!is.na(grp_levels_raw)]
  available_group_levels <- grp_levels_raw
  excluded_group_levels <- character(0)
  n_group_filtered <- 0L

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
    keep_group <- as.character(df_clean[[group_var]]) %in%
      selected_group_levels
    n_group_filtered <- sum(!keep_group)
    excluded_group_levels <- setdiff(
      available_group_levels,
      selected_group_levels
    )
    df_clean <- droplevels(df_clean[keep_group, , drop = FALSE])
    df_clean[[group_var]] <- ifelse(
      as.character(df_clean[[group_var]]) == favored_chr,
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
        paste(grp_levels_raw, collapse = ", "),
        call. = FALSE
      )
    }
    disadvantaged_chr <- setdiff(grp_levels_raw, favored_chr)
    df_clean[[group_var]] <- ifelse(
      as.character(df_clean[[group_var]]) == favored_chr, 0, 1
    )
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
    # No explicit declaration: require 0/1 coding, group 0 is favored by convention.
    df_clean[[group_var]] <- suppressWarnings(as.numeric(as.character(df_clean[[group_var]])))
    if (any(is.na(df_clean[[group_var]]))) {
      stop(
        "Error: group_var is not coded 0/1. Either recode it as 0 and 1, or ",
        "supply group_levels = c(favored, disadvantaged).",
        call. = FALSE
      )
    }
    if (!setequal(unique(df_clean[[group_var]]), c(0, 1))) {
      stop(
        "Error: group_var must be coded as 0 and 1 (or use group_levels).",
        call. = FALSE
      )
    }
    favored_label <- "0"
    disadvantaged_label <- "1"
    group_selection_mode <- "implicit_0_1"
  }

  if (!all(c(0, 1) %in% unique(df_clean[[group_var]]))) {
    stop("Error: group_var must contain observations in both groups.", call. = FALSE)
  }

  if (model_type %in% c("lpm", "logit", "probit")) {
    df_clean[[dep_var]] <- as.numeric(as.character(df_clean[[dep_var]]))

    if (any(is.na(df_clean[[dep_var]]))) {
      stop(
        "Error: dep_var could not be converted safely to numeric values.",
        call. = FALSE
      )
    }

    unique_y <- unique(df_clean[[dep_var]])
    if (!all(unique_y %in% c(0, 1))) {
      stop(
        "Error: for lpm, logit, and probit models, dep_var must be coded as 0 and 1.",
        call. = FALSE
      )
    }
    if (length(unique_y) < 2) {
      stop(
        "Error: dep_var has no variation (all values are ", unique_y[1],
        "). A degenerate outcome cannot be decomposed.",
        call. = FALSE
      )
    }
  }


  N_obs <- nrow(df_clean)
  if (!quiet && n_missing_dropped > 0) {
    cat(sprintf(
      "\n[!] Data Prep: %d rows with NAs were dropped.\n",
      n_missing_dropped
    ))
  }
  if (!quiet && n_group_filtered > 0) {
    cat(sprintf(
      "[!] Group selection: %d rows outside group_levels were filtered. N = %d\n",
      n_group_filtered,
      N_obs
    ))
  }

  # A1: honor documented auto-activation of survey mode when weights are supplied.
  if (!use_svy && !is.null(weight_var)) {
    use_svy <- TRUE
    if (!quiet) message("Note: 'weight_var' supplied with use_svy = FALSE; survey mode activated automatically.")
  }

  df_clean <- df_clean %>%
    dplyr::mutate(
      svy_weight = if (use_svy && !is.null(weight_var)) as.numeric(.data[[weight_var]]) else 1,
      svy_strata = if (use_svy && !is.null(strata_var)) as.character(.data[[strata_var]]) else "1",
      svy_psu    = if (use_svy && !is.null(psu_var)) as.character(.data[[psu_var]]) else as.character(dplyr::row_number())
    )

  if (use_svy) {
    group_weight_totals <- tapply(
      df_clean$svy_weight,
      df_clean[[group_var]],
      sum
    )
    if (
      length(group_weight_totals) != 2 ||
        any(!is.finite(group_weight_totals)) ||
        any(group_weight_totals <= 0)
    ) {
      stop(
        "Error: each selected group must have positive total sampling weight.",
        call. = FALSE
      )
    }
  }

  idx_0 <- which(df_clean[[group_var]] == 0); idx_1 <- which(df_clean[[group_var]] == 1)
  N_0 <- length(idx_0); N_1 <- length(idx_1)
  pop_size <- sum(df_clean$svy_weight)
  n_strata <- if(use_svy) length(unique(df_clean$svy_strata)) else NA
  # A2: count PSUs within strata (nested design), not globally, so re-used PSU
  # labels across strata are not collapsed.
  n_psu    <- if(use_svy) nrow(dplyr::distinct(df_clean, .data$svy_strata, .data$svy_psu)) else NA
  psu_per_stratum <- if (use_svy) {
    table(
      unique(df_clean[c("svy_strata", "svy_psu")])$svy_strata
    )
  } else {
    integer(0)
  }
  if (
    use_svy &&
      vce_method == "bootstrap" &&
      any(psu_per_stratum < 2)
  ) {
    stop(
      "Error: Rao-Wu survey bootstrap requires at least two PSUs in ",
      "every stratum. Collapse or redesign singleton strata, or use an ",
      "appropriate set of externally supplied replicate weights.",
      call. = FALSE
    )
  }
  # design_df is computed after k (number of model parameters) is known; see below.

  # ------------------------------------------------------------------------------
  # DYNAMIC LABELS DICTIONARY: Maps dummies to informative strings and tags numeric
  # ------------------------------------------------------------------------------
  display_dict <- list()
  for (var in model_variables) {
    if (is.factor(df_clean[[var]])) {
      lvls <- levels(df_clean[[var]])
      ref_lvl <- lvls[1]
      displayed_levels <- if (var %in% normalized_factors) lvls else lvls[-1]
      for (l in displayed_levels) {
        dummy_name <- paste0(var, l)
        display_dict[[dummy_name]] <- if (var %in% normalized_factors) {
          sprintf("%s (%s; normalized)", var, l)
        } else {
          sprintf("%s (%s vs. Ref: %s)", var, l, ref_lvl)
        }
      }
    } else if (is.numeric(df_clean[[var]])) {
      display_dict[[var]] <- sprintf("%s (Continuous)", var)
    } else {
      display_dict[[var]] <- var
    }
  }

  # ==============================================================================
  # 2. EXACT ALIGNMENT OF PREDICTOR MATRICES & STRUCTURAL VALIDATION
  # ==============================================================================
  factor_contrasts <- stats::setNames(
    lapply(
      factor_predictors,
      function(v) stats::contr.treatment(levels(df_clean[[v]]), base = 1)
    ),
    factor_predictors
  )
  X_full <- model.matrix(
    form_base,
    data = df_clean,
    contrasts.arg = factor_contrasts
  )
  Y_full <- if (is.factor(df_clean[[dep_var]])) as.numeric(df_clean[[dep_var]]) - 1 else as.numeric(df_clean[[dep_var]])

  # Build the post-estimation deviation-contrast representation used for the
  # detailed decomposition. The fitted models remain in their identified K-1
  # treatment-coded form; selected factors are expanded to all K categories
  # only after fitting, exactly as in Jann's normalize().
  fit_model_terms <- colnames(X_full)
  fit_assign <- attr(X_full, "assign")
  fit_term_labels <- attr(stats::terms(form_base), "term.labels")
  fit_term_columns <- stats::setNames(
    lapply(seq_along(fit_term_labels), function(i) {
      fit_model_terms[fit_assign == i]
    }),
    fit_term_labels
  )
  term_variables <- stats::setNames(
    lapply(
      fit_term_labels,
      function(term) {
        all.vars(stats::as.formula(paste("~", term)))
      }
    ),
    fit_term_labels
  )

  term_column_map <- fit_term_columns
  normalization_blocks <- list()
  normalized_term_specs <- list()

  # Main categorical effects are centered around the intercept.
  for (var in normalized_factors) {
    if (!(var %in% fit_term_labels)) {
      stop(
        "Error: normalized factor '", var,
        "' must enter the model with its main effect. ",
        "Use a hierarchical specification such as '", var,
        " * continuous_variable'.",
        call. = FALSE
      )
    }

    lvls <- levels(df_clean[[var]])
    full_cols <- paste0(var, lvls)
    old_cols <- fit_term_columns[[var]]
    expected_old_cols <- full_cols[-1]
    if (!identical(old_cols, expected_old_cols)) {
      stop(
        "Error: factor '", var,
        "' could not be mapped to treatment-coded indicators for normalization.",
        call. = FALSE
      )
    }

    term_column_map[[var]] <- full_cols
    normalization_blocks[[var]] <- list(
      levels = lvls,
      omitted_level = lvls[[1]],
      fit_columns = old_cols,
      normalized_columns = full_cols,
      interactions = list()
    )
    normalized_term_specs[[var]] <- list(
      kind = "main",
      factor = var,
      anchor_fit = "(Intercept)",
      anchor_decomp = "(Intercept)",
      old_columns = old_cols,
      new_columns = full_cols
    )
  }

  # A factor-by-continuous interaction is centered independently around the
  # corresponding continuous main-effect coefficient, mirroring
  # normalize(d1x ... dKx # x) in Jann's oaxaca.
  for (term in fit_term_labels) {
    vars_in_term <- term_variables[[term]]
    normalized_in_term <- intersect(vars_in_term, normalized_factors)
    if (length(vars_in_term) < 2 || length(normalized_in_term) == 0) {
      next
    }
    if (length(vars_in_term) != 2 || length(normalized_in_term) != 1) {
      stop(
        "Error: normalization currently supports interactions between one ",
        "normalized factor and one continuous variable. Unsupported term: ",
        term,
        call. = FALSE
      )
    }

    factor_var <- normalized_in_term[[1]]
    continuous_var <- setdiff(vars_in_term, factor_var)
    if (
      length(continuous_var) != 1 ||
        !is.numeric(df_clean[[continuous_var]])
    ) {
      stop(
        "Error: normalized interaction '", term,
        "' must pair factor '", factor_var,
        "' with one continuous numeric variable.",
        call. = FALSE
      )
    }
    if (!(continuous_var %in% fit_term_labels)) {
      stop(
        "Error: normalized interaction '", term,
        "' requires the main effect of continuous variable '",
        continuous_var, "'. Use a hierarchical model specification.",
        call. = FALSE
      )
    }
    anchor_columns <- fit_term_columns[[continuous_var]]
    if (
      length(anchor_columns) != 1 ||
        !identical(anchor_columns, continuous_var)
    ) {
      stop(
        "Error: continuous anchor '", continuous_var,
        "' could not be mapped to a single main-effect coefficient.",
        call. = FALSE
      )
    }

    factor_columns <- normalization_blocks[[factor_var]]$normalized_columns
    nonbase_factor_columns <- factor_columns[-1]
    old_cols <- fit_term_columns[[term]]
    factor_first <- paste0(
      nonbase_factor_columns,
      ":",
      continuous_var
    )
    continuous_first <- paste0(
      continuous_var,
      ":",
      nonbase_factor_columns
    )
    ordered_old_cols <- if (all(factor_first %in% old_cols)) {
      factor_first
    } else if (all(continuous_first %in% old_cols)) {
      continuous_first
    } else {
      stop(
        "Error: interaction '", term,
        "' could not be mapped to the treatment-coded levels of factor '",
        factor_var, "'.",
        call. = FALSE
      )
    }
    new_cols <- paste0(
      factor_columns,
      ":",
      continuous_var
    )

    term_column_map[[term]] <- new_cols
    normalized_term_specs[[term]] <- list(
      kind = "interaction",
      factor = factor_var,
      continuous = continuous_var,
      anchor_fit = continuous_var,
      anchor_decomp = continuous_var,
      old_columns = ordered_old_cols,
      new_columns = new_cols
    )
    normalization_blocks[[factor_var]]$interactions[[term]] <- list(
      continuous = continuous_var,
      fit_columns = ordered_old_cols,
      normalized_columns = new_cols
    )
    for (j in seq_along(new_cols)) {
      display_dict[[new_cols[[j]]]] <- sprintf(
        "%s (%s; normalized) \u00d7 %s",
        factor_var,
        normalization_blocks[[factor_var]]$levels[[j]],
        continuous_var
      )
    }
  }

  decomp_terms <- c(
    "(Intercept)",
    unlist(term_column_map[fit_term_labels], use.names = FALSE)
  )
  if (anyDuplicated(decomp_terms)) {
    stop(
      "Error: normalization generated duplicate model-term names. ",
      "Rename the affected predictor or factor levels.",
      call. = FALSE
    )
  }

  coefficient_transform <- matrix(
    0,
    nrow = length(fit_model_terms),
    ncol = length(decomp_terms),
    dimnames = list(fit_model_terms, decomp_terms)
  )
  mean_transform <- coefficient_transform
  coefficient_transform["(Intercept)", "(Intercept)"] <- 1
  mean_transform["(Intercept)", "(Intercept)"] <- 1

  for (term in fit_term_labels) {
    spec <- normalized_term_specs[[term]]
    old_cols <- fit_term_columns[[term]]
    if (is.null(spec)) {
      coefficient_transform[cbind(old_cols, old_cols)] <- 1
      mean_transform[cbind(old_cols, old_cols)] <- 1
      next
    }

    old_cols <- spec$old_columns
    new_cols <- spec$new_columns
    K <- length(new_cols)

    # Coefficients: insert the omitted coefficient as zero, subtract the
    # equal-category mean from every level, and add that mean to the relevant
    # anchor: the intercept for a main effect or the continuous coefficient for
    # an interaction.
    coefficient_transform[old_cols, spec$anchor_decomp] <-
      coefficient_transform[old_cols, spec$anchor_decomp] + 1 / K
    coefficient_transform[old_cols, new_cols] <- -1 / K
    for (j in seq_along(old_cols)) {
      coefficient_transform[old_cols[[j]], new_cols[[j + 1]]] <-
        coefficient_transform[old_cols[[j]], new_cols[[j + 1]]] + 1
    }

    # Means: reconstruct the omitted column from its anchor minus the remaining
    # K-1 columns. For main effects the anchor is one; for interactions it is
    # the continuous variable itself.
    mean_transform[spec$anchor_fit, new_cols[[1]]] <- 1
    mean_transform[old_cols, new_cols[[1]]] <- -1
    for (j in seq_along(old_cols)) {
      mean_transform[old_cols[[j]], new_cols[[j + 1]]] <- 1
    }
  }

  X_decomp <- X_full %*% mean_transform
  colnames(X_decomp) <- decomp_terms
  k_decomp <- ncol(X_decomp)

  pooled_group_mean <- sum(
    df_clean[[group_var]] * df_clean$svy_weight
  ) / sum(df_clean$svy_weight)
  pooled_anchor_value <- switch(
    pooled_anchor,
    "favored" = 0,
    "disadvantaged" = 1,
    "centered" = pooled_group_mean
  )
  pooled_group_column <- df_clean[[group_var]] - pooled_anchor_value
  pooled_group_column_name <- ".oby_pooled_group"
  while (pooled_group_column_name %in% colnames(X_full)) {
    pooled_group_column_name <- paste0(
      pooled_group_column_name,
      "_"
    )
  }
  pooled_anchor_description <- switch(
    pooled_anchor,
    "favored" = "favored group (D = 0)",
    "disadvantaged" = "disadvantaged group (D = 1)",
    "centered" = sprintf(
      "weighted group mean (D = %.17g)",
      pooled_group_mean
    )
  )
  pooled_anchor_info <- list(
    requested = pooled_anchor,
    applicable = identical(ref_method, "pooled"),
    effective = if (identical(ref_method, "pooled")) {
      pooled_anchor
    } else {
      NA_character_
    },
    material_nonlinear = identical(ref_method, "pooled") &&
      model_type %in% c("logit", "probit"),
    reference_value = if (identical(ref_method, "pooled")) {
      pooled_anchor_value
    } else {
      NA_real_
    },
    group_mean = pooled_group_mean,
    model_term = if (identical(ref_method, "pooled")) {
      pooled_group_column_name
    } else {
      NA_character_
    },
    description = if (identical(ref_method, "pooled")) {
      pooled_anchor_description
    } else {
      "not applicable because ref_method is not pooled"
    }
  )

  if (ref_method %in% c("pooled", "neumark")) {
    X_pool <- if (ref_method == "pooled") {
      cbind(
        X_full,
        stats::setNames(
          data.frame(pooled_group_column),
          pooled_group_column_name
        )
      )
    } else {
      X_full
    }
    X_pool <- as.matrix(X_pool)
  }

  X_0 <- X_full[idx_0, , drop = FALSE]; Y_0 <- Y_full[idx_0]
  X_1 <- X_full[idx_1, , drop = FALSE]; Y_1 <- Y_full[idx_1]
  X_decomp_0 <- X_decomp[idx_0, , drop = FALSE]
  X_decomp_1 <- X_decomp[idx_1, , drop = FALSE]
  w_0 <- df_clean$svy_weight[idx_0]; w_1 <- df_clean$svy_weight[idx_1]; w_all <- df_clean$svy_weight
  k <- ncol(X_full)
  # Stata's ordinary linearized and prefix-replication results use asymptotic
  # normal inference. Survey results use the design degrees of freedom.
  design_df <- if (!use_svy && vce_method == "jackknife") {
    max(1, N_obs - 1)
  } else if (!use_svy) {
    Inf
  } else if (use_svy) {
    max(1, n_psu - n_strata)
  }

  family_obj <- switch(
    model_type,
    "ols" = gaussian(),
    "lpm" = gaussian(),
    "logit" = quasibinomial(link = "logit"),
    "probit" = quasibinomial(link = "probit")
  )
  # Nonlinear decompositions are sensitive to small coefficient differences
  # because Yun's detailed weights are ratios of linear-index contributions.
  # A tight convergence criterion also aligns the MLEs with Stata's logit and
  # probit estimates instead of stopping at glm.fit()'s comparatively loose
  # default deviance tolerance.
  fit_control <- stats::glm.control(
    epsilon = if (model_type %in% c("logit", "probit")) 1e-16 else 1e-8,
    maxit = if (model_type %in% c("logit", "probit")) 100 else 25
  )

  diagnostics <- list()

  # Group-specific models must be identified: each group needs more observations
  # than estimated parameters (aliased handling exists, but an unidentified
  # group model makes the decomposition meaningless).
  if (k >= min(N_0, N_1)) {
    diagnostics$group_df <- sprintf("CRITICAL MODEL VALIDATION: Predictors (%d) exceeds or equals the size of the smallest group (%d). Group-specific model is unidentified.", k, min(N_0, N_1))
    if (!relax) stop(diagnostics$group_df) else if (!quiet) warning(diagnostics$group_df, call. = FALSE)
  }

  # A4: standard errors from the covariance diagonal without masking non-positive
  # variances (a symptom of a non-PSD covariance) via abs(); such terms are
  # flagged and returned as NA instead of a fabricated SE.
  safe_se <- function(V) {
    vd <- diag(V)
    neg <- which(vd < -1e-8)
    if (length(neg) > 0) {
      diagnostics$neg_var <<- sprintf(
        "Numerical Warning: %d non-positive variance term(s) detected (non-PSD covariance); their SE set to NA.",
        length(neg)
      )
      if (!relax && !quiet) warning(diagnostics$neg_var)
    }
    vd[vd < 0] <- NA_real_
    sqrt(vd)
  }

  # Scalar version for aggregated (grouped) variances: same policy as safe_se -
  # a negative variance sum (non-PSD covariance) yields NA, never sqrt(abs()).
  safe_se_val <- function(v) {
    if (is.na(v)) return(NA_real_)
    if (v < -1e-8) {
      diagnostics$neg_var_grp <<- "Numerical Warning: negative grouped variance sum detected (non-PSD covariance); grouped SE set to NA."
      if (!relax && !quiet) warning(diagnostics$neg_var_grp)
      return(NA_real_)
    }
    sqrt(max(v, 0))
  }

  # --- PLAUSIBILITY SHIELD 1: Extreme Leverage Weights ---
  if (use_svy) {
    max_w <- max(w_all, na.rm = TRUE)
    med_w <- median(w_all, na.rm = TRUE)
    if (med_w > 0 && (max_w > 50 * med_w)) {
      diagnostics$extreme_weights <- sprintf("Leverage Warning: Extreme sampling weights detected (Max is %.0f times the median).", max_w/med_w)
      if (!relax && !quiet) warning(diagnostics$extreme_weights)
    }
  }

  # --- PLAUSIBILITY SHIELD 2: Perfect Collinearity ---
  if (base::qr(X_full)$rank < ncol(X_full)) {
    diagnostics$collinearity <- "CRITICAL: Perfect collinearity detected in the global model matrix."
    if (!relax) stop(diagnostics$collinearity)
  }

  if (use_svy) {
    lonely_count <- df_clean %>%
      dplyr::group_by(.data$svy_strata) %>%
      dplyr::summarise(
        n = dplyr::n_distinct(.data$svy_psu),
        .groups = "drop"
      ) %>%
      dplyr::filter(.data$n == 1) %>%
      nrow()
    if (lonely_count > 0 && !quiet) cat(sprintf("Warning: %d strata contain a single PSU. Applied adjustment: '%s'.\n", lonely_count, lonely_psu))
  }

  # ==============================================================================
  # 3. BASE MODEL ESTIMATION & COMMON SUPPORT AUDIT
  # ==============================================================================
  # Fast GLM fitting on pre-aligned matrices
  mod_0 <- suppressWarnings(glm.fit(
    x = X_0, y = Y_0, weights = w_0 / mean(w_0),
    family = family_obj, control = fit_control
  ))
  mod_1 <- suppressWarnings(glm.fit(
    x = X_1, y = Y_1, weights = w_1 / mean(w_1),
    family = family_obj, control = fit_control
  ))

  b_0 <- mod_0$coefficients; b_1 <- mod_1$coefficients

  # B2: flag non-convergence (glm.fit warnings are suppressed above).
  if (isFALSE(mod_0$converged) || isFALSE(mod_1$converged)) {
    diagnostics$convergence <- "Warning: a group-specific GLM did not converge (possible separation or sparse data)."
    if (!relax && !quiet) warning(diagnostics$convergence)
  }

  # B4: flag aliased (within-group collinear) coefficients instead of silently zeroing.
  na0 <- is.na(b_0); na1 <- is.na(b_1)
  if (any(na0) || any(na1)) {
    diagnostics$aliased <- sprintf(
      "Note: %d coefficient(s) aliased (collinear within a group) and set to 0.", sum(na0) + sum(na1)
    )
    if (!quiet) message(diagnostics$aliased)
  }
  b_0[na0] <- 0; b_1[na1] <- 0

  if (model_type %in% c("logit", "probit")) {
    b0_clean <- b_0[names(b_0) != "(Intercept)"]
    b1_clean <- b_1[names(b_1) != "(Intercept)"]
    max_beta <- max(abs(c(b0_clean, b1_clean)), na.rm = TRUE)

    # B3: probit coefficients live on a tighter scale than logit ones.
    beta_thr <- if (model_type == "probit") 3 else 5
    if (!is.infinite(max_beta) && max_beta > beta_thr) {
      diagnostics$coef_explosion <- sprintf("Plausibility Warning: Extremely large coefficients detected (|beta| Max = %.1f, threshold = %g). Risk of Quasi-Complete Separation due to sparse data.", max_beta, beta_thr)
      if (!relax && !quiet) warning(diagnostics$coef_explosion)
    }
  }

  mcfadden_r2 <- function(m) 1 - (m$deviance / m$null.deviance)
  fit_metrics <- tibble::tibble(Model = c("Group 0", "Group 1"), Pseudo_R2 = c(mcfadden_r2(mod_0), mcfadden_r2(mod_1)), Deviance = c(mod_0$deviance, mod_1$deviance))
  models_list <- list(group_0 = mod_0, group_1 = mod_1)

  if (ref_method %in% c("pooled", "neumark")) {
    mod_ref <- suppressWarnings(glm.fit(
      x = X_pool, y = Y_full, weights = w_all / mean(w_all),
      family = family_obj, control = fit_control
    ))
    b_star <- mod_ref$coefficients[colnames(X_full)]
    b_star[is.na(b_star)] <- 0

    fit_metrics <- fit_metrics %>% tibble::add_row(Model = paste("Reference (", ref_method, ")"), Pseudo_R2 = mcfadden_r2(mod_ref), Deviance = mod_ref$deviance)
    models_list$reference <- mod_ref
  } else {
    W_ref <- switch(ref_method, "group1" = 1, "group2" = 0, "reimers" = 0.5, "cotton" = sum(w_0)/sum(w_all))
    b_star <- W_ref * b_0 + (1 - W_ref) * b_1
  }

  names(b_0) <- names(b_1) <- names(b_star) <- colnames(X_full)
  x0_m <- as.numeric(colSums(X_0 * w_0) / sum(w_0))
  x1_m <- as.numeric(colSums(X_1 * w_1) / sum(w_1))
  names(x0_m) <- names(x1_m) <- colnames(X_full)
  base_estimates_fit <- c(
    stats::setNames(b_0, paste0("b0:", colnames(X_full))),
    stats::setNames(b_1, paste0("b1:", colnames(X_full))),
    stats::setNames(b_star, paste0("b_ref:", colnames(X_full))),
    stats::setNames(x0_m, paste0("x0:", colnames(X_full))),
    stats::setNames(x1_m, paste0("x1:", colnames(X_full)))
  )

  b_0_decomp <- as.numeric(b_0 %*% coefficient_transform)
  b_1_decomp <- as.numeric(b_1 %*% coefficient_transform)
  b_star_decomp <- as.numeric(b_star %*% coefficient_transform)
  x0_decomp <- as.numeric(x0_m %*% mean_transform)
  x1_decomp <- as.numeric(x1_m %*% mean_transform)
  names(b_0_decomp) <- names(b_1_decomp) <-
    names(b_star_decomp) <- decomp_terms
  names(x0_decomp) <- names(x1_decomp) <- decomp_terms

  base_estimates <- c(
    stats::setNames(b_0_decomp, paste0("b0:", decomp_terms)),
    stats::setNames(b_1_decomp, paste0("b1:", decomp_terms)),
    stats::setNames(b_star_decomp, paste0("b_ref:", decomp_terms)),
    stats::setNames(x0_decomp, paste0("x0:", decomp_terms)),
    stats::setNames(x1_decomp, paste0("x1:", decomp_terms))
  )

  expected_names <- c(
    "Group_0", "Group_1", "Difference", "Explained", "Unexplained",
    paste0("Exp_", decomp_terms),
    paste0("Unexp_", decomp_terms)
  )

  # ==============================================================================
  # 4. VARIANCE ESTIMATION ENGINES
  # ==============================================================================
  replication_details <- NULL
  if (vce_method == "linearized") {

    likelihood_parts <- function(m, w) {
      prior_weights <- as.numeric(m$prior.weights)
      response <- as.numeric(m$y)
      eta <- as.numeric(m$linear.predictors)
      mu <- as.numeric(m$fitted.values)
      if (model_type %in% c("ols", "lpm")) {
        score_values <- prior_weights * (response - mu)
        hessian_weights <- prior_weights
      } else {
        mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
        variance <- mu * (1 - mu)
        q <- if (model_type == "logit") variance else stats::dnorm(eta)
        response_residual <- response - mu
        score_values <- prior_weights * response_residual * q / variance
        hessian_weights <- prior_weights * q^2 / variance
      }
      if (model_type == "probit") {
        # glm.fit() uses Fisher scoring. Stata's probit sandwich instead uses
        # the observed-information Hessian (OIM), which differs for this
        # non-canonical link. For eta = X beta, q = d Phi(eta)/d eta:
        # -d2l/deta2 = q^2/v - (y-mu)q'/v
        #               + (y-mu)q^2(1-2mu)/v^2.
        q_prime <- -eta * q
        hessian_weights <- prior_weights * (
          q^2 / variance -
            response_residual * q_prime / variance +
            response_residual * q^2 * (1 - 2 * mu) / variance^2
        )
      }
      list(
        score_values = score_values * mean(w),
        hessian_weights = hessian_weights
      )
    }
    get_svy_scores <- function(m, X, w) {
      X * likelihood_parts(m, w)$score_values
    }
    get_svy_bread <- function(m, X, w) {
      hessian_weights <- likelihood_parts(m, w)$hessian_weights
      MASS::ginv(crossprod(X, X * hessian_weights)) / mean(w)
    }
    # B1: honor all documented lonely_psu options in the linearized meat.
    #  - "fail"      : stop on a single-PSU stratum
    #  - "remove"/"certainty": drop its (zero) between-PSU contribution
    #  - "adjust"    : center the single PSU on the grand mean
    #  - "average"   : impute the average contribution of non-lonely strata
    get_svy_meat <- function(U, strata, psu) {
      Meat <- matrix(0, ncol(U), ncol(U))
      U_mean <- colMeans(U)
      regular_stratum_meats <- list()
      average_lonely_count <- 0L
      for (h in unique(strata)) {
        idx <- which(strata == h); zh <- rowsum(U[idx, , drop = FALSE], group = psu[idx]); nh <- nrow(zh)
        if (nh > 1) {
          stratum_meat <- (nh / (nh - 1)) *
            crossprod(sweep(zh, 2, colMeans(zh), "-"))
          Meat <- Meat + stratum_meat
          regular_stratum_meats[[length(regular_stratum_meats) + 1L]] <-
            stratum_meat
        } else {
          if (lonely_psu == "fail") {
            stop(sprintf("Stratum '%s' contains a single PSU (lonely_psu = 'fail').", h), call. = FALSE)
          }
          if (lonely_psu %in% c("remove", "certainty")) next
          if (lonely_psu == "adjust") {
            Meat <- Meat + crossprod(sweep(zh, 2, U_mean, "-"))
          }
          if (lonely_psu == "average") {
            average_lonely_count <- average_lonely_count + 1L
          }
        }
      }
      if (average_lonely_count > 0L) {
        if (length(regular_stratum_meats) == 0L) {
          stop(
            "Cannot apply lonely_psu = 'average': no non-lonely stratum is available.",
            call. = FALSE
          )
        }
        average_meat <- Reduce(`+`, regular_stratum_meats) /
          length(regular_stratum_meats)
        Meat <- Meat + average_lonely_count * average_meat
      }
      return(Meat)
    }
    pad_U <- function(U_sub, idx) { U <- matrix(0, nrow(df_clean), ncol(U_sub)); U[idx, ] <- U_sub; return(U) }

    if (ref_method %in% c("pooled", "neumark")) {
      Bread_mod <- as.matrix(Matrix::bdiag(get_svy_bread(mod_0,X_0,w_0), get_svy_bread(mod_1,X_1,w_1), get_svy_bread(mod_ref,X_pool,w_all)))
      U_mod <- cbind(pad_U(get_svy_scores(mod_0,X_0,w_0), idx_0), pad_U(get_svy_scores(mod_1,X_1,w_1), idx_1), get_svy_scores(mod_ref,X_pool,w_all))
      Meat_mod <- if(use_svy) get_svy_meat(U_mod, df_clean$svy_strata, df_clean$svy_psu) else crossprod(U_mod)

      V_mod <- Bread_mod %*% Meat_mod %*% t(Bread_mod)
      if (!use_svy) {
        # Match Stata suest's small-sample adjustment for the joint robust
        # covariance used by oaxaca with pooled/omega references.
        V_mod <- (N_obs / (N_obs - 1)) * V_mod
      }
      idx_bs_in_V <- 2 * k + match(colnames(X_full), colnames(X_pool))
      V0_coef <- V_mod[c(1:(2*k), idx_bs_in_V), c(1:(2*k), idx_bs_in_V)]
    } else {
      Bread_mod <- as.matrix(Matrix::bdiag(get_svy_bread(mod_0,X_0,w_0), get_svy_bread(mod_1,X_1,w_1)))
      U_mod <- cbind(pad_U(get_svy_scores(mod_0,X_0,w_0), idx_0), pad_U(get_svy_scores(mod_1,X_1,w_1), idx_1))
      Meat_mod <- if(use_svy) get_svy_meat(U_mod, df_clean$svy_strata, df_clean$svy_psu) else crossprod(U_mod)

      T_mat <- rbind(cbind(diag(k), matrix(0,k,k)), cbind(matrix(0,k,k), diag(k)), cbind(W_ref*diag(k), (1-W_ref)*diag(k)))
      V_group_coef <- Bread_mod %*% Meat_mod %*% t(Bread_mod)
      if (!use_svy) {
        # Stata uses different finite-sample corrections for the native robust
        # covariance: regress uses N_g/(N_g-k), whereas logit/probit use
        # N_g/(N_g-1). Apply the correction equation by equation before
        # constructing a weighted reference coefficient vector.
        adjustment_denominator_0 <- if (
          model_type %in% c("ols", "lpm")
        ) N_0 - k else N_0 - 1
        adjustment_denominator_1 <- if (
          model_type %in% c("ols", "lpm")
        ) N_1 - k else N_1 - 1
        coef_scale <- diag(c(
          rep(sqrt(N_0 / adjustment_denominator_0), k),
          rep(sqrt(N_1 / adjustment_denominator_1), k)
        ))
        V_group_coef <- coef_scale %*% V_group_coef %*% coef_scale
      }
      V0_coef <- T_mat %*% V_group_coef %*% t(T_mat)
    }

    UX0_pad <- pad_U(sweep(X_0, 2, x0_m, "-") * w_0, idx_0); UX1_pad <- pad_U(sweep(X_1, 2, x1_m, "-") * w_1, idx_1)
    Bread_x <- as.matrix(Matrix::bdiag(diag(1/sum(w_0), k), diag(1/sum(w_1), k)))
    Meat_x <- if(use_svy) get_svy_meat(cbind(UX0_pad, UX1_pad), df_clean$svy_strata, df_clean$svy_psu) else crossprod(cbind(UX0_pad, UX1_pad))
    V_x <- Bread_x %*% Meat_x %*% t(Bread_x)
    if (!use_svy) {
      # Match Stata mean's finite-sample covariance adjustment.
      mean_scale <- diag(c(
        rep(sqrt(N_0 / (N_0 - 1)), k),
        rep(sqrt(N_1 / (N_1 - 1)), k)
      ))
      V_x <- mean_scale %*% V_x %*% mean_scale
    }
    V_total_fit <- as.matrix(Matrix::bdiag(V0_coef, V_x))
    rownames(V_total_fit) <- colnames(V_total_fit) <-
      names(base_estimates_fit)
    joint_transform <- as.matrix(Matrix::bdiag(
      coefficient_transform,
      coefficient_transform,
      coefficient_transform,
      mean_transform,
      mean_transform
    ))
    V_total <- t(joint_transform) %*% V_total_fit %*% joint_transform
    rownames(V_total) <- colnames(V_total) <- names(base_estimates)

    cdf_f <- switch(
      model_type,
      "ols" = function(z) z,
      "lpm" = function(z) z,
      "logit" = function(z) 1 / (1 + exp(-z)),
      "probit" = function(z) stats::pnorm(z)
    )
    pdf_f <- switch(
      model_type,
      "ols" = function(z) rep(1, length(z)),
      "lpm" = function(z) rep(1, length(z)),
      # Numerically stable logistic density (exp(z)/(1+exp(z))^2 overflows to NaN for z > ~709)
      "logit" = function(z) stats::dlogis(z),
      "probit" = function(z) stats::dnorm(z)
    )

    xb0 <- as.vector(X_decomp_0 %*% b_0_decomp)
    xb1 <- as.vector(X_decomp_1 %*% b_1_decomp)
    xb0_s <- as.vector(X_decomp_0 %*% b_star_decomp)
    xb1_s <- as.vector(X_decomp_1 %*% b_star_decomp)

    P0 <- sum(cdf_f(xb0)*w_0)/sum(w_0); P1 <- sum(cdf_f(xb1)*w_1)/sum(w_1)
    P0s <- sum(cdf_f(xb0_s)*w_0)/sum(w_0); P1s <- sum(cdf_f(xb1_s)*w_1)/sum(w_1)
    est_ov <- c("Group_0"=P0, "Group_1"=P1, "Difference"=P0-P1, "Explained"=P0s-P1s, "Unexplained"=(P0-P0s)+(P1s-P1))

    b0_n <- as.numeric(b_0_decomp)
    b1_n <- as.numeric(b_1_decomp)
    bs_n <- as.numeric(b_star_decomp)
    WEn <- (x0_decomp - x1_decomp) * bs_n
    WE <- if(abs(sum(WEn))>1e-15) {
      WEn / sum(WEn)
    } else {
      rep(0, k_decomp)
    }
    WUn <- x0_decomp * (b0_n - bs_n) +
      x1_decomp * (bs_n - b1_n)
    WU <- if(abs(sum(WUn))>1e-15) {
      WUn / sum(WUn)
    } else {
      rep(0, k_decomp)
    }

    est_det_E <- est_ov["Explained"] * WE; est_det_U <- est_ov["Unexplained"] * WU

    fx0b0 <- as.numeric(colSums(
      X_decomp_0 * as.numeric(pdf_f(xb0)) * w_0
    ) / sum(w_0))
    fx1b1 <- as.numeric(colSums(
      X_decomp_1 * as.numeric(pdf_f(xb1)) * w_1
    ) / sum(w_1))
    fx0B <- as.numeric(colSums(
      X_decomp_0 * as.numeric(pdf_f(xb0_s)) * w_0
    ) / sum(w_0))
    fx1B <- as.numeric(colSums(
      X_decomp_1 * as.numeric(pdf_f(xb1_s)) * w_1
    ) / sum(w_1))
    fx0b0_m <- sum(pdf_f(xb0)*w_0)/sum(w_0); fx1b1_m <- sum(pdf_f(xb1)*w_1)/sum(w_1)
    fx0B_m  <- sum(pdf_f(xb0_s)*w_0)/sum(w_0); fx1B_m  <- sum(pdf_f(xb1_s)*w_1)/sum(w_1)

    z_k <- rep(0, k_decomp)
    G_ov <- rbind(
      c(fx0b0, z_k, z_k, fx0b0_m * b0_n, z_k),
      c(z_k, fx1b1, z_k, z_k, fx1b1_m * b1_n),
      c(fx0b0, -fx1b1, z_k, fx0b0_m * b0_n, -fx1b1_m * b1_n),
      c(z_k, z_k, fx0B - fx1B, fx0B_m * bs_n, -fx1B_m * bs_n),
      c(fx0b0, -fx1b1, -fx0B + fx1B, fx0b0_m * b0_n - fx0B_m * bs_n, fx1B_m * bs_n - fx1b1_m * b1_n)
    )

    SE_val <- sum(WEn); SU_val <- sum(WUn)
    dWE <- matrix(
      0,
      nrow = k_decomp,
      ncol = 5 * k_decomp
    )
    dWU <- matrix(
      0,
      nrow = k_decomp,
      ncol = 5 * k_decomp
    )

    if(abs(SE_val) > 1e-15) {
      dWE[, (2*k_decomp+1):(3*k_decomp)] <-
        diag((x0_decomp - x1_decomp)/SE_val) -
        outer(WE, x0_decomp - x1_decomp)/SE_val
      dWE[, (3*k_decomp+1):(4*k_decomp)] <-
        diag(bs_n/SE_val) - outer(WE, bs_n)/SE_val
      dWE[, (4*k_decomp+1):(5*k_decomp)] <-
        diag(-bs_n/SE_val) - outer(WE, -bs_n)/SE_val
    }
    if(abs(SU_val) > 1e-15) {
      dWU[, 1:k_decomp] <-
        diag(x0_decomp/SU_val) - outer(WU, x0_decomp)/SU_val
      dWU[, (k_decomp+1):(2*k_decomp)] <-
        diag(-x1_decomp/SU_val) - outer(WU, -x1_decomp)/SU_val
      dWU[, (2*k_decomp+1):(3*k_decomp)] <-
        diag((-x0_decomp + x1_decomp)/SU_val) -
        outer(WU, -x0_decomp + x1_decomp)/SU_val
      dWU[, (3*k_decomp+1):(4*k_decomp)] <-
        diag((b0_n - bs_n)/SU_val) - outer(WU, b0_n - bs_n)/SU_val
      dWU[, (4*k_decomp+1):(5*k_decomp)] <-
        diag((bs_n - b1_n)/SU_val) - outer(WU, bs_n - b1_n)/SU_val
    }

    G_Det_E <- as.numeric(est_ov["Explained"]) * dWE + outer(WE, as.numeric(G_ov[4, ]))
    G_Det_U <- as.numeric(est_ov["Unexplained"]) * dWU + outer(WU, as.numeric(G_ov[5, ]))

    est_flat <- c(est_ov, est_det_E, est_det_U); names(est_flat) <- expected_names
    G_full <- rbind(G_ov, G_Det_E, G_Det_U)
    V_full <- G_full %*% V_total %*% t(G_full)
    rownames(V_full) <- colnames(V_full) <- names(est_flat)
    se_flat <- safe_se(V_full)

  } else {
    # ----------------------------------------------------------------------------
    # ENGINE B: REPLICATION WITH FAILSAFE SHIELD & ASYNC PROGRESS BAR
    # ----------------------------------------------------------------------------
    standard_bootstrap <- vce_method == "bootstrap" && !use_svy
    standard_jackknife <- vce_method == "jackknife" && !use_svy
    # C1: restore the caller's RNG state on exit instead of leaving it altered.
    # resample::jackknife() initializes .Random.seed when it does not exist,
    # despite being deterministic, so restore that state as well.
    if (!is.null(seed) || standard_jackknife) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      if (!is.null(seed)) set.seed(seed)
    }

    replicate_center_groups <- NULL
    if (standard_bootstrap) {
      # boot::boot() performs ordinary observation-level resampling below.
      # Each index sample is converted to integer frequency weights before
      # fitting. simple = TRUE generates one sample at a time and preserves the
      # established seed-to-replicate mapping without retaining an N x R index
      # array. The conventional covariance scale remains 1/(R-1).
      replicate_scale <- 1 / (boot_reps - 1)
      replicate_rscales <- rep(1, boot_reps)
      n_replicates <- boot_reps
    } else if (standard_jackknife) {
      # resample::jackknife() supplies the N deterministic leave-one-out
      # samples. Retained rows receive the conventional JK1 factor so the
      # replicate inputs remain identical to Stata/survey JK1.
      replicate_scale <- (N_obs - 1) / N_obs
      replicate_rscales <- rep(1, N_obs)
      n_replicates <- N_obs
    } else {
      form_w <- as.formula(paste0("~", "svy_weight"))
      form_strata <- if (use_svy && !is.null(strata_var)) {
        as.formula(paste0("~", "svy_strata"))
      } else {
        NULL
      }
      form_psu <- if (use_svy && !is.null(psu_var)) {
        as.formula(paste0("~", "svy_psu"))
      } else {
        ~1
      }

      des_base <- svydesign(
        ids = form_psu,
        strata = form_strata,
        weights = form_w,
        data = df_clean,
        nest = use_svy
      )
      # C1: restore the caller's global option on exit.
      old_opts <- options(survey.lonely.psu = lonely_psu)
      on.exit(options(old_opts), add = TRUE)
      if (vce_method == "bootstrap") {
        # Rao--Wu rescaled bootstrap: sample n_h - 1 PSUs with replacement
        # within each stratum and multiply their frequency by n_h/(n_h - 1).
        des_rep <- as.svrepdesign(
          des_base,
          type = "subbootstrap",
          replicates = boot_reps
        )
      } else {
        jackknife_type <- if (is.null(form_strata)) "JK1" else "JKn"
        des_rep <- as.svrepdesign(des_base, type = jackknife_type)
      }
      survey_replicate_type <- if (vce_method == "bootstrap") {
        "Rao-Wu rescaled bootstrap"
      } else {
        jackknife_type
      }
      survey_generator_scale <- des_rep$scale
      replicate_rscales <- des_rep$rscales
      n_replicates <- length(replicate_rscales)
      # Stata's svy bootstrap with one bootstrap sample represented by each
      # bsrweight() variable uses bsn(1), hence scale 1/B. The survey package
      # supplies the Rao--Wu replicate weights; only its finite-B scale
      # convention 1/(B-1) is replaced here so results are directly comparable
      # with Stata's documented 1/B convention.
      replicate_scale <- if (vce_method == "bootstrap") {
        1 / n_replicates
      } else {
        survey_generator_scale
      }
      if (vce_method == "jackknife" && !is.null(form_strata)) {
        replicate_factors <- as.matrix(des_rep$repweights)
        candidate_groups <- vapply(
          seq_len(ncol(replicate_factors)),
          function(rep_index) {
            changed <- abs(
              replicate_factors[, rep_index] - 1
            ) > 1e-12
            changed_strata <- unique(
              df_clean$svy_strata[changed]
            )
            if (length(changed_strata) == 1L) {
              changed_strata
            } else {
              NA_character_
            }
          },
          character(1)
        )
        # Standard JKn replicates each alter exactly one stratum. Stata's
        # default pseudovalue formula is equivalent to centering those
        # replicates on their stratum-specific mean. Unusual lonely-PSU
        # adjustments that alter several strata fall back to the global mean.
        if (!anyNA(candidate_groups)) {
          replicate_center_groups <- candidate_groups
        }
      }
    }

    n_total_runs <- n_replicates + 1
    if (!quiet) {
      cat(sprintf("\n[*] Initializing %s variance estimation (%d replications)...\n", toupper(vce_method), n_total_runs - 1))
      pb <- utils::txtProgressBar(min = 0, max = n_total_runs, style = 3)
      env_pb <- new.env()
      env_pb$count <- 0
    }

    calc_yun_rep <- function(w, data) {
      if (!quiet) {
        env_pb$count <- env_pb$count + 1
        utils::setTxtProgressBar(pb, env_pb$count)
      }

      tryCatch({
        w0 <- w[idx_0]; w1 <- w[idx_1]
        if (sum(w0) == 0 || sum(w1) == 0) return(stats::setNames(rep(NA_real_, length(expected_names)), expected_names))

        m0 <- suppressWarnings(glm.fit(
          x = X_0, y = Y_0, weights = w0,
          family = family_obj, control = fit_control
        ))
        m1 <- suppressWarnings(glm.fit(
          x = X_1, y = Y_1, weights = w1,
          family = family_obj, control = fit_control
        ))
        b0 <- m0$coefficients; b1 <- m1$coefficients
        b0[is.na(b0)] <- 0; b1[is.na(b1)] <- 0

        if (ref_method %in% c("pooled", "neumark")) {
          X_pool_rep <- X_pool
          if (
            ref_method == "pooled" &&
              pooled_anchor == "centered"
          ) {
            replicate_group_mean <- sum(
              w * df_clean[[group_var]]
            ) / sum(w)
            X_pool_rep[, pooled_group_column_name] <-
              df_clean[[group_var]] - replicate_group_mean
          }
          mp <- suppressWarnings(glm.fit(
            x = X_pool_rep, y = Y_full, weights = w,
            family = family_obj, control = fit_control
          ))
          bs <- if(ref_method == "pooled") mp$coefficients[colnames(X_full)] else mp$coefficients
        } else {
          # The reference weight is a full-sample decomposition setting.
          # In particular, Stata's weight(#) remains fixed across prefix
          # replications, so Cotton's full-sample group share must not be
          # re-estimated inside each replicate.
          W <- W_ref
          bs <- W * b0 + (1 - W) * b1
        }
        bs[is.na(bs)] <- 0

        linv <- family_obj$linkinv
        P0 <- weighted.mean(linv(as.vector(X_0 %*% b0)), w0); P1 <- weighted.mean(linv(as.vector(X_1 %*% b1)), w1)
        P0s <- weighted.mean(linv(as.vector(X_0 %*% bs)), w0); P1s <- weighted.mean(linv(as.vector(X_1 %*% bs)), w1)

        oE <- P0s - P1s; oU <- (P0 - P0s) + (P1s - P1)
        x0m <- colSums(X_0 * w0)/sum(w0); x1m <- colSums(X_1 * w1)/sum(w1)

        b0d <- as.numeric(b0 %*% coefficient_transform)
        b1d <- as.numeric(b1 %*% coefficient_transform)
        bsd <- as.numeric(bs %*% coefficient_transform)
        x0d <- as.numeric(x0m %*% mean_transform)
        x1d <- as.numeric(x1m %*% mean_transform)

        sum_WE <- sum((x0d - x1d) * bsd)
        sum_WU <- sum(x0d * (b0d - bsd) + x1d * (bsd - b1d))
        WE <- if(abs(sum_WE) > 1e-15) {
          ((x0d - x1d) * bsd) / sum_WE
        } else {
          rep(0, k_decomp)
        }
        WU <- if(abs(sum_WU) > 1e-15) {
          (x0d * (b0d - bsd) + x1d * (bsd - b1d)) / sum_WU
        } else {
          rep(0, k_decomp)
        }

        dE <- as.numeric(oE) * WE; dU <- as.numeric(oU) * WU
        c("Group_0"=P0, "Group_1"=P1, "Difference"=P0-P1, "Explained"=oE, "Unexplained"=oU,
          stats::setNames(dE, paste0("Exp_", decomp_terms)),
          stats::setNames(dU, paste0("Unexp_", decomp_terms)))
      }, error = function(e) { return(stats::setNames(rep(NA_real_, length(expected_names)), expected_names)) })
    }

    if (standard_bootstrap) {
      boot_result <- boot::boot(
        data = df_clean,
        statistic = function(data, indices) {
          frequencies <- tabulate(indices, nbins = N_obs)
          calc_yun_rep(frequencies, data)
        },
        R = n_replicates,
        sim = "ordinary",
        stype = "i",
        simple = TRUE
      )
      res_matrix <- list(
        theta = boot_result$t0,
        replicates = boot_result$t
      )
    } else if (standard_jackknife) {
      row_ids <- seq_len(N_obs)
      jackknife_result <- resample::jackknife(
        data = row_ids,
        statistic = function(retained_rows) {
          if (length(retained_rows) == N_obs) {
            jackknife_weights <- w_all
          } else {
            jackknife_weights <- rep(
              N_obs / (N_obs - 1),
              N_obs
            )
            jackknife_weights[
              setdiff(row_ids, retained_rows)
            ] <- 0
          }
          calc_yun_rep(jackknife_weights, df_clean)
        },
        statisticNames = expected_names,
        trace = FALSE
      )
      res_matrix <- list(
        theta = jackknife_result$observed,
        replicates = jackknife_result$replicates
      )
    } else {
      res_matrix <- suppressWarnings(
        withReplicates(
          des_rep,
          calc_yun_rep,
          return.replicates = TRUE
        )
      )
    }
    colnames(res_matrix$replicates) <- expected_names
    if (!quiet) close(pb)

    valid_idx <- complete.cases(res_matrix$replicates)
    na_reps <- sum(!valid_idx)

    if (sum(valid_idx) < 2) {
      diagnostics$jackknife_fail <- "CRITICAL: Too many replicates failed due to empty clusters. Variance cannot be estimated."
      if (!relax) stop(diagnostics$jackknife_fail)
    }

    if (na_reps > 0 && !quiet) cat(sprintf("\n   > Info Audit: %d replication steps failed due to sample sparsity and were automatically skipped.\n", na_reps))

    est_flat <- res_matrix$theta; names(est_flat) <- expected_names

    if (sum(valid_idx) >= 2) {
      reps <- res_matrix$replicates[valid_idx, , drop=FALSE]
      valid_center_groups <- if (
        !is.null(replicate_center_groups)
      ) {
        replicate_center_groups[valid_idx]
      } else {
        NULL
      }
      if (is.null(valid_center_groups)) {
        replicate_center <- colMeans(reps)
        diffs <- sweep(reps, 2, replicate_center, "-")
        replicate_center_method <- "replicate mean"
      } else {
        group_centers <- lapply(
          unique(valid_center_groups),
          function(group) {
            colMeans(
              reps[valid_center_groups == group, , drop = FALSE]
            )
          }
        )
        names(group_centers) <- unique(valid_center_groups)
        replicate_center <- do.call(
          rbind,
          lapply(valid_center_groups, function(group) {
            group_centers[[group]]
          })
        )
        colnames(replicate_center) <- expected_names
        diffs <- reps - replicate_center
        replicate_center_method <-
          "stratum-specific replicate mean"
      }
      rscales <- if (
        length(replicate_rscales) == length(valid_idx)
      ) {
        replicate_rscales[valid_idx]
      } else {
        replicate_rscales
      }
      # Dropping failed replicates removes terms from the variance sum and biases it
      # downward; rescale by (total replicates / valid replicates) to compensate.
      rep_adjust <- length(valid_idx) / sum(valid_idx)
      if (rep_adjust > 1) {
        diagnostics$rep_adjust <- sprintf("Variance Note: %d failed replicate(s) dropped; variance rescaled by %.3f to compensate.", na_reps, rep_adjust)
      }
      effective_scale <- replicate_scale * rep_adjust
      V_full <- effective_scale * crossprod(diffs * sqrt(rscales))
      se_flat <- safe_se(V_full)
      rownames(V_full) <- colnames(V_full) <- names(est_flat)
      replication_details <- list(
        center_method = replicate_center_method,
        center = if (is.matrix(replicate_center)) {
          replicate_center
        } else {
          stats::setNames(replicate_center, expected_names)
        },
        replicates = reps,
        requested_replicates = length(valid_idx),
        valid_replicates = sum(valid_idx),
        failed_replicates = na_reps,
        scale = effective_scale,
        generator_scale = if (
          exists("survey_generator_scale", inherits = FALSE)
        ) {
          survey_generator_scale
        } else {
          replicate_scale
        },
        rscales = rscales,
        ordinary_bootstrap = standard_bootstrap,
        ordinary_jackknife = standard_jackknife,
        engine = if (standard_bootstrap) {
          "boot::boot"
        } else if (standard_jackknife) {
          "resample::jackknife"
        } else {
          "survey::withReplicates"
        },
        simulation = if (standard_bootstrap) {
          "ordinary"
        } else if (standard_jackknife) {
          "delete-one"
        } else {
          NULL
        },
        statistic_type = if (standard_bootstrap) {
          "indices"
        } else if (standard_jackknife) {
          "retained row indices"
        } else {
          NULL
        },
        fit_weight_type = if (standard_bootstrap) {
          "frequency"
        } else if (standard_jackknife) {
          "JK1 replicate factor"
        } else {
          NULL
        },
        survey_replicate_type = if (
          exists("survey_replicate_type", inherits = FALSE)
        ) {
          survey_replicate_type
        } else {
          NULL
        }
      )
    } else {
      V_full <- matrix(NA_real_, length(est_flat), length(est_flat))
      se_flat <- rep(NA_real_, length(est_flat))
    }
  }

  # Match Jann's e(b) layout by omitting the structurally zero explained
  # contribution of the intercept.
  explained_intercept <- "Exp_(Intercept)"
  if (explained_intercept %in% names(est_flat)) {
    keep_output <- names(est_flat) != explained_intercept
    est_flat <- est_flat[keep_output]
    se_flat <- se_flat[keep_output]
    V_full <- V_full[keep_output, keep_output, drop = FALSE]
    if (!is.null(replication_details)) {
      replication_details$center <- if (
        is.matrix(replication_details$center)
      ) {
        replication_details$center[, keep_output, drop = FALSE]
      } else {
        replication_details$center[keep_output]
      }
      replication_details$replicates <-
        replication_details$replicates[, keep_output, drop = FALSE]
    }
  }

  # --- PLAUSIBILITY SHIELD 4: SE Explosion (relative to the outcome scale) ---
  # B3: compare the maximum SE to the scale of the group predictions rather than
  # to a fixed absolute threshold (which is meaningless across outcome scales).
  max_se <- suppressWarnings(max(se_flat, na.rm = TRUE))
  scale_ref <- max(abs(est_flat[c("Group_0", "Group_1")]), na.rm = TRUE)
  scale_ref <- if (is.finite(scale_ref) && scale_ref > 0) scale_ref else 1
  if (!is.infinite(max_se) && max_se > 10 * scale_ref) {
    diagnostics$se_explosion <- sprintf("Instability Warning: Standard Error is %.1f times the outcome scale (Max SE = %.1f). Likely sparse categories.", max_se / scale_ref, max_se)
    if (!relax && !quiet) warning(diagnostics$se_explosion)
  }

  # --- PLAUSIBILITY SHIELD 5: Yun weight degeneracy ---
  # When the linear-scale weight denominator of a component is ~0 the detailed
  # weights are zeroed, but in logit/probit the (probability-scale) component can
  # still be nonzero: detailed rows then do not sum to the component. Flag it.
  chk_tol <- 1e-8 * max(1, abs(est_flat[["Difference"]]))
  bad_E <- abs(sum(est_flat[grep("^Exp_", names(est_flat))]) - est_flat[["Explained"]]) > chk_tol
  bad_U <- abs(sum(est_flat[grep("^Unexp_", names(est_flat))]) - est_flat[["Unexplained"]]) > chk_tol
  if (isTRUE(bad_E) || isTRUE(bad_U)) {
    diagnostics$yun_weights <- "Yun Weights Warning: the linear-scale weight denominator of a component is ~0, so its detailed contributions were zeroed and do NOT sum to the (nonlinear) component. Interpret detailed rows for that component with caution."
    if (!relax && !quiet) warning(diagnostics$yun_weights, call. = FALSE)
  }

  raw_estimates <- est_flat
  raw_vcov <- V_full
  raw_standard_errors <- se_flat

  # ==============================================================================
  # 5. DYNAMIC VARIABLE GROUPING
  # ==============================================================================
  names_res <- names(est_flat)
  grp_exp <- numeric(); grp_unexp <- numeric(); se_grp_E <- numeric(); se_grp_U <- numeric()

  if (length(groupings) > 0) {
    # Match priority: exact model term > exact column name > prefix (literal, not
    # regex, so names with special characters are safe). Prefix matching can
    # over-capture (e.g. "edu" also grabs "education..."), so flag it when a
    # prefix spans columns from more than one model term.
    x_cols <- decomp_terms
    column_to_term <- unlist(lapply(names(term_column_map), function(term) {
      stats::setNames(
        rep(term, length(term_column_map[[term]])),
        term_column_map[[term]]
      )
    }))
    ambiguous_prefixes <- character(0)
    match_group_cols <- function(v) {
      if (v %in% names(term_column_map)) return(term_column_map[[v]])
      if (v %in% x_cols) return(v)
      hits <- x_cols[startsWith(x_cols, v)]
      if (length(hits) > 0) {
        terms_hit <- unique(unname(column_to_term[hits]))
        if (length(terms_hit) > 1) ambiguous_prefixes <<- c(ambiguous_prefixes, v)
      }
      hits
    }

    for (i in seq_along(groupings)) {
      grp_name <- names(groupings)[i]

      g_cols <- unique(unlist(lapply(groupings[[i]], match_group_cols)))
      idx_E <- match(paste0("Exp_", g_cols), names_res); idx_E <- idx_E[!is.na(idx_E)]
      idx_U <- match(paste0("Unexp_", g_cols), names_res); idx_U <- idx_U[!is.na(idx_U)]

      if (length(idx_E) > 0) {
        grp_exp[grp_name] <- sum(est_flat[idx_E]); se_grp_E[grp_name] <- safe_se_val(sum(V_full[idx_E, idx_E]))
      } else { grp_exp[grp_name] <- NA; se_grp_E[grp_name] <- NA }

      if (length(idx_U) > 0) {
        grp_unexp[grp_name] <- sum(est_flat[idx_U]); se_grp_U[grp_name] <- safe_se_val(sum(V_full[idx_U, idx_U]))
      } else { grp_unexp[grp_name] <- NA; se_grp_U[grp_name] <- NA }
    }
    if (length(ambiguous_prefixes) > 0) {
      diagnostics$groupings <- sprintf("Grouping Warning: prefix(es) %s matched columns from more than one model term. Verify the intended variables; use exact term names to disambiguate.", paste(sQuote(unique(ambiguous_prefixes)), collapse = ", "))
      if (!quiet) warning(diagnostics$groupings, call. = FALSE)
    }
    names(grp_exp) <- paste0("GrpExp_", names(grp_exp)); names(grp_unexp) <- paste0("GrpUnexp_", names(grp_unexp))
    names(se_grp_E) <- names(grp_exp); names(se_grp_U) <- names(grp_unexp)

    est_flat <- c(est_flat, grp_exp, grp_unexp); se_flat <- c(se_flat, se_grp_E, se_grp_U)
  }

  # ==============================================================================
  # 6. FORMATTING EXPORT TIBBLES (WITH DYNAMIC DICTIONARY MAPPING)
  # ==============================================================================
  q_val <- 1 - ((1 - level) / 2)
  make_tibble <- function(pattern, rm_str) {
    idx <- if(pattern == "Overall") 1:5 else grep(pattern, names(est_flat))
    if(length(idx) == 0) return(NULL)

    est <- est_flat[idx]
    se <- se_flat[idx]
    stat <- est / se

    # Evaluate the upper tail directly. Subtracting pt() from one rounds very
    # small but nonzero p-values to zero (unlike Stata's coefficient table).
    p_val <- 2 * pt(abs(stat), df = design_df, lower.tail = FALSE)
    ci_low <- est - qt(q_val, df = design_df) * se
    ci_high <- est + qt(q_val, df = design_df) * se

    clean_terms <- gsub(rm_str, "", names(est))

    display_terms <- unname(sapply(clean_terms, function(x) {
      if (x %in% names(display_dict)) return(display_dict[[x]])
      return(x)
    }))

    if (pattern == "Overall") {
      return(
        tibble::tibble(
          Term = display_terms,
          Estimate = unname(est),
          Std_Error = unname(se),
          Statistic = unname(stat),
          P_Value = unname(p_val),
          Conf_Low = unname(ci_low),
          Conf_High = unname(ci_high)
        )
      )
    }

    total_gap <- as.numeric(est_flat["Difference"])
    pct_total <- if (abs(total_gap) > 1e-12) (est / total_gap) * 100 else rep(NA_real_, length(est))

    is_explained_component <- grepl("Exp_", pattern)

    if (is_explained_component) {
      explained_gap <- as.numeric(est_flat["Explained"])
      pct_explained <- if (abs(explained_gap) > 1e-12) {
        (est / explained_gap) * 100
      } else {
        rep(NA_real_, length(est))
      }

      return(
        tibble::tibble(
          Term = display_terms,
          Estimate = unname(est),
          `% Contribution Explained` = unname(pct_explained),
          `% Contribution Total` = unname(pct_total),
          Std_Error = unname(se),
          Statistic = unname(stat),
          P_Value = unname(p_val),
          Conf_Low = unname(ci_low),
          Conf_High = unname(ci_high)
        )
      )
    }

    unexplained_gap <- as.numeric(est_flat["Unexplained"])
    pct_unexplained <- if (abs(unexplained_gap) > 1e-12) {
      (est / unexplained_gap) * 100
    } else {
      rep(NA_real_, length(est))
    }

    tibble::tibble(
      Term = display_terms,
      Estimate = unname(est),
      `% Contribution Unexplained` = unname(pct_unexplained),
      `% Contribution Total` = unname(pct_total),
      Std_Error = unname(se),
      Statistic = unname(stat),
      P_Value = unname(p_val),
      Conf_Low = unname(ci_low),
      Conf_High = unname(ci_high)
    )
  }

  tb_overall <- make_tibble("Overall", "")
  tb_det_E <- make_tibble("^Exp_", "^Exp_"); tb_det_U <- make_tibble("^Unexp_", "^Unexp_")
  tb_grp_E <- if(length(groupings)>0) make_tibble("^GrpExp_", "^GrpExp_") else NULL
  tb_grp_U <- if(length(groupings)>0) make_tibble("^GrpUnexp_", "^GrpUnexp_") else NULL

  # ==============================================================================
  # 7. CONSOLE OUTPUT
  # ==============================================================================
  if (!quiet) {
    vce_label <- if (vce_method == "linearized") {
      if (use_svy) "linearized (survey design)" else "linearized (Huber-White)"
    } else {
      vce_method
    }
    cat("\n")
    cat(rep("-", 75), "\n", sep="")
    cat("OAXACA-BLINDER YUN DECOMPOSITION\n")
    cat(rep("-", 75), "\n", sep="")
    cat(sprintf("Number of strata = %-15s Number of obs   = %d\n", ifelse(use_svy, n_strata, "N/A"), N_obs))
    cat(sprintf("Number of PSUs   = %-15s Population size = %.4f\n", ifelse(use_svy, n_psu, "N/A"), pop_size))
    cat(sprintf("                                   Design df       = %s\n", ifelse(is.infinite(design_df), "Inf", design_df)))
    cat(sprintf("Method           = %-15s Model           = %s\n", ref_method, model_type))
    if (ref_method == "pooled") {
      cat(sprintf(
        "Pooled anchor    = %-15s Reference point = %s\n",
        pooled_anchor,
        pooled_anchor_description
      ))
    }
    cat(sprintf("VCE Engine       = %s\n", vce_label))
    cat(sprintf(
      "Normalized factors = %s\n",
      if (length(normalized_factors) > 0) {
        paste(normalized_factors, collapse = ", ")
      } else {
        "none"
      }
    ))
    cat(sprintf("Group 0 (favored)      : %-16s N of obs 0 = %d\n", paste0(group_var, " = ", favored_label), N_0))
    cat(sprintf("Group 1 (disadvantaged): %-16s N of obs 1 = %d\n", paste0(group_var, " = ", disadvantaged_label), N_1))
    cat("Gap orientation  : Difference = E(Y|favored) - E(Y|disadvantaged)\n")
    cat(rep("-", 75), "\n", sep="")

    if (length(diagnostics) > 0) {
      cat("\n[DIAGNOSTICS & WARNINGS]\n")
      for (warn in diagnostics) cat(" *", warn, "\n")
    }

    cat("\n[MODEL FIT METRICS]\n")
    print(as.data.frame(fit_metrics), row.names = FALSE)
    cat("\n[DECOMPOSITION SUMMARY]\n")
    print(as.data.frame(tb_overall), row.names = FALSE)

    cat("\nNote: Detailed categorical terms are labeled as reference contrasts or\n")
    cat("      normalized deviations; continuous predictors indicate '(Continuous)'.\n")
  }

  resultados <- list(
    summary_stats = list(
      N = N_obs,
      Pop = pop_size,
      Strata = n_strata,
      PSUs = n_psu,
      DF = design_df,
      N0 = N_0,
      N1 = N_1,
      favored = favored_label,
      disadvantaged = disadvantaged_label,
      group_selection = group_selection_mode,
      filtered_groups = excluded_group_levels,
      filtered_observations = n_group_filtered
    ),
    model_metrics = fit_metrics,
    diagnostics = diagnostics,
    results_overall = tb_overall,
    results_detailed_explained = tb_det_E,
    results_detailed_unexplained = tb_det_U,
    results_grouped_explained = tb_grp_E,
    results_grouped_unexplained = tb_grp_U,
    models = models_list,
    raw = list(
      estimates = raw_estimates,
      standard_errors = raw_standard_errors,
      vcov = raw_vcov,
      group_coefficients = list(
        group_0 = b_0_decomp,
        group_1 = b_1_decomp
      ),
      reference_coefficients = b_star_decomp,
      covariate_means = list(
        group_0 = x0_decomp,
        group_1 = x1_decomp
      ),
      analytic_sample = df_clean$.oby_source_row,
      model_terms = decomp_terms,
      fit_model_terms = fit_model_terms,
      base_estimates = base_estimates,
      base_vcov = if (exists("V_total", inherits = FALSE)) V_total else NULL,
      normalization = list(
        requested = normalize,
        factors = normalized_factors,
        active = length(normalized_factors) > 0,
        method = "equal-category deviation contrasts",
        blocks = normalization_blocks
      ),
      pooled_anchor = pooled_anchor_info,
      replication = replication_details,
      group_mapping = list(
        favored = favored_label,
        disadvantaged = disadvantaged_label,
        selected = c(favored_label, disadvantaged_label),
        available = available_group_levels,
        excluded = excluded_group_levels,
        mode = group_selection_mode
      ),
      vce_method = vce_method,
      survey_mode = use_svy
    )
  )
  return(invisible(resultados))
}
