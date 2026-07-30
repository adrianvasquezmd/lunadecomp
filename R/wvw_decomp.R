#' Wagstaff-Type Decomposition of Socioeconomic Health Inequality
#'
#' @description
#' Decomposes a rank-dependent health inequality index into the contributions of
#' observed determinants using the regression-based framework proposed by
#' Wagstaff, van Doorslaer, and Watanabe.
#'
#' The method expresses the concentration index of a health variable as the sum
#' of determinant-specific contributions plus a residual component. Each
#' contribution is constructed from the determinant's marginal effect on the
#' outcome, its mean, and its socioeconomic concentration index. The function
#' supports linear models and nonlinear binary-response models through average
#' marginal effects.
#'
#' @details
#' Let `y` denote the health variable and `R` the fractional socioeconomic rank.
#' The standard concentration index can be written as:
#'
#' \deqn{
#'   C_y = \frac{2 \operatorname{cov}(y, R)}{\mu_y}.
#' }
#'
#' For a linear model,
#'
#' \deqn{
#'   y_i = \alpha + \sum_k \beta_k x_{ki} + \varepsilon_i,
#' }
#'
#' the decomposition can be written as:
#'
#' \deqn{
#'   C_y = \sum_k \left(\frac{\beta_k \bar{x}_k}{\mu_y}\right) C_k
#'         + \frac{GC_\varepsilon}{\mu_y},
#' }
#'
#' where `C_k` is the concentration index of determinant `x_k`,
#' `GC_epsilon` is the generalized concentration index of the residual, and
#' `beta_k * mean(x_k) / mean(y)` is the elasticity of `y` with respect to
#' `x_k` evaluated at the sample mean.
#'
#' For `model_type = "ols"`, the marginal effects are the fitted regression
#' coefficients. For `model_type = "logit"` or `model_type = "probit"`, the
#' function estimates average marginal effects. Continuous predictors use the
#' weighted average derivative of the inverse-link function. A separate binary
#' predictor is handled as the weighted average discrete change from zero to
#' one. A factor level is handled as the discrete change from the omitted level,
#' setting the complete factor block to the target category.
#'
#' Mutually exclusive categories should therefore be supplied as a factor.
#' Supplying their one-hot columns as separate predictors estimates a different
#' counterfactual: each dummy is toggled independently while the other dummy
#' columns retain their observed values. Although that parameterization has the
#' same likelihood and fitted values, its nonlinear marginal effects and
#' decomposition need not equal the factor-based result.
#'
#' Logit and probit models use a tightened Fisher-scoring criterion
#' (`epsilon = 1e-12`, `maxit = 100`) to avoid preventable differences from
#' Stata caused by the looser [stats::glm.fit()] default. Samples exhibiting
#' perfect classification with extreme coefficients are rejected because the
#' finite maximum-likelihood coefficients do not exist.
#'
#' Unlike the exact OLS identity, the nonlinear marginal-effect decomposition is
#' generally an approximation. The returned overall table therefore includes
#' `Non-linear Approximation Error`, defined so that total inequality equals the
#' explained component plus the concentration of the prediction residual plus
#' this approximation error.
#'
#' @section Categorical predictors and groupings:
#' Each entry in `indep_vars` must be the name of a column in `data`;
#' transformations and interaction expressions are not accepted. Character
#' predictors are converted to factors and factors are expanded by
#' [stats::model.matrix()]. The first factor level is the omitted reference.
#' Changing that reference changes individual dummy coefficients and
#' contributions, but not fitted values or the contribution of the complete
#' factor block.
#'
#' Custom `groupings` must have unique, non-empty names and may refer to exact
#' input terms, exact model-matrix columns, or unambiguous column prefixes. A
#' model column cannot occur in more than one group. Unknown, ambiguous, and
#' overlapping specifications are rejected. Predictor columns omitted from a
#' custom grouping are retained automatically under their original terms, so
#' grouped contributions always form an exact partition of the explained part.
#'
#' @section Analytic sample and sampling weights:
#' A single complete-case sample is formed from every variable required by the
#' selected model, rank, and survey design. Missing required values are excluded
#' before fitting and their original row positions are recorded in `raw`.
#' Sampling weights must be numeric or losslessly coercible to numeric, finite,
#' and non-negative. Zero-weight rows are excluded from the analytic sample, as
#' in Stata's weighted estimation, and are recorded separately. Negative,
#' non-finite, nonnumeric, and all-zero weights are rejected.
#'
#' @section Correction methods:
#' The argument `correction` determines the scale of the inequality index being
#' decomposed. Internally, the function starts from the weighted generalized
#' concentration covariance term `GC = 2 * cov(y, R)` and applies the following
#' scaling:
#'
#' \describe{
#'   \item{`"standard"`}{Standard concentration index. Uses
#'   `GC / mean(y)`. This is the classic relative concentration index and is most
#'   appropriate for positive, ratio-scale, effectively unbounded outcomes.}
#'
#'   \item{`"generalized"`}{Generalized concentration index. Uses `GC` without
#'   division by the outcome mean. This represents an absolute rank-dependent
#'   inequality measure.}
#'
#'   \item{`"erreygers"`}{Erreygers-corrected concentration index. Uses
#'   `4 * GC / (dep_max - dep_min)`. This correction is intended for bounded
#'   cardinal outcomes and reflects an absolute inequality value judgment. It
#'   requires meaningful lower and upper bounds.}
#'
#'   \item{`"wagstaff"`}{Wagstaff-normalized concentration index. Uses
#'   `(dep_max - dep_min) * GC / ((dep_max - mean(y)) * (mean(y) - dep_min))`.
#'   This correction rescales the concentration index by its feasible bounds and
#'   is often discussed for bounded or binary variables. It reflects a different
#'   value judgment from the Erreygers correction and can be sensitive when the
#'   outcome mean is close to its bounds.}
#' }
#'
#' The choice between Erreygers and Wagstaff corrections is not purely technical.
#' For bounded and binary outcomes, the standard concentration index has
#' mean-dependent bounds. The Erreygers and Wagstaff corrections address this
#' issue using different normative assumptions. Users should report the chosen
#' correction and justify whether the outcome is interpreted as an attainment
#' or as a shortfall. For OLS outcomes, `dep_min` and `dep_max` must be supplied
#' as theoretical bounds; the function never substitutes the observed sample
#' minimum and maximum. Logit/probit outcomes use their theoretical bounds
#' `[0, 1]`. Values outside the declared bounds are rejected.
#'
#' @section Outcome orientation:
#' If the outcome is bounded, interpretation depends on whether larger values
#' indicate a favorable attainment, such as coverage or health, or an adverse
#' shortfall, such as disease or malnutrition. Set `is_shortfall = TRUE` when
#' larger values represent a worse outcome. This does not recode the data; it
#' records the orientation in diagnostics so that interpretation of the sign is
#' transparent.
#'
#' @section Reported decomposition columns:
#' The detailed table reports the following core quantities:
#'
#' \describe{
#'   \item{`Marginal Effects`}{Estimated marginal effect of the determinant on
#'   the health outcome. For OLS this corresponds to the regression coefficient;
#'   for logit/probit it is an average marginal effect.}
#'   \item{`Elasticity`}{`Marginal Effects * Mean_X / mean(y)`. This is the
#'   elasticity used in the Wagstaff decomposition for the standard concentration
#'   index.}
#'   \item{`Concentration Index`}{Socioeconomic concentration index of the
#'   determinant.}
#'   \item{`Estimate`}{Absolute contribution of the determinant to the selected
#'   inequality index (named `Estimate` for consistency with [oby_decomp()] and
#'   [f_decomp()]; it is the determinant's contribution).}
#'   \item{`% Contribution Explained`}{Contribution as a percentage of the sum
#'   of all explained contributions.}
#'   \item{`% Contribution Total`}{Contribution as a percentage of the total
#'   selected inequality index.}
#' }
#'
#' @section Survey designs and variance estimation:
#' If `use_svy = TRUE`, the function uses the supplied final sampling weight,
#' strata, and primary sampling unit (PSU). Supplying any of `weight_var`,
#' `strata_var`, or `psu_var` activates survey mode automatically, so design
#' variables cannot be ignored silently. Reused PSU labels are nested within
#' strata.
#'
#' Standard errors can be estimated using analytic linearization, jackknife, or
#' bootstrap. In `"linearized"` mode the function uses a joint Huber--White /
#' Taylor sandwich for the complete vector of contributions, total inequality,
#' residual inequality, and nonlinear approximation error. The estimating
#' equations jointly propagate the fitted coefficients, average marginal
#' effects, determinant concentration covariances, the outcome mean, and the
#' selected correction factor. When `ses_var` is used, the influence of the
#' empirically estimated fractional rank is included, with weighted ties handled
#' as a block. A rank supplied through `precalc_rank_var` is instead treated as
#' fixed, because its generating data and design are not available to the
#' function.
#'
#' For auditability, `raw$linearization$conindex_robust_benchmark` also reports
#' the robust convenient-regression variance used by Stata's `conindex, robust`.
#' That historical benchmark conditions on the computed rank and normalization
#' scalars, so it need not equal the function's fuller rank-aware variance.
#' `raw$linearization$conindex_survey_benchmark` analogously stores the
#' convenient-regression estimate and Taylor standard error under the declared
#' survey design.
#'
#' With only a sampling weight, `"linearized"` retains the observation-level
#' Huber--White convention and asymptotic-normal inference (`DF = Inf`). With
#' strata or PSUs it applies Taylor ultimate-cluster covariance and uses
#' `number of PSUs - number of strata` design degrees of freedom, matching
#' Stata's `svy` convention.
#'
#' Outside survey mode, [resample::jackknife()] generates one leave-one-out
#' sample per analytic observation and the covariance uses scale
#' `(N - 1) / N`. [boot::boot()] generates ordinary samples of `N`
#' observations with replacement and the covariance uses scale `1 / (B - 1)`.
#' Both methods center on the replicate mean, matching Stata's default non-MSE
#' `jackknife:`/`bootstrap:` prefixes. Every replicate recomputes the fitted
#' model, fractional ranks, marginal effects, correction factor, contributions,
#' total, and residual. Replicate estimates, weights, centers, scales, engine
#' metadata, and the complete covariance matrix are retained in
#' `raw$replication`.
#' If an ordinary replicate fails, the non-MSE covariance uses the number of
#' valid replicates: `(R - 1) / R` for jackknife and `1 / (R - 1)` for
#' bootstrap, matching Stata's handling of missing replicate statistics.
#'
#' In survey mode, `"jackknife"` uses JK1 when there is no stratification and
#' JKn when strata are declared. JKn deletes one PSU at a time, inflates the
#' remaining PSUs in that stratum, and centers each replicate on the mean of
#' the replicates belonging to the same stratum. `"bootstrap"` uses the
#' Rao--Wu rescaled bootstrap: within each stratum it samples `n_h - 1` PSUs
#' with replacement, applies the `n_h / (n_h - 1)` rescaling, centers on the
#' replicate mean, and uses final scale `1 / B`. Both survey-replicate methods
#' use `number of PSUs - number of strata` degrees of freedom. A sampling
#' weight is never interpreted as an ordinary frequency-bootstrap weight.
#'
#' Survey jackknife requires at least two PSUs in every stratum. Survey
#' bootstrap has the same requirement by default, but
#' `bootstrap_singleton = "certainty"` uses
#' [svrep::as_bootstrap_design()] to hold a singleton PSU fixed with replicate
#' factor 1 and assign it zero first-stage variance. Regular strata retain
#' Rao--Wu \eqn{n_h-1} resampling. This should be selected only when treating
#' the singleton as a certainty PSU is substantively defensible. Survey
#' replication stops if any model replicate fails; `relax = TRUE` permits only
#' an explicitly exploratory failed-replicate adjustment. Taylor
#' linearization supports lonely-PSU rules through `lonely_psu`: `"adjust"`
#' corresponds to Stata `singleunit(centered)`, `"certainty"` and `"remove"`
#' to `singleunit(certainty)`, and `"average"` to `singleunit(scaled)`.
#'
#' The built-in survey contract is an ultimate-cluster, with-replacement
#' approximation with one final weight, optional strata, and one PSU stage.
#' Finite-population corrections, lower-stage identifiers, first-stage PPS
#' resampling, recalibration, and externally supplied replicate weights are not
#' represented. Such designs require a future external-design interface rather
#' than being approximated silently.
#'
#' @param data A data frame containing the outcome, determinants, socioeconomic
#' ranking variable, and optional survey design variables.
#' @param dep_var Character string. Name of the health outcome variable.
#' @param indep_vars Character vector. Bare column names of the explanatory
#' variables to include in the decomposition model. Character variables are
#' converted to factors; factor variables are expanded using `model.matrix()`.
#' Formula transformations and interactions are not supported.
#' @param ses_var Optional character string. Name of the socioeconomic variable
#' used to rank observations from most disadvantaged to most advantaged. Required
#' unless `precalc_rank_var` is supplied.
#' @param groupings Optional named list. Each element contains exact variable
#' names, exact model-matrix column names, or unambiguous model-matrix prefixes
#' to aggregate detailed dummy-level contributions into conceptual groups.
#' Unknown, ambiguous, overlapping, or duplicate specifications are rejected.
#' Unlisted determinant columns are retained automatically by original term. If
#' empty, all terms are auto-grouped by their original model terms.
#' @param correction Character string. Inequality index correction to decompose.
#' Options are `"standard"`, `"generalized"`, `"erreygers"`, and `"wagstaff"`.
#' @param model_type Character string. Model family. Options are `"ols"`,
#' `"logit"`, and `"probit"`.
#' @param precalc_rank_var Optional character string. Name of a pre-calculated
#' fractional socioeconomic rank with finite values in the closed interval
#' `[0, 1]`. If supplied, its values are used directly instead of computing
#' ranks from `ses_var`.
#' @param dep_min Optional numeric value. Theoretical lower bound of the outcome.
#' Required for Erreygers/Wagstaff corrections with OLS outcomes; inferred as
#' zero for logit/probit outcomes.
#' @param dep_max Optional numeric value. Theoretical upper bound of the outcome.
#' Required for Erreygers/Wagstaff corrections with OLS outcomes; inferred as
#' one for logit/probit outcomes.
#' @param is_shortfall Logical. If `TRUE`, indicates that larger values of the
#' outcome represent worse health or a deficit/shortfall.
#' @param use_svy Logical. If `TRUE`, applies survey design settings. Supplying
#' any survey design variable also activates this mode automatically.
#' @param weight_var Optional character string. Name of the sampling weight
#' variable. Values must be finite and non-negative. Rows with zero weight are
#' excluded from the analytic sample and recorded in `raw`.
#' @param strata_var Optional character string. Name of the stratification
#' variable.
#' @param psu_var Optional character string. Name of the primary sampling unit
#' variable.
#' @param vce_method Character string. Variance estimation method. Options are
#' `"linearized"`, `"jackknife"`, and `"bootstrap"`.
#' @param boot_reps Integer. Number of bootstrap replications when
#' `vce_method = "bootstrap"`. Ordinary observation-level samples are generated
#' by [boot::boot()]; survey bootstrap continues to use replicate weights from
#' the declared survey design.
#' @param lonely_psu Character string. Taylor handling of strata with one PSU.
#' Options are `"fail"`, `"remove"`, `"certainty"`, `"adjust"`, and
#' `"average"`. Survey jackknife requires at least two PSUs in every stratum;
#' survey-bootstrap singleton handling is controlled separately by
#' `bootstrap_singleton`.
#' @param bootstrap_singleton Character string. Handling of a stratum with one
#' PSU in survey bootstrap. Options are `"fail"` (the default) and
#' `"certainty"`. `"fail"` requires at least two PSUs per stratum and stops
#' rather than making an assumption. `"certainty"` fixes the singleton at
#' replicate factor 1 through [svrep::as_bootstrap_design()], assigning zero
#' first-stage variance to that stratum.
#' @param level Numeric. Confidence level for confidence intervals. Default is
#' `0.95`.
#' @param seed Optional integer. Random seed for reproducibility of replicate
#' resampling. Defaults to `NULL` (non-deterministic); set an explicit value to
#' reproduce results. The prior state of the global random number generator is
#' restored on exit.
#' @param relax Logical. If `TRUE`, selected validation problems are stored as
#' diagnostics or warnings instead of stopping execution.
#' @param quiet Logical. If `TRUE`, suppresses console output.
#'
#' @return
#' An invisible named list. The inferential columns (`Term`, `Estimate`,
#' `Std_Error`, `Statistic`, `P_Value`, `Conf_Low`, `Conf_High`) are harmonized
#' with [oby_decomp()] and [f_decomp()]; method-specific columns are preserved.
#' Numeric values are returned at full machine precision. Rounding is applied
#' only to console display.
#' \describe{
#'   \item{summary_stats}{Sample size, population size, number of strata, number
#'   of PSUs, and design degrees of freedom.}
#'   \item{model_metrics}{Model fit information.}
#'   \item{diagnostics}{Methodological and numerical diagnostics.}
#'   \item{results_overall}{Total inequality index, explained component, and
#'   residual/unexplained component (plus a non-linear approximation error row for
#'   logit/probit), each with standard error, test statistic, p-value, and
#'   confidence interval, and a `Percentage` column giving the share of the total
#'   index.}
#'   \item{results_grouped}{Grouped determinant contributions with the shared
#'   inferential columns and a `Dummies` count of aggregated model-matrix terms.}
#'   \item{results_detailed}{Detailed determinant-level table: `Marginal
#'   Effects`, `Mean_X`, `Elasticity`, `Concentration Index`, the contribution
#'   (`Estimate`), percent contributions, and inferential statistics.}
#'   \item{raw}{Unrounded audit information: analytic source rows, exclusion
#'   causes (including missing required values and zero sampling weights), weights,
#'   design identifiers, fractional ranks, model matrix and coefficients,
#'   determinant means and generalized concentration covariances, correction
#'   scale, unrounded estimate vectors, covariance matrices and, for replicate
#'   methods, replicate estimates, factors, centers, scales, engine metadata,
#'   and validity indicators. For `linearized`, it also contains the complete joint
#'   influence-function linearization. The `complete` flags identify covariance
#'   matrices that contain all cross-covariances.}
#' }
#'
#' @examples
#' data(lunadecomp_example)
#'
#' fit <- wvw_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "health_score",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   ses_var = "ses",
#'   model_type = "ols",
#'   correction = "standard",
#'   vce_method = "linearized",
#'   quiet = TRUE
#' )
#'
#' fit$results_overall
#' fit$results_detailed
#'
#' \donttest{
#' fit_erreygers <- wvw_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "health_score",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   ses_var = "ses",
#'   model_type = "ols",
#'   correction = "erreygers",
#'   dep_min = 0,
#'   dep_max = 1,
#'   vce_method = "linearized",
#'   quiet = TRUE
#' )
#'
#' fit_erreygers$results_overall
#' }
#'
#' @author
#' Adrian Vasquez-Mejia, MD, MSc \cr
#' Oscar J. Mujica, MD, MPH, PHE, FACE \cr
#' Antonio Sanhueza, MPH, MSc, PhD
#'
#' @references
#' Wagstaff, A., van Doorslaer, E., & Watanabe, N. (2003). On decomposing the
#' causes of health sector inequalities with an application to malnutrition
#' inequalities in Vietnam. \emph{Journal of Econometrics}, 112(1), 207--223.
#' \doi{10.1016/S0304-4076(02)00161-6}
#'
#' Wagstaff, A. (2005). The bounds of the concentration index when the variable
#' of interest is binary, with an application to immunization inequality.
#' \emph{Health Economics}, 14(4), 429--432. \doi{10.1002/hec.953}
#'
#' Erreygers, G. (2009). Correcting the concentration index.
#' \emph{Journal of Health Economics}, 28(2), 504--515.
#' \doi{10.1016/j.jhealeco.2008.02.003}
#'
#' Kjellsson, G., & Gerdtham, U.-G. (2013). On correcting the concentration
#' index for binary variables. \emph{Journal of Health Economics}, 32(4),
#' 659--670. \doi{10.1016/j.jhealeco.2012.10.012}
#'
#' @importFrom dplyr %>% select all_of mutate mutate_if row_number arrange desc case_when
#' @importFrom tidyr drop_na
#' @importFrom tibble tibble add_row as_tibble
#' @importFrom rlang .data
#' @importFrom survey svydesign svyglm as.svrepdesign withReplicates
#' @importFrom stats as.formula complete.cases cov.wt dlogis dnorm gaussian glm.fit median model.matrix na.omit plogis pnorm pt qt quasibinomial terms var vcov update
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
wvw_decomp <- function(
    data, dep_var, indep_vars, ses_var = NULL, groupings = list(),
    correction = "standard", model_type = "ols",
    precalc_rank_var = NULL,
    dep_min = NULL, dep_max = NULL, is_shortfall = FALSE,
    use_svy = FALSE, weight_var = NULL, strata_var = NULL, psu_var = NULL,
    vce_method = "linearized", boot_reps = 500, lonely_psu = "adjust",
    level = 0.95, seed = NULL, relax = FALSE, quiet = FALSE,
    bootstrap_singleton = "fail"
) {

  # C1: restore the caller's RNG state on exit instead of leaving it altered.
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
  if (!is.logical(use_svy) || length(use_svy) != 1L ||
      is.na(use_svy)) {
    stop("Error: use_svy must be exactly TRUE or FALSE.")
  }
  if (
    !use_svy &&
      any(!vapply(
        list(weight_var, strata_var, psu_var),
        is.null,
        logical(1)
      ))
  ) {
    use_svy <- TRUE
  }

  # ==============================================================================
  # 1. STRUCTURAL VALIDATION & DATA CONFIGURATION
  # ==============================================================================
  valid_vce <- c("linearized", "jackknife", "bootstrap")
  valid_corr <- c("erreygers", "wagstaff", "standard", "generalized")
  valid_models <- c("ols", "logit", "probit")
  valid_lonely_psu <- c(
    "fail", "remove", "certainty", "adjust", "average"
  )
  valid_bootstrap_singleton <- c("fail", "certainty")

  if (!(vce_method %in% valid_vce)) stop("Error: vce_method must be 'linearized', 'jackknife', or 'bootstrap'.")
  if (!(correction %in% valid_corr)) stop("Error: correction must be 'erreygers', 'wagstaff', 'standard', or 'generalized'.")
  if (!(model_type %in% valid_models)) stop("Error: model_type must be 'ols', 'logit', or 'probit'.")
  if (
    !is.character(lonely_psu) ||
      length(lonely_psu) != 1L ||
      is.na(lonely_psu) ||
      !(lonely_psu %in% valid_lonely_psu)
  ) {
    stop(
      "Error: lonely_psu must be 'fail', 'remove', 'certainty', ",
      "'adjust', or 'average'."
    )
  }
  if (
    !is.character(bootstrap_singleton) ||
      length(bootstrap_singleton) != 1L ||
      is.na(bootstrap_singleton) ||
      !(bootstrap_singleton %in% valid_bootstrap_singleton)
  ) {
    stop(
      "Error: bootstrap_singleton must be 'fail' or 'certainty'."
    )
  }
  if (
    identical(vce_method, "bootstrap") &&
      (
        !is.numeric(boot_reps) ||
          length(boot_reps) != 1L ||
          is.na(boot_reps) ||
          !is.finite(boot_reps) ||
          boot_reps < 2 ||
          boot_reps != as.integer(boot_reps)
      )
  ) {
    stop("Error: boot_reps must be one integer greater than or equal to 2.")
  }
  if (!is.logical(is_shortfall) || length(is_shortfall) != 1L ||
      is.na(is_shortfall)) {
    stop("Error: is_shortfall must be exactly TRUE or FALSE.")
  }
  if (!is.data.frame(data)) {
    stop("Error: data must be a data frame.")
  }
  if (!is.character(dep_var) || length(dep_var) != 1L ||
      is.na(dep_var) || dep_var == "") {
    stop("Error: dep_var must name exactly one variable.")
  }
  if (!is.character(indep_vars) || length(indep_vars) == 0L ||
      anyNA(indep_vars) || any(indep_vars == "")) {
    stop("Error: indep_vars must be a non-empty character vector of variable names.")
  }
  if (anyDuplicated(indep_vars)) {
    stop("Error: indep_vars contains duplicated variable names.")
  }
  if (!is.list(groupings)) {
    stop("Error: groupings must be a named list.")
  }

  if (is.null(ses_var) && is.null(precalc_rank_var)) {
    stop("Error: You must provide either 'ses_var' or 'precalc_rank_var' to establish the socioeconomic rank.")
  }

  old_lonely_psu <- getOption("survey.lonely.psu")
  options(survey.lonely.psu = lonely_psu)
  on.exit(options(survey.lonely.psu = old_lonely_psu), add = TRUE)

  req_vars <- unique(c(dep_var, ses_var, precalc_rank_var, indep_vars, weight_var, strata_var, psu_var))
  req_vars <- req_vars[!is.na(req_vars) & req_vars != ""]
  absent_variables <- setdiff(req_vars, names(data))
  if (length(absent_variables) > 0L) {
    stop(
      "Error: variables not found in data: ",
      paste(absent_variables, collapse = ", "),
      "."
    )
  }

  # Preserve the original row positions so every rank, model-matrix row, and
  # replicate result can be traced back to the caller's data after complete-case
  # filtering and SES sorting.
  source_row_name <- ".__wvw_source_row__"
  while (source_row_name %in% names(data)) {
    source_row_name <- paste0(".", source_row_name)
  }
  data_audit <- data
  data_audit[[source_row_name]] <- seq_len(nrow(data_audit))
  input_n <- nrow(data_audit)

  df_clean <- data_audit %>%
    dplyr::select(dplyr::all_of(c(source_row_name, req_vars))) %>%
    tidyr::drop_na() %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    droplevels()

  complete_case_rows <- df_clean[[source_row_name]]
  missing_excluded_rows <- setdiff(seq_len(input_n), complete_case_rows)
  zero_weight_excluded_rows <- integer(0)

  # Preserve the full binary representation of values that are already
  # numeric. Routing them through as.character() silently rounded externally
  # supplied fractional ranks before converting them back to double.
  numeric_preserving_coercion <- function(x) {
    if (is.numeric(x)) return(as.numeric(x))
    suppressWarnings(as.numeric(as.character(x)))
  }
  df_clean[[dep_var]] <- numeric_preserving_coercion(df_clean[[dep_var]])
  if (!is.null(ses_var)) {
    df_clean[[ses_var]] <- numeric_preserving_coercion(df_clean[[ses_var]])
  }
  if (!is.null(precalc_rank_var)) {
    df_clean[[precalc_rank_var]] <-
      numeric_preserving_coercion(df_clean[[precalc_rank_var]])
  }
  if (!is.null(weight_var)) {
    df_clean[[weight_var]] <-
      numeric_preserving_coercion(df_clean[[weight_var]])
  }

  numeric_check_vars <- unique(c(
    dep_var,
    ses_var,
    precalc_rank_var,
    weight_var
  ))
  numeric_check_vars <- numeric_check_vars[!is.na(numeric_check_vars) & numeric_check_vars != ""]
  bad_numeric <- numeric_check_vars[vapply(numeric_check_vars, function(v) any(is.na(df_clean[[v]])), logical(1))]

  if (length(bad_numeric) > 0) {
    stop("CRITICAL DATA ERROR: Numeric conversion introduced NA values in: ", paste(bad_numeric, collapse = ", "), ". Check variable coding.")
  }

  nonfinite_numeric <- numeric_check_vars[
    vapply(
      numeric_check_vars,
      function(v) any(!is.finite(df_clean[[v]])),
      logical(1)
    )
  ]
  if (length(nonfinite_numeric) > 0L) {
    stop(
      "CRITICAL DATA ERROR: Non-finite values detected in: ",
      paste(nonfinite_numeric, collapse = ", "),
      "."
    )
  }

  if (!is.null(precalc_rank_var)) {
    precalculated_rank <- df_clean[[precalc_rank_var]]
    if (any(!is.finite(precalculated_rank))) {
      stop(
        "CRITICAL DATA ERROR: precalc_rank_var must contain only finite values."
      )
    }
    if (any(precalculated_rank < 0 | precalculated_rank > 1)) {
      stop(
        "CRITICAL DATA ERROR: precalc_rank_var must be a fractional rank in the closed interval [0, 1]."
      )
    }
  }

  if (!is.null(weight_var)) {
    supplied_weight <- df_clean[[weight_var]]
    if (any(supplied_weight < 0)) {
      stop("CRITICAL DATA ERROR: Sampling weights cannot be negative.")
    }
    zero_weight <- supplied_weight == 0
    if (any(zero_weight)) {
      zero_weight_excluded_rows <-
        df_clean[[source_row_name]][zero_weight]
      df_clean <- droplevels(df_clean[!zero_weight, , drop = FALSE])
    }
    if (nrow(df_clean) == 0L) {
      stop(
        "CRITICAL DATA ERROR: No positive sampling weight remains in the analytic sample."
      )
    }
  }

  N_obs <- nrow(df_clean)
  if (!quiet && N_obs < nrow(data)) {
    cat(sprintf(
      "\n[!] Data Prep: %d row(s) with missing required values or zero weights were excluded. N = %d\n",
      nrow(data) - N_obs,
      N_obs
    ))
  }
  if (N_obs == 0) {
    stop(
      "CRITICAL: Zero observations remaining after required-value and weight filtering."
    )
  }

  if (model_type %in% c("logit", "probit")) {
    unique_y <- unique(df_clean[[dep_var]])
    if (!all(unique_y %in% c(0, 1))) stop("CRITICAL: For logit/probit, dep_var must be strictly binary (0 and 1).")
    if (length(unique_y) < 2) stop("CRITICAL: dep_var has no variation (all values are ", unique_y[1], "). A degenerate outcome cannot be decomposed.")
    dep_min <- 0; dep_max <- 1
  }

  sort_var <- if (!is.null(precalc_rank_var)) precalc_rank_var else ses_var
  if (!is.null(sort_var)) df_clean <- df_clean[order(df_clean[[sort_var]]), ]

  df_clean <- df_clean %>%
    dplyr::mutate(
      svy_weight = if (use_svy && !is.null(weight_var)) as.numeric(.data[[weight_var]]) else 1,
      svy_strata = if (use_svy && !is.null(strata_var)) as.character(.data[[strata_var]]) else "1",
      svy_psu    = if (use_svy && !is.null(psu_var)) as.character(.data[[psu_var]]) else as.character(dplyr::row_number())
    )

  pop_size <- sum(df_clean$svy_weight)
  n_strata <- if(use_svy) length(unique(df_clean$svy_strata)) else NA
  # W2: count PSUs within strata (nested), not globally.
  n_psu    <- if(use_svy) nrow(unique(df_clean[, c("svy_strata", "svy_psu")])) else NA
  psu_per_stratum <- if (use_svy) {
    table(
      unique(
        df_clean[, c("svy_strata", "svy_psu")]
      )$svy_strata
    )
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
  diagnostics <- list()
  if (length(zero_weight_excluded_rows) > 0L) {
    diagnostics$zero_weights <- sprintf(
      "Data Note: %d row(s) with zero sampling weight were excluded from the analytic sample.",
      length(zero_weight_excluded_rows)
    )
  }
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

  # W3: SEs from a variance without masking non-positive terms via abs().
  safe_se <- function(v) {
    neg <- which(v < -1e-8)
    if (length(neg) > 0) {
      diagnostics$neg_var <<- sprintf(
        "Numerical Warning: %d non-positive variance term(s) detected (non-PSD covariance); their SE set to NA.",
        length(neg)
      )
      if (!relax && !quiet) warning(diagnostics$neg_var, call. = FALSE)
    }
    v[v < 0] <- NA_real_
    sqrt(v)
  }

  # === 2. METHODOLOGICAL ALERTS ===
  if (vce_method == "linearized") {
    diagnostics$vce <- paste(
      "Methodology Note: 'linearized' uses a joint Huber-White/Taylor sandwich",
      "for coefficients, marginal effects, concentration covariances, correction",
      "scale, total, residual, and approximation error. Empirical ranks are",
      "rank-aware; a user-supplied precalculated rank is treated as fixed."
    )
    if (!quiet) message(diagnostics$vce)
  }

  if (correction %in% c("erreygers", "wagstaff")) {
    if (is.null(dep_min) || is.null(dep_max)) {
      stop(
        "CRITICAL: Erreygers and Wagstaff corrections require explicit ",
        "theoretical dep_min and dep_max bounds for OLS outcomes. ",
        "Observed sample minima and maxima are not valid substitutes."
      )
    }
    valid_bound <- function(x) {
      is.numeric(x) && length(x) == 1L && is.finite(x)
    }
    if (!valid_bound(dep_min) || !valid_bound(dep_max) ||
        dep_min >= dep_max) {
      stop(
        "CRITICAL: dep_var bounds must be finite numeric scalars with ",
        "dep_min < dep_max."
      )
    }
    observed_min <- min(df_clean[[dep_var]], na.rm = TRUE)
    observed_max <- max(df_clean[[dep_var]], na.rm = TRUE)
    bound_tolerance <- sqrt(.Machine$double.eps) *
      max(1, abs(dep_min), abs(dep_max))
    if (observed_min < dep_min - bound_tolerance ||
        observed_max > dep_max + bound_tolerance) {
      stop(
        sprintf(
          paste0(
            "CRITICAL: dep_var values [%.17g, %.17g] fall outside ",
            "the declared theoretical bounds [%.17g, %.17g]."
          ),
          observed_min,
          observed_max,
          dep_min,
          dep_max
        )
      )
    }

    orient_msg <- ifelse(
      is_shortfall,
      paste(
        "Interpretation Alert: The outcome is coded as a SHORTFALL",
        "(higher values = worse outcome). Interpret the sign as inequality",
        "in the adverse outcome. Consider recoding to attainment if you",
        "want positive values to represent better health/coverage."
      ),
      paste(
        "Interpretation Alert: The outcome is coded as an ATTAINMENT",
        "(higher values = better outcome). Interpret the sign as inequality",
        "in the favorable outcome."
      )
    )
    diagnostics$orientation <- orient_msg
    if (!quiet) warning(diagnostics$orientation, call. = FALSE)
  } else {
    dep_min <- NA; dep_max <- NA
  }

  if (use_svy) {
    max_w <- max(df_clean$svy_weight, na.rm = TRUE)
    med_w <- median(df_clean$svy_weight, na.rm = TRUE)
    if (med_w > 0 && (max_w > 50 * med_w)) {
      diagnostics$extreme_weights <- sprintf("Leverage Warning: Extreme sampling weights detected (Max is %.0f times the median). Resampling VCE strongly advised.", max_w/med_w)
      if (!quiet) warning(diagnostics$extreme_weights, call. = FALSE)
    }
    if (
      vce_method == "jackknife" &&
        any(psu_per_stratum < 2)
    ) {
      stop(
        "Error: survey replication requires at least two PSUs in every ",
        "stratum. Use Taylor linearization with an explicit lonely_psu ",
        "rule, collapse or redesign singleton strata, or supply an ",
        "appropriate external replicate-weight design.",
        call. = FALSE
      )
    }
    if (
      vce_method == "bootstrap" &&
        bootstrap_has_singletons &&
        bootstrap_singleton == "fail"
    ) {
      stop(
        "Error: Rao-Wu survey bootstrap requires at least two PSUs in ",
        "every stratum. Set bootstrap_singleton = \"certainty\" only when ",
        "assigning zero first-stage variance to singleton strata is ",
        "defensible; otherwise collapse or redesign singleton strata, or ",
        "supply an appropriate external replicate-weight design.",
        call. = FALSE
      )
    }
    if (
      vce_method == "bootstrap" &&
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
    if (vce_method %in% c("jackknife", "bootstrap")) {
      diagnostics$survey_replication_scope <- paste0(
        "Survey Replication Scope: ultimate-cluster variance using one ",
        "final weight, optional strata, and one PSU stage; no FPC, ",
        "lower-stage identifiers, calibration replicate weights, or ",
        "first-stage PPS resampling is represented."
      )
      if (!quiet) {
        message(diagnostics$survey_replication_scope)
      }
    }
  }

  # === 3. WAGSTAFF MATH CALCULATORS ===
  calc_cov_r <- function(y_var, r_var, w) {
    w_norm <- w / sum(w)
    return(2 * (sum(w_norm * y_var * r_var) - sum(w_norm * y_var) * sum(w_norm * r_var)))
  }

  calc_cov_r_vec <- function(X_mat, r_var, w) {
    w_norm <- w / sum(w)
    return(2 * (colSums(w_norm * X_mat * r_var) - colSums(w_norm * X_mat) * sum(w_norm * r_var)))
  }

  get_scale_m <- function(mu_y, corr_type, min_h, max_h) {
    if (corr_type == "erreygers") return(4 / (max_h - min_h))
    if (corr_type == "wagstaff") {
      denom <- (max_h - mu_y) * (mu_y - min_h)
      if (denom <= 1e-12) stop("CRITICAL: Wagstaff correction is unstable because the outcome mean is too close to its theoretical bounds.")
      return((max_h - min_h) / denom)
    }
    if (corr_type == "standard") {
      if (mu_y <= 1e-12) {
        stop(
          "CRITICAL: Mean of dep_var is zero, negative, or nearly zero. ",
          "The 'standard' index requires a strictly positive mean."
        )
      }
      return(1 / mu_y)
    }
    return(1)
  }

  get_ames <- function(betas, mm, wt, m_type, assign_idx) {
    if (m_type == "ols") return(betas)
    xb <- as.vector(mm %*% betas)
    ames <- numeric(length(betas)); names(ames) <- names(betas)
    w_norm <- wt / sum(wt)

    for (k in seq_along(betas)) {
      if (k == 1) { ames[k] <- 0; next }
      term_id <- assign_idx[k]
      cols_in_term <- which(assign_idx == term_id)
      is_bin <- all(stats::na.omit(mm[, k]) %in% c(0, 1))

      if (is_bin) {
        mm1 <- mm; mm0 <- mm
        if (length(cols_in_term) > 1) { mm1[, cols_in_term] <- 0; mm0[, cols_in_term] <- 0 }
        mm1[, k] <- 1; mm0[, k] <- 0
        if (m_type == "logit") {
          p1 <- stats::plogis(as.vector(mm1 %*% betas)); p0 <- stats::plogis(as.vector(mm0 %*% betas))
        } else {
          p1 <- stats::pnorm(as.vector(mm1 %*% betas)); p0 <- stats::pnorm(as.vector(mm0 %*% betas))
        }
        ames[k] <- sum(w_norm * (p1 - p0))
      } else {
        # Numerically stable logistic density (exp(xb)/(1+exp(xb))^2 overflows to NaN for xb > ~709)
        pdf_val <- if (m_type == "logit") stats::dlogis(xb) else stats::dnorm(xb)
        ames[k] <- sum(w_norm * betas[k] * pdf_val)
      }
    }
    return(ames)
  }

  get_scale_derivative <- function(mu_y, corr_type, min_h, max_h) {
    if (corr_type == "standard") return(-1 / mu_y^2)
    if (corr_type == "wagstaff") {
      range_h <- max_h - min_h
      denom <- (max_h - mu_y) * (mu_y - min_h)
      denom_derivative <- max_h + min_h - 2 * mu_y
      return(-range_h * denom_derivative / denom^2)
    }
    0
  }

  # Influence function of GC(v) = 2 cov(v, R). If R is generated from the
  # empirical SES distribution, GC is a weighted order-2 functional and both
  # appearances of an observation (health value and rank distribution) matter.
  # For an externally supplied rank, only the ordinary covariance influence is
  # identifiable and the supplied rank is conditioned on.
  get_gc_influence <- function(
      values,
      ranks,
      norm_weights,
      empirical_rank_values = NULL
  ) {
    values <- as.matrix(values)
    n <- nrow(values)
    k_values <- ncol(values)
    means <- colSums(norm_weights * values)
    rank_mean <- sum(norm_weights * ranks)
    gc <- 2 * (
      colSums(norm_weights * values * ranks) -
        means * rank_mean
    )

    if (is.null(empirical_rank_values)) {
      centered_values <- sweep(values, 2, means, "-")
      influence <- 2 * centered_values * as.numeric(ranks - rank_mean)
      influence <- sweep(influence, 2, gc, "-")
    } else {
      # df_clean is sorted by SES, so equal values form contiguous blocks.
      new_group <- c(
        TRUE,
        empirical_rank_values[-1L] !=
          empirical_rank_values[-length(empirical_rank_values)]
      )
      group_id <- cumsum(new_group)
      weighted_values <- values * norm_weights
      group_totals <- rowsum(
        weighted_values,
        group = group_id,
        reorder = FALSE
      )
      cumulative_totals <- apply(group_totals, 2, cumsum)
      if (is.null(dim(cumulative_totals))) {
        cumulative_totals <- matrix(
          cumulative_totals,
          ncol = k_values
        )
      }
      all_totals <- matrix(
        colSums(group_totals),
        nrow = nrow(group_totals),
        ncol = k_values,
        byrow = TRUE
      )
      higher_group_totals <- all_totals - cumulative_totals
      second_argument_term <- (
        higher_group_totals + 0.5 * group_totals
      )[group_id, , drop = FALSE]
      a_value <- colSums(norm_weights * values * ranks)
      influence_a <- values * ranks + second_argument_term
      influence_a <- sweep(influence_a, 2, 2 * a_value, "-")
      influence <- 2 * influence_a -
        sweep(values, 2, means, "-")
    }

    # The exact weighted mean is zero analytically. Centering only removes
    # floating-point accumulation error and stabilizes covariance identities.
    influence <- sweep(
      influence,
      2,
      colSums(norm_weights * influence),
      "-"
    )
    colnames(influence) <- colnames(values)
    list(
      estimate = gc,
      influence = influence,
      rank_treated_as = if (is.null(empirical_rank_values)) {
        "fixed precalculated rank"
      } else {
        "estimated weighted empirical fractional rank"
      }
    )
  }

  get_ame_linearization <- function(
      betas,
      mm,
      norm_weights,
      m_type,
      assign_idx,
      beta_influence
  ) {
    k_beta <- length(betas)
    n <- nrow(mm)
    gradient <- matrix(
      0,
      nrow = k_beta,
      ncol = k_beta,
      dimnames = list(names(betas), names(betas))
    )
    direct <- matrix(
      0,
      nrow = n,
      ncol = k_beta,
      dimnames = list(NULL, names(betas))
    )

    if (m_type == "ols") {
      diag(gradient) <- 1
      ame <- betas
      influence <- beta_influence
      colnames(influence) <- names(betas)
      return(list(
        estimate = ame,
        gradient_beta = gradient,
        direct_influence = direct,
        influence = influence
      ))
    }

    inverse_link <- if (m_type == "logit") {
      stats::plogis
    } else {
      stats::pnorm
    }
    link_density <- if (m_type == "logit") {
      stats::dlogis
    } else {
      stats::dnorm
    }
    eta <- as.vector(mm %*% betas)
    ame <- numeric(k_beta)
    names(ame) <- names(betas)

    for (column in seq_len(k_beta)) {
      if (column == 1L) next
      term_id <- assign_idx[column]
      term_columns <- which(assign_idx == term_id)
      is_binary <- all(
        stats::na.omit(mm[, column]) %in% c(0, 1)
      )

      if (is_binary) {
        mm_one <- mm
        mm_zero <- mm
        if (length(term_columns) > 1L) {
          mm_one[, term_columns] <- 0
          mm_zero[, term_columns] <- 0
        }
        mm_one[, column] <- 1
        mm_zero[, column] <- 0
        eta_one <- as.vector(mm_one %*% betas)
        eta_zero <- as.vector(mm_zero %*% betas)
        individual_effect <- inverse_link(eta_one) -
          inverse_link(eta_zero)
        gradient[column, ] <- colSums(
          norm_weights * (
            mm_one * link_density(eta_one) -
              mm_zero * link_density(eta_zero)
          )
        )
      } else {
        density <- link_density(eta)
        density_derivative <- if (m_type == "logit") {
          probability <- stats::plogis(eta)
          density * (1 - 2 * probability)
        } else {
          -eta * density
        }
        individual_effect <- betas[column] * density
        gradient[column, ] <- colSums(
          norm_weights * (
            betas[column] * mm * density_derivative
          )
        )
        gradient[column, column] <-
          gradient[column, column] +
          sum(norm_weights * density)
      }

      ame[column] <- sum(norm_weights * individual_effect)
      direct[, column] <- individual_effect - ame[column]
    }

    influence <- direct + beta_influence %*% t(gradient)
    influence <- sweep(
      influence,
      2,
      colSums(norm_weights * influence),
      "-"
    )
    colnames(influence) <- names(betas)
    list(
      estimate = ame,
      gradient_beta = gradient,
      direct_influence = direct,
      influence = influence
    )
  }

  # === 4. DESIGN MATRIX, COLLINEARITY & DYNAMIC DICTIONARY ===
  one_level_factors <- indep_vars[
    vapply(
      indep_vars,
      function(variable) {
        is.factor(df_clean[[variable]]) &&
          nlevels(df_clean[[variable]]) < 2L
      },
      logical(1)
    )
  ]
  if (length(one_level_factors) > 0L) {
    stop(
      "CRITICAL MODEL VALIDATION: Factor predictor(s) have fewer than two observed levels: ",
      paste(one_level_factors, collapse = ", "),
      "."
    )
  }

  form_base <- stats::as.formula(paste("~", paste(indep_vars, collapse = " + ")))
  if (any(attr(stats::terms(form_base), "order") > 1)) {
    stop("CRITICAL: Interaction terms detected. Analytical AMEs with interactions break Wagstaff's aditivity and are disabled.")
  }

  X_full <- stats::model.matrix(form_base, data = df_clean)
  c_names <- colnames(X_full)

  if (any(!is.finite(X_full))) {
    bad_columns <- colnames(X_full)[
      apply(X_full, 2, function(column) any(!is.finite(column)))
    ]
    stop(
      "CRITICAL DATA ERROR: Non-finite values detected in model-matrix columns: ",
      paste(bad_columns, collapse = ", "),
      "."
    )
  }

  if (ncol(X_full) >= N_obs) {
    stop(sprintf("CRITICAL MODEL VALIDATION: Predictors (%d) exceeds or equals observations (%d). Unidentified model.", ncol(X_full), N_obs))
  }
  if (use_svy && !is.null(psu_var) && ncol(X_full) >= n_psu) {
    diagnostics$df_psu <- sprintf("Model Validation Warning: Predictors (%d) exceeds or equals PSUs (%d). VCE matrix is not full rank.", ncol(X_full), n_psu)
    if (!relax && !quiet) warning(diagnostics$df_psu, call. = FALSE)
  }

  weighted_matrix <- X_full * sqrt(df_clean$svy_weight)
  weighted_qr <- qr(weighted_matrix, tol = 1e-10)
  if (weighted_qr$rank < ncol(X_full)) {
    aliased_columns <- colnames(X_full)[
      weighted_qr$pivot[
        seq.int(weighted_qr$rank + 1L, ncol(X_full))
      ]
    ]
    stop(
      "CRITICAL MODEL VALIDATION: Exact collinearity detected in model-matrix column(s): ",
      paste(aliased_columns, collapse = ", "),
      "."
    )
  }

  var_check <- apply(X_full[, -1, drop=FALSE], 2, function(x) stats::var(x, na.rm = TRUE))
  if (any(var_check == 0)) {
    stop("CRITICAL: Zero variance detected in predictors: ", paste(names(var_check)[var_check == 0], collapse=", "), ". Please remove them.")
  }

  assign_attr <- attr(X_full, "assign"); term_labels <- attr(stats::terms(form_base), "term.labels")
  display_dict <- list(); auto_groups <- list()

  for (i in seq_along(term_labels)) {
    term <- term_labels[i]; cols <- which(assign_attr == i)
    if(length(cols) > 0) {
      auto_groups[[term]] <- cols
      if (is.factor(df_clean[[term]]) || is.character(df_clean[[term]])) {
        ref_level <- levels(as.factor(df_clean[[term]]))[1]
        for (cn in c_names[cols]) display_dict[[cn]] <- sprintf("%s (%s vs. Ref: %s)", term, substring(cn, nchar(term) + 1), ref_level)
      } else {
        for (cn in c_names[cols]) display_dict[[cn]] <- sprintf("%s (Continuous)", term)
      }
    }
  }
  display_dict[["(Intercept)"]] <- "Baseline Reference (Intercept)"

  final_groups <- list()
  if (length(groupings) == 0) {
    for (term in names(auto_groups)) {
      cols <- auto_groups[[term]]
      if (length(cols) == 1) {
        col_name <- c_names[cols]
        g_name <- if(!is.null(display_dict[[col_name]])) display_dict[[col_name]] else term
      } else {
        ref_level <- levels(as.factor(df_clean[[term]]))[1]
        g_name <- sprintf("%s (Factor block: %d categories vs. Ref: %s)", term, length(cols), ref_level)
      }
      final_groups[[g_name]] <- c_names[cols]
    }
  } else {
    if (is.null(names(groupings)) || any(names(groupings) == "") ||
        anyDuplicated(names(groupings))) {
      stop("CRITICAL GROUPING ERROR: groupings must have unique, non-empty names.")
    }
    ambiguous_prefixes <- character(0)
    unmatched_tokens <- character(0)
    for (g in names(groupings)) {
      if (!is.character(groupings[[g]]) ||
          length(groupings[[g]]) == 0L ||
          anyNA(groupings[[g]]) ||
          any(groupings[[g]] == "")) {
        stop(
          "CRITICAL GROUPING ERROR: group '",
          g,
          "' must contain one or more variable, term, or column names."
        )
      }
      matched_cols <- unlist(lapply(groupings[[g]], function(v) {
        # Match priority: exact term name > exact column name > prefix. Prefix matching
        # can over-capture (e.g. "edu" also grabs "education..."), so flag it when the
        # prefix spans columns from more than one model term.
        if (v %in% names(auto_groups)) return(c_names[auto_groups[[v]]])
        if (v %in% c_names) return(v)
        hits <- c_names[startsWith(c_names, v)]
        hits <- setdiff(hits, "(Intercept)")
        if (length(hits) == 0L) {
          unmatched_tokens <<- c(unmatched_tokens, v)
          return(character(0))
        }
        if (length(hits) > 0) {
          terms_hit <- unique(assign_attr[match(hits, c_names)])
          if (length(terms_hit) > 1) ambiguous_prefixes <<- c(ambiguous_prefixes, v)
        }
        return(hits)
      }))
      matched_cols <- unique(intersect(matched_cols, c_names[c_names != "(Intercept)"]))
      if (length(matched_cols) > 0) final_groups[[g]] <- matched_cols
    }
    if (length(unmatched_tokens) > 0L) {
      stop(
        "CRITICAL GROUPING ERROR: unknown grouping token(s): ",
        paste(unique(unmatched_tokens), collapse = ", "),
        "."
      )
    }
    if (length(ambiguous_prefixes) > 0) {
      stop(
        "CRITICAL GROUPING ERROR: ambiguous prefix(es) matched columns from more than one model term: ",
        paste(unique(ambiguous_prefixes), collapse = ", "),
        ". Use exact term or model-column names."
      )
    }

    assigned_columns <- unlist(final_groups, use.names = FALSE)
    duplicated_columns <- unique(
      assigned_columns[duplicated(assigned_columns)]
    )
    if (length(duplicated_columns) > 0L) {
      stop(
        "CRITICAL GROUPING ERROR: model columns assigned to more than one group: ",
        paste(duplicated_columns, collapse = ", "),
        "."
      )
    }

    determinant_columns <- setdiff(c_names, "(Intercept)")
    uncovered_columns <- setdiff(
      determinant_columns,
      assigned_columns
    )
    if (length(uncovered_columns) > 0L) {
      diagnostics$groupings_uncovered <- sprintf(
        paste0(
          "Grouping Note: %d uncovered model column(s) were retained ",
          "automatically in their original model terms."
        ),
        length(uncovered_columns)
      )
      for (term in names(auto_groups)) {
        remaining <- intersect(
          c_names[auto_groups[[term]]],
          uncovered_columns
        )
        if (length(remaining) == 0L) next
        if (length(remaining) == 1L) {
          candidate <- if (!is.null(display_dict[[remaining]])) {
            display_dict[[remaining]]
          } else {
            term
          }
        } else {
          candidate <- paste0(term, " (Ungrouped remainder)")
        }
        candidate <- make.unique(
          c(names(final_groups), candidate),
          sep = " #"
        )[length(final_groups) + 1L]
        final_groups[[candidate]] <- remaining
      }
    }

    final_assignment <- unlist(final_groups, use.names = FALSE)
    if (!setequal(final_assignment, determinant_columns) ||
        anyDuplicated(final_assignment)) {
      stop(
        "INTERNAL GROUPING ERROR: grouped columns do not form an exact partition of the determinant matrix."
      )
    }
  }

  if (ncol(X_full) > 2) {
    cor_X <- suppressWarnings(stats::cov.wt(X_full[, -1, drop = FALSE], wt = df_clean$svy_weight, cor = TRUE)$cor)
    if (!is.null(cor_X) && !any(is.na(cor_X))) {
      inv_cor <- tryCatch(solve(cor_X), error = function(e) NULL)
      if (is.null(inv_cor)) {
        diagnostics$collinearity <- "CRITICAL MODEL VALIDATION: Perfect collinearity (Singular Matrix) detected."
        if (!relax) stop(diagnostics$collinearity) else if (!quiet) warning(diagnostics$collinearity, call. = FALSE)
      } else if (any(diag(inv_cor) > 15)) {
        diagnostics$vif <- "Model Validation Warning: Severe Multicollinearity detected (VIF > 15). SEs may be unstable."
        if (!relax && !quiet) warning(diagnostics$vif, call. = FALSE)
      }
    }
  }

  # === 5. MOTOR PRINCIPAL ===
  form_mod <- stats::as.formula(paste(dep_var, "~", paste(indep_vars, collapse = " + ")))
  form_w <- stats::as.formula("~svy_weight")
  form_strata <- if (use_svy && !is.null(strata_var)) stats::as.formula("~svy_strata") else NULL
  form_psu <- if (use_svy && !is.null(psu_var)) stats::as.formula("~svy_psu") else ~1
  des_base <- survey::svydesign(ids = form_psu, strata = form_strata, weights = form_w, data = df_clean, nest = use_svy)

  fam_link <- if (model_type == "ols") stats::gaussian(link = "identity")
              else if (model_type == "logit") stats::quasibinomial(link = "logit")
              else stats::quasibinomial(link = "probit")
  # Match the tightened nonlinear convergence policy already validated in
  # oby_decomp() and f_decomp(). The glm.fit() default deviance tolerance can
  # stop probit several Fisher-scoring steps before Stata's solution.
  fit_control <- stats::glm.control(
    epsilon = if (model_type %in% c("logit", "probit")) 1e-12 else 1e-8,
    maxit = if (model_type %in% c("logit", "probit")) 100 else 25
  )

  get_r_i <- function(data_df, w, precalc_val, ses_name) {
    if (!is.null(precalc_val)) return(precalc_val)
    sigma <- w / sum(w)
    vals <- rle(as.numeric(data_df[[ses_name]]))
    agg_sigma <- as.numeric(tapply(sigma, rep(seq_along(vals$lengths), vals$lengths), sum))
    cum_sigma <- cumsum(agg_sigma)
    return(rep(c(0, cum_sigma[-length(cum_sigma)]) + 0.5 * agg_sigma, vals$lengths))
  }

  W_base <- df_clean$svy_weight; w_norm <- W_base / sum(W_base)
  Y_base <- df_clean[[dep_var]]; mu_y_base <- sum(w_norm * Y_base)

  R_base <- get_r_i(df_clean, W_base, if(!is.null(precalc_rank_var)) df_clean[[precalc_rank_var]] else NULL, if(!is.null(ses_var)) ses_var else names(df_clean)[1])
  M_base <- get_scale_m(mu_y_base, correction, dep_min, dep_max)
  Total_Index <- M_base * calc_cov_r(Y_base, R_base, W_base)

  # Stata conindex's directly comparable benchmark. Its robust convenient
  # regression conditions on the estimated rank, rank variance, mean, and
  # normalization factor. It is retained as an audit target, while the public
  # linearized VCE below uses the fuller rank-aware joint influence function.
  des_base <- update(des_base, .wvw_rank = R_base)
  s2R <- sum(w_norm * (R_base - sum(w_norm * R_base))^2)
  conindex_robust_benchmark <- function(vals, scale_M) {
    if (N_obs <= 2L) {
      return(list(
        estimate = NA_real_,
        std_error = NA_real_,
        vcov = matrix(NA_real_, 2, 2)
      ))
    }
    square_root_weight <- sqrt(W_base)
    benchmark_matrix <- cbind(
      `(Intercept)` = square_root_weight,
      rank = R_base * square_root_weight
    )
    benchmark_response <-
      2 * s2R * scale_M * vals * square_root_weight
    if (qr(benchmark_matrix)$rank < ncol(benchmark_matrix)) {
      return(list(
        estimate = 0,
        std_error = NA_real_,
        coefficients = stats::setNames(
          c(NA_real_, 0),
          colnames(benchmark_matrix)
        ),
        vcov = matrix(
          NA_real_,
          nrow = 2,
          ncol = 2,
          dimnames = list(
            colnames(benchmark_matrix),
            colnames(benchmark_matrix)
          )
        ),
        finite_sample_correction =
          "undefined because the fractional rank has no variation",
        conditioned_on =
          "fractional ranks, rank variance, outcome mean, and correction scale"
      ))
    }
    benchmark_bread <- solve(crossprod(benchmark_matrix))
    benchmark_coefficients <- as.vector(
      benchmark_bread %*%
        crossprod(benchmark_matrix, benchmark_response)
    )
    benchmark_residual <- benchmark_response -
      as.vector(benchmark_matrix %*% benchmark_coefficients)
    benchmark_vcov <- (
      N_obs / (N_obs - ncol(benchmark_matrix))
    ) * benchmark_bread %*%
      crossprod(benchmark_matrix * benchmark_residual) %*%
      benchmark_bread
    dimnames(benchmark_vcov) <- list(
      colnames(benchmark_matrix),
      colnames(benchmark_matrix)
    )
    list(
      estimate = unname(benchmark_coefficients[2]),
      std_error = sqrt(benchmark_vcov[2, 2]),
      coefficients = stats::setNames(
        benchmark_coefficients,
        colnames(benchmark_matrix)
      ),
      vcov = benchmark_vcov,
      finite_sample_correction =
        "HC1 N/(N-2), matching conindex, robust",
      conditioned_on =
        "fractional ranks, rank variance, outcome mean, and correction scale"
    )
  }
  conindex_total_benchmark <-
    conindex_robust_benchmark(Y_base, M_base)
  conindex_survey_benchmark <- function(vals, scale_M) {
    if (!use_svy) return(NULL)
    benchmark_terms <- c(
      "(Intercept)",
      ".wvw_benchmark_rank"
    )
    if (
      !is.finite(s2R) ||
        s2R <= .Machine$double.eps ||
        length(unique(R_base)) < 2L
    ) {
      return(list(
        estimate = 0,
        std_error = NA_real_,
        coefficients = stats::setNames(
          c(NA_real_, 0),
          benchmark_terms
        ),
        vcov = matrix(
          NA_real_,
          nrow = 2L,
          ncol = 2L,
          dimnames = list(benchmark_terms, benchmark_terms)
        ),
        degrees_of_freedom = survey::degf(des_base),
        conditioned_on = paste(
          "undefined survey convenient regression because the",
          "fractional rank has no variation"
        )
      ))
    }
    benchmark_design <- update(
      des_base,
      .wvw_benchmark_response = 2 * s2R * vals,
      .wvw_benchmark_rank = R_base
    )
    benchmark_fit <- survey::svyglm(
      .wvw_benchmark_response ~ .wvw_benchmark_rank,
      design = benchmark_design,
      family = stats::gaussian()
    )
    benchmark_coefficients <- stats::coef(benchmark_fit)
    benchmark_vcov <- stats::vcov(benchmark_fit)
    if (
      !(".wvw_benchmark_rank" %in% names(
        benchmark_coefficients
      )) ||
        !(".wvw_benchmark_rank" %in% rownames(
          benchmark_vcov
        ))
    ) {
      return(list(
        estimate = NA_real_,
        std_error = NA_real_,
        coefficients = benchmark_coefficients,
        vcov = benchmark_vcov,
        degrees_of_freedom = survey::degf(benchmark_design),
        conditioned_on = paste(
          "survey convenient regression was rank deficient"
        )
      ))
    }
    list(
      estimate = unname(
        scale_M *
          benchmark_coefficients[".wvw_benchmark_rank"]
      ),
      std_error = unname(
        abs(scale_M) * sqrt(
          benchmark_vcov[
            ".wvw_benchmark_rank",
            ".wvw_benchmark_rank"
          ]
        )
      ),
      coefficients = benchmark_coefficients,
      vcov = benchmark_vcov,
      degrees_of_freedom = survey::degf(benchmark_design),
      conditioned_on = paste(
        "fractional ranks, rank variance, outcome mean, and",
        "correction scale under the declared survey design"
      )
    )
  }
  conindex_total_survey_benchmark <-
    conindex_survey_benchmark(Y_base, M_base)

  mod_base <- suppressWarnings(stats::glm.fit(
    x = X_full,
    y = Y_base,
    weights = W_base,
    family = fam_link,
    control = fit_control
  ))
  if (!mod_base$converged) {
    diagnostics$convergence <- "CRITICAL MODEL VALIDATION: GLM algorithm did not converge."
    if (!relax) stop(diagnostics$convergence) else if (!quiet) warning(diagnostics$convergence, call. = FALSE)
  }

  b_star <- mod_base$coefficients; b_star[is.na(b_star)] <- 0; names(b_star) <- c_names

  if (model_type %in% c("logit", "probit")) {
    perfectly_classified <-
      any(abs(b_star[is.finite(b_star)]) > 15) &&
      all(
        (mod_base$fitted.values >= 0.5) ==
          (Y_base == 1)
      )
    if (perfectly_classified) {
      stop(
        paste0(
          "CRITICAL MODEL VALIDATION: perfect separation detected. ",
          "Finite logit/probit coefficients do not exist for this ",
          "estimation sample."
        ),
        call. = FALSE
      )
    }
  }

  if (model_type %in% c("logit", "probit") && any(abs(b_star) > 15)) {
    diagnostics$separation <- "Validation Warning: Extreme coefficients detected (|beta| > 15). Possible perfect separation (Hauck-Donner effect)."
    if (!relax && !quiet) warning(diagnostics$separation, call. = FALSE)
  }

  AME_base <- get_ames(b_star, X_full, W_base, model_type, assign_attr)
  mu_x_base <- colSums(w_norm * X_full)
  GC_X_base <- calc_cov_r_vec(X_full, R_base, W_base)

  # NA (not 0) when the determinant mean is ~0: the concentration index is undefined
  # there (display only; contributions use GC_X directly and are unaffected).
  CI_X_base <- ifelse(abs(mu_x_base) < 1e-12, NA_real_, GC_X_base / mu_x_base)

  Contribs <- M_base * AME_base * GC_X_base; Contribs[1] <- 0

  y_hat <- as.vector(X_full %*% b_star)
  if (model_type == "logit") y_hat <- stats::plogis(y_hat) else if (model_type == "probit") y_hat <- stats::pnorm(y_hat)
  res_fit <- Y_base - y_hat
  True_Resid <- M_base * calc_cov_r(res_fit, R_base, W_base)
  Approx_Err <- Total_Index - sum(Contribs) - True_Resid

  # === 6. VARIANCE ESTIMATION (VCE) ===
  complex_survey_design <- use_svy &&
    (!is.null(psu_var) || !is.null(strata_var))
  design_df <- if (complex_survey_design) {
    max(1, survey::degf(des_base))
  } else {
    Inf
  }
  V_CC <- matrix(0, nrow = length(c_names), ncol = length(c_names), dimnames = list(c_names, c_names))
  model_vcov_audit <- NULL
  model_design_coefficients_audit <- NULL
  model_vcov_type_audit <- NULL
  replication_audit <- NULL
  linearization_audit <- NULL
  V_all <- NULL

  if (vce_method == "linearized") {
    eta_base <- as.vector(X_full %*% b_star)
    if (model_type == "ols") {
      score_scalar <- Y_base - y_hat
      information_scalar <- rep(1, N_obs)
    } else {
      probability <- if (model_type == "logit") {
        stats::plogis(eta_base)
      } else {
        stats::pnorm(eta_base)
      }
      probability <- pmin(pmax(probability, 1e-15), 1 - 1e-15)
      binary_variance <- probability * (1 - probability)
      link_derivative <- if (model_type == "logit") {
        binary_variance
      } else {
        stats::dnorm(eta_base)
      }
      response_residual <- Y_base - probability
      score_scalar <- response_residual * link_derivative /
        binary_variance
      information_scalar <- link_derivative^2 / binary_variance
      if (model_type == "probit") {
        derivative_prime <- -eta_base * link_derivative
        information_scalar <-
          link_derivative^2 / binary_variance -
          response_residual * derivative_prime / binary_variance +
          response_residual * link_derivative^2 *
            (1 - 2 * probability) / binary_variance^2
      }
    }

    information <- crossprod(
      X_full,
      X_full * (w_norm * information_scalar)
    )
    beta_bread <- MASS::ginv(information)
    beta_influence <- (X_full * score_scalar) %*% beta_bread
    beta_influence <- sweep(
      beta_influence,
      2,
      colSums(w_norm * beta_influence),
      "-"
    )
    colnames(beta_influence) <- c_names

    ame_linearization <- get_ame_linearization(
      betas = b_star,
      mm = X_full,
      norm_weights = w_norm,
      m_type = model_type,
      assign_idx = assign_attr,
      beta_influence = beta_influence
    )
    empirical_rank_values <- if (is.null(precalc_rank_var)) {
      as.numeric(df_clean[[ses_var]])
    } else {
      NULL
    }
    determinant_gc_linearization <- get_gc_influence(
      values = X_full,
      ranks = R_base,
      norm_weights = w_norm,
      empirical_rank_values = empirical_rank_values
    )
    outcome_gc_linearization <- get_gc_influence(
      values = matrix(
        Y_base,
        ncol = 1,
        dimnames = list(NULL, dep_var)
      ),
      ranks = R_base,
      norm_weights = w_norm,
      empirical_rank_values = empirical_rank_values
    )
    residual_gc_linearization <- get_gc_influence(
      values = matrix(
        res_fit,
        ncol = 1,
        dimnames = list(NULL, "Residual")
      ),
      ranks = R_base,
      norm_weights = w_norm,
      empirical_rank_values = empirical_rank_values
    )

    scale_derivative <- get_scale_derivative(
      mu_y_base,
      correction,
      dep_min,
      dep_max
    )
    mean_y_influence <- Y_base - mu_y_base
    contribution_influence <-
      sweep(
        ame_linearization$influence,
        2,
        M_base * GC_X_base,
        "*"
      ) +
      sweep(
        determinant_gc_linearization$influence,
        2,
        M_base * AME_base,
        "*"
      ) +
      outer(
        mean_y_influence,
        scale_derivative * AME_base * GC_X_base
      )
    contribution_influence[, 1] <- 0
    contribution_influence <- sweep(
      contribution_influence,
      2,
      colSums(w_norm * contribution_influence),
      "-"
    )
    colnames(contribution_influence) <- c_names

    total_influence <-
      M_base * outcome_gc_linearization$influence[, 1] +
      scale_derivative *
        outcome_gc_linearization$estimate[1] *
        mean_y_influence

    prediction_derivative <- if (model_type == "ols") {
      X_full
    } else {
      prediction_density <- if (model_type == "logit") {
        stats::dlogis(eta_base)
      } else {
        stats::dnorm(eta_base)
      }
      X_full * prediction_density
    }
    residual_beta_gradient <- -calc_cov_r_vec(
      prediction_derivative,
      R_base,
      W_base
    )
    residual_influence <-
      M_base * (
        residual_gc_linearization$influence[, 1] +
          as.vector(beta_influence %*% residual_beta_gradient)
      ) +
      scale_derivative *
        residual_gc_linearization$estimate[1] *
        mean_y_influence

    component_names <- c(c_names, "Total Index", "Residual")
    component_influence <- cbind(
      contribution_influence,
      `Total Index` = total_influence,
      Residual = residual_influence
    )
    component_influence <- sweep(
      component_influence,
      2,
      colSums(w_norm * component_influence),
      "-"
    )
    colnames(component_influence) <- component_names

    vcov_from_influence <- function(influence, component_labels) {
      influence <- as.matrix(influence)
      linearized_design <- des_base
      generated_names <- paste0(
        ".wvw_if_",
        seq_len(ncol(influence))
      )
      for (column in seq_along(generated_names)) {
        linearized_design$variables[[generated_names[column]]] <-
          influence[, column] / sum(W_base)
      }
      total_formula <- stats::reformulate(generated_names)
      total_fit <- survey::svytotal(
        total_formula,
        linearized_design
      )
      covariance <- stats::vcov(total_fit)
      dimnames(covariance) <- list(
        component_labels,
        component_labels
      )
      covariance
    }

    V_all <- vcov_from_influence(
      component_influence,
      component_names
    )
    V_CC <- V_all[
      seq_along(c_names),
      seq_along(c_names),
      drop = FALSE
    ]
    SE_Cont <- safe_se(diag(V_CC))
    SE_Cont[1] <- 0
    se_expl_ov <- safe_se(sum(V_CC))
    se_total_ov <- safe_se(
      V_all["Total Index", "Total Index"]
    )
    se_resid_ov <- safe_se(V_all["Residual", "Residual"])

    model_vcov_audit <- vcov_from_influence(
      beta_influence,
      c_names
    )
    if (!complex_survey_design && model_type == "ols") {
      model_vcov_audit <- (
        (N_obs - 1) / (N_obs - ncol(X_full))
      ) * model_vcov_audit
      model_vcov_type_audit <- paste(
        "Huber-White HC1 N/(N-k), matching Stata regress,",
        "from the joint score influence"
      )
    } else if (!complex_survey_design) {
      model_vcov_type_audit <- paste(
        "Huber-White ML sandwich N/(N-1), matching Stata",
        paste0(model_type, ", vce(robust)")
      )
    } else {
      model_vcov_type_audit <- paste(
        "Taylor-linearized observed-score sandwich under the",
        "declared survey design"
      )
    }
    model_design_coefficients_audit <- b_star

    conindex_residual_benchmark <-
      conindex_robust_benchmark(res_fit, M_base)
    approximation_influence <- total_influence -
      rowSums(contribution_influence) -
      residual_influence
    linearization_row_names <- as.character(
      df_clean[[source_row_name]]
    )
    rownames(component_influence) <- linearization_row_names
    rownames(beta_influence) <- linearization_row_names
    rownames(ame_linearization$direct_influence) <-
      linearization_row_names
    rownames(ame_linearization$influence) <-
      linearization_row_names
    rownames(determinant_gc_linearization$influence) <-
      linearization_row_names
    names(approximation_influence) <-
      linearization_row_names
    linearization_audit <- list(
      rank_treated_as =
        determinant_gc_linearization$rank_treated_as,
      finite_sample_correction = if (complex_survey_design) {
        "survey Taylor design covariance"
      } else {
        "N/(N-1) covariance of the joint influence function"
      },
      correction_scale_derivative = scale_derivative,
      component_names = component_names,
      influence_function = component_influence,
      linearized_total_variable =
        component_influence / sum(W_base),
      coefficient_influence = beta_influence,
      marginal_effect_gradient_beta =
        ame_linearization$gradient_beta,
      marginal_effect_direct_influence =
        ame_linearization$direct_influence,
      marginal_effect_influence =
        ame_linearization$influence,
      determinant_gc_influence =
        determinant_gc_linearization$influence,
      approximation_error_influence =
        approximation_influence,
      conindex_robust_benchmark = list(
        total = conindex_total_benchmark,
        residual = conindex_residual_benchmark
      ),
      conindex_survey_benchmark = list(
        total = conindex_total_survey_benchmark,
        residual = conindex_survey_benchmark(
          res_fit,
          M_base
        )
      ),
      vcov = V_all
    )

  } else {
    ordinary_resampling <- !use_svy
    des_rep_type <- NA_character_
    survey_generator_scale <- NA_real_
    survey_effective_scale <- NA_real_
    survey_replicate_type <- NA_character_
    survey_replicate_engine <- NA_character_
    replicate_center_groups <- NULL

    if (ordinary_resampling && vce_method == "bootstrap") {
      # Generate ordinary observation-level samples with boot::boot(). The
      # existing WVW replicate pipeline below continues to recompute the model,
      # ranks, marginal effects, correction, and complete result vector from
      # the resulting integer frequency weights. simple = TRUE preserves the
      # established seed-to-replicate mapping.
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
      replicate_factors <- replicate_weights
      generator_scale <- 1 / (boot_reps - 1)
      replicate_rscales <- rep(1, boot_reps)
      generator_type <-
        "ordinary multinomial bootstrap by analytic observation"
      design_df <- Inf
    } else if (ordinary_resampling) {
      # resample::jackknife() supplies the deterministic leave-one-out
      # samples. Preserve the caller's RNG state around the call because the
      # package initializes .Random.seed when it does not already exist.
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
      replicate_factors <- replicate_weights
      generator_scale <- (N_obs - 1) / N_obs
      replicate_rscales <- rep(1, N_obs)
      generator_type <- "ordinary delete-one jackknife"
      design_df <- max(1, N_obs - 1)
    } else {
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
          survey_replicate_engine <- "svrep::as_bootstrap_design"
          svrep::as_bootstrap_design(
            des_base,
            type = "Rao-Wu-Yue-Beaumont",
            replicates = boot_reps,
            mse = FALSE,
            samp_method_by_stage = "SRSWR"
          )
        } else {
          survey_replicate_engine <- "survey::as.svrepdesign"
          survey::as.svrepdesign(
            des_base,
            type = "subbootstrap",
            replicates = boot_reps
          )
        }
      } else {
        survey_replicate_engine <- "survey::as.svrepdesign"
        survey::as.svrepdesign(
          des_base,
          type = des_rep_type
        )
      }
      replicate_factors <- as.matrix(des_rep$repweights)
      replicate_weights <- replicate_factors * W_base
      generator_scale <- des_rep$scale
      replicate_rscales <- des_rep$rscales
      survey_generator_scale <- des_rep$scale
      survey_replicate_type <- if (
        vce_method == "bootstrap"
      ) {
        if (
          bootstrap_has_singletons &&
            bootstrap_singleton == "certainty"
        ) {
          "Rao-Wu-Yue-Beaumont bootstrap with singleton certainty"
        } else {
          "Rao-Wu rescaled bootstrap"
        }
      } else {
        des_rep_type
      }
      generator_type <- survey_replicate_type
      survey_effective_scale <- if (
        vce_method == "bootstrap"
      ) {
        1 / ncol(replicate_factors)
      } else {
        generator_scale
      }
      if (vce_method == "jackknife" && !is.null(strata_var)) {
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
        if (!anyNA(candidate_groups)) {
          replicate_center_groups <- candidate_groups
        }
      }
      design_df <- max(1, n_psu - n_strata)
    }
    n_replicates <- ncol(replicate_weights)
    colnames(replicate_weights) <- colnames(replicate_factors) <-
      paste0("replicate_", seq_len(n_replicates))
    rownames(replicate_weights) <- rownames(replicate_factors) <-
      as.character(df_clean[[source_row_name]])
    if (!quiet) {
      cat(sprintf(
        "\n[*] Running %s variance estimation...\n",
        toupper(vce_method)
      ))
      pb <- utils::txtProgressBar(
        min = 0,
        max = n_replicates,
        style = 3
      )
      env_pb <- new.env()
      env_pb$count <- 0
    }

    n_extra <- 2L  # append per-replicate Total index and Residual for overall SEs
    calc_rep <- function(w, data_r) {
      if (!quiet) {
        env_pb$count <- env_pb$count + 1
        utils::setTxtProgressBar(
          pb,
          min(env_pb$count, n_replicates)
        )
      }
      # NA (not 0) so degenerate replicates are dropped instead of entering the variance as (0 - estimate)^2
      if (sum(w) == 0) return(rep(NA_real_, length(c_names) + n_extra))
      idx_v <- w > 0
      w_sub <- w[idx_v]; y_sub <- Y_base[idx_v]; X_sub <- X_full[idx_v, , drop = FALSE]

      tryCatch({
        m_r <- suppressWarnings(stats::glm.fit(
          x = X_sub,
          y = y_sub,
          weights = w_sub,
          family = fam_link,
          control = fit_control
        ))
        if (!m_r$converged) return(rep(NA, length(c_names) + n_extra))

        R_rep <- get_r_i(data_r[idx_v, , drop=FALSE], w_sub, if(!is.null(precalc_rank_var)) data_r[[precalc_rank_var]][idx_v] else NULL, if(!is.null(ses_var)) ses_var else names(data_r)[1])
        M_rep <- get_scale_m(sum(w_sub * y_sub) / sum(w_sub), correction, dep_min, dep_max)
        GC_X_rep <- calc_cov_r_vec(X_sub, R_rep, w_sub)

        cfs <- m_r$coefficients; cfs[is.na(cfs)] <- 0
        if (
          model_type %in% c("logit", "probit") &&
            any(abs(cfs[is.finite(cfs)]) > 15) &&
            all(
              (m_r$fitted.values >= 0.5) ==
                (y_sub == 1)
            )
        ) {
          return(rep(NA_real_, length(c_names) + n_extra))
        }
        Cont_rep <- M_rep * get_ames(cfs, X_sub, w_sub, model_type, assign_attr) * GC_X_rep
        Cont_rep[1] <- 0

        Total_rep <- M_rep * calc_cov_r(y_sub, R_rep, w_sub)
        yhat_r <- as.vector(X_sub %*% cfs)
        if (model_type == "logit") yhat_r <- stats::plogis(yhat_r) else if (model_type == "probit") yhat_r <- stats::pnorm(yhat_r)
        Resid_rep <- M_rep * calc_cov_r(y_sub - yhat_r, R_rep, w_sub)
        return(c(Cont_rep, .__total__ = Total_rep, .__resid__ = Resid_rep))
      }, error = function(e) return(rep(NA, length(c_names) + n_extra)))
    }

    replicate_estimates <- vapply(
      seq_len(n_replicates),
      function(rep_index) {
        calc_rep(
          replicate_weights[, rep_index],
          df_clean
        )
      },
      numeric(length(c_names) + n_extra)
    )
    res_matrix <- list(
      theta = c(Contribs, Total_Index, True_Resid),
      replicates = t(replicate_estimates)
    )
    if (!quiet) close(pb)
    replicate_component_names <- c(
      c_names,
      "Total Index",
      "Residual"
    )
    colnames(res_matrix$replicates) <- replicate_component_names
    valid_idx <- stats::complete.cases(res_matrix$replicates)

    fail_rate <- sum(!valid_idx) / length(valid_idx)
    if (!ordinary_resampling && any(!valid_idx) && !relax) {
      stop(
        "CRITICAL: ",
        sum(!valid_idx),
        " survey replicate(s) failed. Dropping a survey replicate can ",
        "invalidate the design variance. Inspect the sparse or separated ",
        "replicates, or rerun with relax = TRUE for an explicitly ",
        "exploratory approximation.",
        call. = FALSE
      )
    }
    if (fail_rate > 0.05) {
      diagnostics$resampling <- sprintf("Statistical Alert: %.1f%% of resampling iterations failed to converge. SEs might be less reliable.", fail_rate * 100)
      if (!relax && fail_rate > 0.2) stop(diagnostics$resampling) else if (!quiet) warning(diagnostics$resampling, call. = FALSE)
    }
    if (sum(valid_idx) < 2) stop("CRITICAL: Resampling variance estimation completely failed due to sparsity.")

    kc <- length(c_names)
    centers_all <- c(Contribs, Total_Index, True_Resid)
    names(centers_all) <- replicate_component_names
    valid_replicates <- res_matrix$replicates[
      valid_idx,
      ,
      drop = FALSE
    ]
    valid_center_groups <- if (
      !is.null(replicate_center_groups)
    ) {
      replicate_center_groups[valid_idx]
    } else {
      NULL
    }
    if (is.null(valid_center_groups)) {
      replicate_center <- colMeans(valid_replicates)
      replicate_center_method <- "replicate mean"
      diffs_all <- sweep(
        valid_replicates,
        2,
        replicate_center,
        "-"
      )
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
      replicate_center <- do.call(
        rbind,
        lapply(valid_center_groups, function(group) {
          group_centers[[group]]
        })
      )
      colnames(replicate_center) <-
        replicate_component_names
      replicate_center_method <-
        "stratum-specific replicate mean"
      diffs_all <- valid_replicates - replicate_center
    }
    replicate_mse <- FALSE
    rscales <- if (
      length(replicate_rscales) == length(valid_idx)
    ) {
      replicate_rscales[valid_idx]
    } else {
      replicate_rscales
    }
    if (ordinary_resampling) {
      # Stata's non-MSE prefixes calculate the sample covariance of the
      # nonmissing replicate estimates. Consequently their scale is determined
      # by the number of valid, rather than requested, replications.
      valid_count <- sum(valid_idx)
      effective_scale <- if (vce_method == "bootstrap") {
        1 / (valid_count - 1)
      } else {
        (valid_count - 1) / valid_count
      }
      rep_adjust <- effective_scale / generator_scale
      if (vce_method == "jackknife") {
        design_df <- max(1, valid_count - 1)
      }
    } else {
      # With relax = TRUE this is only an exploratory survey safeguard. The
      # default survey path stops above instead of silently dropping a PSU
      # replicate.
      rep_adjust <- length(valid_idx) / sum(valid_idx)
      effective_scale <- survey_effective_scale * rep_adjust
    }
    if (any(!valid_idx)) {
      diagnostics$rep_adjust <- sprintf(
        paste0(
          "Variance Note: %d failed replicate(s) dropped; the effective ",
          "replicate scale is %.8g."
        ),
        sum(!valid_idx),
        effective_scale
      )
    }
    V_all <- effective_scale *
      crossprod(diffs_all * sqrt(rscales))
    V_CC <- V_all[seq_len(kc), seq_len(kc), drop = FALSE]
    se_total_ov <- safe_se(V_all[kc + 1, kc + 1])
    se_resid_ov <- safe_se(V_all[kc + 2, kc + 2])
    se_expl_ov  <- safe_se(sum(V_CC))
    rownames(V_CC) <- colnames(V_CC) <- c_names
    SE_Cont <- safe_se(diag(V_CC)); SE_Cont[1] <- 0

    rownames(V_all) <- colnames(V_all) <- replicate_component_names
    replication_audit <- list(
      method = vce_method,
      generator_type = generator_type,
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
      current_center_method = replicate_center_method,
      center_method = replicate_center_method,
      center = replicate_center,
      mse = replicate_mse,
      full_sample = stats::setNames(centers_all, replicate_component_names),
      estimates = res_matrix$replicates,
      replicates = res_matrix$replicates,
      valid = valid_idx,
      requested_replicates = n_replicates,
      valid_replicates = sum(valid_idx),
      failed_replicates = sum(!valid_idx),
      replicate_factors = replicate_factors,
      replicate_weights = replicate_weights,
      generator_scale = generator_scale,
      effective_scale = effective_scale,
      scale = effective_scale,
      rscales = replicate_rscales,
      valid_rscales = rscales,
      failed_replicate_adjustment = rep_adjust,
      vcov = V_all,
      pending_methodological_review = FALSE
    )
  }

  # In replicate modes the point estimator is based on glm.fit(). Fit a
  # design-aware model only for auditability; failures do not alter the
  # decomposition or its replicate variance.
  if (is.null(model_vcov_audit)) {
    model_audit_fit <- tryCatch(
      suppressWarnings(
        survey::svyglm(form_mod, design = des_base, family = fam_link)
      ),
      error = function(e) NULL
    )
    if (!is.null(model_audit_fit)) {
      model_vcov_audit <- stats::vcov(model_audit_fit)
      model_design_coefficients_audit <- stats::coef(model_audit_fit)
      model_vcov_type_audit <-
        "survey::svyglm design-based covariance (audit only)"
    }
  }

  # === 7. INFERENTIAL METRICS & DYNAMIC GROUPING ===
  stat_t <- Contribs / ifelse(SE_Cont == 0, 1e-12, SE_Cont)
  prob_t <- 2 * stats::pt(
    abs(stat_t),
    df = design_df,
    lower.tail = FALSE
  )

  # Grouped table harmonized with oby_decomp/f_decomp: 'Term' label and 'Estimate'
  # as the contribution, plus the shared inferential columns. 'Dummies' (number of
  # model-matrix columns aggregated) is a wvw-specific extra that is preserved.
  tb_grouped <- tibble::tibble(Term = character(), Estimate = numeric(), `% Contribution Explained` = numeric(), `% Contribution Total` = numeric(), Std_Error = numeric(), Statistic = numeric(), P_Value = numeric(), Conf_Low = numeric(), Conf_High = numeric(), Dummies = integer())

  q_val <- 1 - ((1 - level) / 2)

  for (i in seq_along(final_groups)) {
    var_name <- names(final_groups)[i]
    g_cols <- final_groups[[var_name]]
    g_contrib <- sum(Contribs[g_cols])

    se_g <- safe_se(sum(V_CC[g_cols, g_cols]))

    g_stat_t <- g_contrib / ifelse(is.na(se_g) | se_g == 0, NA_real_, se_g)
    g_prob_t <- 2 * stats::pt(
      abs(g_stat_t),
      df = design_df,
      lower.tail = FALSE
    )

    tb_grouped <- rbind(tb_grouped, data.frame(
        Term = var_name,
        Estimate = g_contrib,
        `% Contribution Explained` = ifelse(abs(sum(Contribs)) > 1e-12, (g_contrib / sum(Contribs)) * 100, 0),
        `% Contribution Total` = ifelse(abs(Total_Index) > 1e-12, (g_contrib / Total_Index) * 100, 0),
        Std_Error = se_g, Statistic = g_stat_t, P_Value = g_prob_t,
        Conf_Low = g_contrib - stats::qt(q_val, df = design_df) * se_g,
        Conf_High = g_contrib + stats::qt(q_val, df = design_df) * se_g,
        Dummies = length(g_cols),
        check.names = FALSE
    ))
  }

  tb_grouped <- tibble::as_tibble(tb_grouped)
  tb_grouped <- tb_grouped %>% dplyr::arrange(dplyr::desc(abs(.data$Estimate)))

  grouped_transform <- matrix(
    0,
    nrow = length(final_groups),
    ncol = length(c_names),
    dimnames = list(names(final_groups), c_names)
  )
  if (length(final_groups) > 0) {
    for (group_name in names(final_groups)) {
      grouped_transform[group_name, final_groups[[group_name]]] <- 1
    }
  }
  grouped_vcov_audit <- grouped_transform %*% V_CC %*%
    t(grouped_transform)

  # === 8. RESULTS ASSEMBLY & OUTPUT (CRAN STYLE) ===
  display_terms <- unname(sapply(c_names, function(x) ifelse(x %in% names(display_dict), display_dict[[x]], x)))
  explained_index <- sum(Contribs)
  elasticity_vec <- if (abs(mu_y_base) > 1e-12) (AME_base * mu_x_base) / mu_y_base else rep(NA_real_, length(Contribs))
  perc_explained_vec <- if (abs(explained_index) > 1e-12) (Contribs / explained_index) * 100 else rep(0, length(Contribs))
  perc_total_vec <- if (abs(Total_Index) > 1e-12) (Contribs / Total_Index) * 100 else rep(0, length(Contribs))

  # 'Estimate' is the determinant's contribution to the index (harmonized column
  # name); the wvw-specific chain (Marginal Effects, Mean_X, Elasticity,
  # Concentration Index) is preserved ahead of it.
  tb_detailed <- tibble::tibble(
    Term = display_terms,
    `Marginal Effects` = as.numeric(AME_base),
    Mean_X = as.numeric(mu_x_base),
    Elasticity = as.numeric(elasticity_vec),
    `Concentration Index` = as.numeric(CI_X_base),
    Estimate = as.numeric(Contribs),
    `% Contribution Explained` = as.numeric(perc_explained_vec),
    `% Contribution Total` = as.numeric(perc_total_vec),
    Std_Error = as.numeric(SE_Cont),
    Statistic = as.numeric(stat_t),
    P_Value = as.numeric(prob_t),
    Conf_Low = as.numeric(
      Contribs - stats::qt(q_val, df = design_df) * SE_Cont
    ),
    Conf_High = as.numeric(
      Contribs + stats::qt(q_val, df = design_df) * SE_Cont
    )
  )

  Metric_Name <- dplyr::case_when(correction == "erreygers" ~ "Total Index: Erreygers (Bounded Absolute)", correction == "wagstaff" ~ "Total Index: Wagstaff (Bounded Relative)", correction == "generalized" ~ "Total Index: Generalized (Unbounded Absolute)", correction == "standard" ~ "Total Index: Standard (Unbounded Relative)")
  Resid_Label <- if (model_type == "ols") "Unexplained (Residual)" else "Unexplained (Residual from prediction errors)"

  # results_overall harmonized with oby_decomp/f_decomp: full inference columns,
  # plus a wvw-specific 'Percentage' (share of the total index) that is preserved.
  ov_terms <- c(Metric_Name, "Explained (Sum of Contributions)", Resid_Label)
  ov_est   <- c(Total_Index, sum(Contribs), True_Resid)
  ov_se    <- as.numeric(c(se_total_ov, se_expl_ov, se_resid_ov))
  ov_pct   <- c(100, ifelse(abs(Total_Index) > 1e-12, (sum(Contribs)/Total_Index)*100, 0), ifelse(abs(Total_Index) > 1e-12, (True_Resid/Total_Index)*100, 0))
  if (model_type != "ols") {
    ov_terms <- c(ov_terms, "Non-linear Approximation Error")
    ov_est   <- c(ov_est, Approx_Err)
    ov_se    <- c(ov_se, NA_real_)
    ov_pct   <- c(ov_pct, ifelse(abs(Total_Index) > 1e-12, (Approx_Err/Total_Index)*100, 0))
  }
  ov_stat <- ov_est / ifelse(is.na(ov_se) | ov_se == 0, NA_real_, ov_se)
  tb_overall <- tibble::tibble(
    Term = ov_terms,
    Estimate = ov_est,
    Percentage = ov_pct,
    Std_Error = ov_se,
    Statistic = ov_stat,
    P_Value = 2 * stats::pt(
      abs(ov_stat),
      df = design_df,
      lower.tail = FALSE
    ),
    Conf_Low = ov_est - stats::qt(q_val, df = design_df) * ov_se,
    Conf_High = ov_est + stats::qt(q_val, df = design_df) * ov_se
  )

  SS_tot <- sum(W_base * (Y_base - mu_y_base)^2); r2 <- if (SS_tot > 0) 1 - (sum(W_base * res_fit^2) / SS_tot) else NA
  fit_metrics <- tibble::tibble(Model = paste0(if (use_svy) "Survey-weighted Wagstaff decomp" else "Wagstaff decomp (Unweighted)", " - ", toupper(model_type)), Pseudo_R_Squared = r2)

  if (!is.null(V_all)) {
    kc <- length(c_names)
    replicate_component_names <- colnames(V_all)
    overall_transform <- matrix(
      0,
      nrow = length(ov_terms),
      ncol = length(replicate_component_names),
      dimnames = list(ov_terms, replicate_component_names)
    )
    overall_transform[ov_terms[1], kc + 1] <- 1
    overall_transform[ov_terms[2], seq_len(kc)] <- 1
    overall_transform[ov_terms[3], kc + 2] <- 1
    if (length(ov_terms) == 4) {
      overall_transform[ov_terms[4], kc + 1] <- 1
      overall_transform[ov_terms[4], seq_len(kc)] <- -1
      overall_transform[ov_terms[4], kc + 2] <- -1
    }
    overall_vcov_audit <- overall_transform %*% V_all %*%
      t(overall_transform)
    overall_vcov_complete <- TRUE
  } else {
    overall_vcov_audit <- matrix(
      NA_real_,
      nrow = length(ov_terms),
      ncol = length(ov_terms),
      dimnames = list(ov_terms, ov_terms)
    )
    diag(overall_vcov_audit) <- ov_se^2
    overall_vcov_complete <- FALSE
  }

  if (overall_vcov_complete) {
    ov_se <- as.numeric(safe_se(diag(overall_vcov_audit)))
    ov_stat <- ov_est /
      ifelse(is.na(ov_se) | ov_se == 0, NA_real_, ov_se)
    tb_overall$Std_Error <- ov_se
    tb_overall$Statistic <- ov_stat
    tb_overall$P_Value <- 2 * stats::pt(
      abs(ov_stat),
      df = design_df,
      lower.tail = FALSE
    )
    tb_overall$Conf_Low <- ov_est -
      stats::qt(q_val, df = design_df) * ov_se
    tb_overall$Conf_High <- ov_est +
      stats::qt(q_val, df = design_df) * ov_se
  }

  source_rows_sorted <- df_clean[[source_row_name]]
  model_matrix_audit <- X_full
  rownames(model_matrix_audit) <- as.character(source_rows_sorted)
  determinant_terms <- vapply(
    assign_attr,
    function(term_index) {
      if (term_index == 0) "(Intercept)" else term_labels[term_index]
    },
    character(1)
  )
  determinant_audit <- data.frame(
    column = c_names,
    source_term = determinant_terms,
    display_term = display_terms,
    marginal_effect = as.numeric(AME_base),
    mean_x = as.numeric(mu_x_base),
    generalized_concentration_covariance = as.numeric(GC_X_base),
    concentration_index = as.numeric(CI_X_base),
    elasticity = as.numeric(elasticity_vec),
    contribution = as.numeric(Contribs),
    stringsAsFactors = FALSE
  )
  grouped_estimates_audit <- vapply(
    final_groups,
    function(columns) sum(Contribs[columns]),
    numeric(1)
  )
  rank_method_audit <- if (is.null(precalc_rank_var)) {
    "generalized weighted fractional midrank (ties grouped)"
  } else {
    "precalculated rank supplied by user"
  }
  rank_audit <- data.frame(
    source_row = source_rows_sorted,
    sort_value = as.numeric(df_clean[[sort_var]]),
    fractional_rank = as.numeric(R_base),
    weight = as.numeric(W_base),
    stringsAsFactors = FALSE
  )

  raw <- list(
    schema_version = "1.0",
    settings = list(
      dep_var = dep_var,
      indep_vars = indep_vars,
      ses_var = ses_var,
      precalc_rank_var = precalc_rank_var,
      correction = correction,
      model_type = model_type,
      dep_min = dep_min,
      dep_max = dep_max,
      is_shortfall = is_shortfall,
      use_svy = use_svy,
      weight_var = weight_var,
      strata_var = strata_var,
      psu_var = psu_var,
      vce_method = vce_method,
      boot_reps = boot_reps,
      lonely_psu = lonely_psu,
      bootstrap_singleton = bootstrap_singleton,
      level = level,
      seed = seed,
      relax = relax,
      rank_method = rank_method_audit,
      numeric_rounding = "none"
    ),
    sample = list(
      input_n = input_n,
      analytic_n = N_obs,
      analytic_source_rows = source_rows_sorted,
      complete_case_source_rows_before_sort = complete_case_rows,
      excluded_source_rows = list(
        missing_required_values = missing_excluded_rows,
        zero_sampling_weight = zero_weight_excluded_rows
      ),
      sort_variable = sort_var
    ),
    design = list(
      weights = stats::setNames(
        as.numeric(W_base),
        as.character(source_rows_sorted)
      ),
      normalized_weights = stats::setNames(
        as.numeric(w_norm),
        as.character(source_rows_sorted)
      ),
      strata = stats::setNames(
        as.character(df_clean$svy_strata),
        as.character(source_rows_sorted)
      ),
      psu = stats::setNames(
        as.character(df_clean$svy_psu),
        as.character(source_rows_sorted)
      ),
      population_weight = pop_size,
      strata_n = n_strata,
      psu_n = n_psu,
      degrees_of_freedom = design_df
    ),
    rank = rank_audit,
    model = list(
      formula = form_mod,
      matrix = model_matrix_audit,
      assign = assign_attr,
      term_labels = term_labels,
      coefficients = b_star,
      design_model_coefficients = model_design_coefficients_audit,
      vcov = model_vcov_audit,
      vcov_type = model_vcov_type_audit,
      fitted = stats::setNames(
        as.numeric(y_hat),
        as.character(source_rows_sorted)
      ),
      residuals = stats::setNames(
        as.numeric(res_fit),
        as.character(source_rows_sorted)
      ),
      converged = isTRUE(mod_base$converged)
    ),
    index = list(
      outcome_mean = mu_y_base,
      rank_mean = sum(w_norm * R_base),
      rank_variance = s2R,
      generalized_concentration_covariance_y =
        calc_cov_r(Y_base, R_base, W_base),
      correction_scale = M_base,
      total = Total_Index,
      explained = explained_index,
      residual = True_Resid,
      nonlinear_approximation_error = Approx_Err
    ),
    determinants = determinant_audit,
    groups = final_groups,
    estimates = list(
      overall = stats::setNames(ov_est, ov_terms),
      detailed = stats::setNames(Contribs, c_names),
      grouped = grouped_estimates_audit
    ),
    vcov = list(
      model = model_vcov_audit,
      detailed = V_CC,
      grouped = grouped_vcov_audit,
      overall = overall_vcov_audit,
      detailed_complete = !is.null(V_all),
      grouped_complete = !is.null(V_all),
      overall_complete = overall_vcov_complete
    ),
    linearization = linearization_audit,
    replication = replication_audit,
    known_limitations = list(
      precalculated_rank = if (
        vce_method == "linearized" &&
          !is.null(precalc_rank_var)
      ) {
        paste(
          "The supplied precalculated rank is conditioned on because its",
          "generating variables and sampling design are unavailable."
        )
      } else {
        NULL
      },
      survey_design = if (use_svy) {
        paste(
          "The built-in design represents one final weight, optional strata,",
          "and one PSU stage under an ultimate-cluster/with-replacement",
          "variance approximation. FPC, lower-stage identifiers, first-stage",
          "PPS resampling, recalibration, and external replicate weights are",
          "not represented."
        )
      } else {
        NULL
      }
    )
  )

  if (!quiet) {
    cat("\n")
    cat(rep("-", 80), "\n", sep="")
    cat("WAGSTAFF DECOMPOSITION OF HEALTH INEQUALITY (RANK-BASED)\n")
    cat(rep("-", 80), "\n", sep="")
    cat(sprintf("Number of obs   = %-15d Metric Corrected = %s\n", N_obs, toupper(correction)))
    cat(sprintf("Population size = %-15.4f VCE Engine       = %s\n", pop_size, toupper(vce_method)))
    cat(sprintf("Design df       = %-15s Model Link       = %s\n", ifelse(is.infinite(design_df), "Inf", design_df), toupper(model_type)))
    cat(sprintf("Number of Strata= %-15s Index Value      = %.5f\n", ifelse(use_svy, n_strata, "N/A"), Total_Index))
    cat(sprintf("Number of PSUs  = %-15s \n", ifelse(use_svy, n_psu, "N/A")))
    cat(rep("-", 80), "\n", sep="")

    if (length(diagnostics) > 0) {
        cat("\n[DIAGNOSTICS & WARNINGS]\n")
        for (diag_name in names(diagnostics)) cat(" *", diagnostics[[diag_name]], "\n")
    } else {
        cat("\n[DIAGNOSTICS: ALL ECONOMETRIC ASSUMPTIONS PASSED]\n")
    }

    round_for_print <- function(table, digits = 4) {
      table <- as.data.frame(table)
      numeric_columns <- vapply(table, is.numeric, logical(1))
      table[numeric_columns] <- lapply(
        table[numeric_columns],
        round,
        digits = digits
      )
      table
    }

    cat("\n[OVERALL INEQUALITY SUMMARY]\n"); print(round_for_print(tb_overall, 5), row.names = FALSE)
    cat("\n[MODEL FIT METRICS]\n"); print(round_for_print(fit_metrics, 4), row.names = FALSE)
    cat(sprintf("\n[%s STRUCTURAL CONTRIBUTIONS]\n", ifelse(length(groupings) > 0, "CUSTOM GROUPED", "AUTO-GROUPED")))
    print(round_for_print(tb_grouped, 4), row.names = FALSE)

    cat("\n[DETAILED DECOMPOSITION]\n")
    print(round_for_print(tb_detailed[-1, c("Term", "Marginal Effects", "Elasticity", "Concentration Index", "Estimate", "% Contribution Explained", "% Contribution Total", "Std_Error", "Statistic", "P_Value")], 4), row.names = FALSE)

    cat("\nNote: In detailed tables, categorical predictors indicate '(Category vs. Ref: BaseCategory)'\n")
    cat("      and continuous predictors indicate '(Continuous)'.\n")
    cat(rep("-", 80), "\n\n", sep="")
  }

  return(invisible(list(
    summary_stats = list(
      N = N_obs,
      Pop = pop_size,
      Strata = n_strata,
      PSUs = n_psu,
      DF = design_df
    ),
    model_metrics = fit_metrics,
    diagnostics = diagnostics,
    results_overall = tb_overall,
    results_grouped = tb_grouped,
    results_detailed = tb_detailed,
    raw = raw
  )))
}
