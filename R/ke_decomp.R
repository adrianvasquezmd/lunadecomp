#' Direct Regression Decomposition of Socioeconomic Health Inequality
#'
#' @description
#' Implements the direct regression approach proposed by Kessels and Erreygers
#' for decomposing socioeconomic inequality of health. The method reformulates a
#' rank-dependent or level-dependent bivariate inequality index as the mean of an
#' individual-level composite variable, then regresses that composite variable on
#' explanatory variables.
#'
#' Unlike the conventional Wagstaff decomposition, this function does not
#' decompose the index by multiplying health-regression coefficients by
#' concentration indices of the covariates. Instead, it follows the direct
#' regression logic: the estimated coefficients are interpreted as marginal
#' effects of the explanatory variables on the individual composite outcome that
#' underlies the inequality index.
#'
#' @details
#' Let `dep_var` denote the health variable and `ses_var` denote the
#' socioeconomic ranking or socioeconomic level variable. The function first
#' constructs an individual-level target variable whose weighted mean equals the
#' selected socioeconomic inequality index. It then estimates a direct regression
#' model:
#'
#' \deqn{
#'   T_i = \alpha + X_i\beta + \varepsilon_i
#' }
#'
#' where `T_i` is the constructed rank- or level-dependent target variable.
#'
#' For `index_type = "rank"`, the target is based on the relative fractional
#' socioeconomic rank. This corresponds to the rank-dependent composite variable
#' used by Kessels and Erreygers for indices such as the concentration index and
#' corrected rank-dependent indices.
#'
#' For `index_type = "level"`, the target is based on relative socioeconomic
#' levels rather than ranks. This corresponds to a level-dependent bivariate
#' inequality index, where the weights depend on socioeconomic levels directly.
#'
#' @section Interpretation:
#' The coefficients reported in `results_detailed` are **marginal effects** from
#' the direct regression of the composite target variable. They should not be
#' interpreted as additive percentage contributions to the observed index in the
#' same way as in the Wagstaff decomposition. In the Kessels-Erreygers framework,
#' variable importance is assessed primarily through statistical evidence from
#' the direct regression, especially F tests and logworth values.
#'
#' The logworth statistic is defined as:
#'
#' \deqn{
#'   -\log_{10}(p)
#' }
#'
#' where `p` is the p value from the corresponding F test. Larger logworth values
#' indicate stronger evidence that a variable or group of dummy variables is
#' important in explaining variation in the direct-regression target.
#'
#' @section Index corrections:
#' The argument `correction` controls the scaling applied to the underlying
#' bivariate index:
#'
#' \itemize{
#'   \item `"generalized"`: generalized absolute index.
#'   \item `"standard"`: relative index obtained by scaling by the mean of the
#'   health variable.
#'   \item `"erreygers"`: bounded correction for cardinal bounded outcomes. For
#'   rank-dependent indices the scaling constant is 4; for level-dependent
#'   indices it is 1, as specified by Kessels and Erreygers.
#'   \item `"wagstaff"`: rank-dependent bounded relative normalization motivated
#'   by the mean-dependent feasible bounds of the concentration index. It is not
#'   available for `index_type = "level"` because the published Wagstaff
#'   normalization does not define a level-dependent index and the analogous
#'   scalar does not normalize its feasible bounds.
#' }
#'
#' When using `"erreygers"` or `"wagstaff"`, both theoretical health bounds must
#' be supplied through `dep_min` and `dep_max`. Sample minima and maxima are not
#' substituted for theoretical bounds because doing so changes the estimand
#' across samples. All analyzed outcome values must lie within the declared
#' bounds.
#'
#' The aggregate Erreygers index is invariant to positive affine
#' transformations when the theoretical bounds are transformed consistently,
#' and bounded rank indices satisfy the usual mirror property. The individual
#' KE target, however, acquires an additional mean-zero rank or level term under
#' a translation or outcome complementation. Consequently, the aggregate index
#' can remain invariant (or change sign under mirroring) while its direct
#' regression coefficients change. This is a property of the KE coefficient
#' decomposition, not a numerical discrepancy.
#'
#' @section Survey designs:
#' If `use_svy = TRUE`, the function uses the supplied survey weights, strata,
#' and primary sampling units. Supplying any of `weight_var`, `strata_var`, or
#' `psu_var` while `use_svy = FALSE` activates survey mode automatically.
#' Without a weight variable, unit sampling weights are used. Without a PSU
#' variable, each analytic observation is treated as its own PSU. Design-based
#' t and F inference uses the design degrees of freedom, PSUs minus strata,
#' matching Stata `svy`.
#'
#' `weight_var` is interpreted as a sampling/design weight. The same normalized
#' weights define the weighted socioeconomic ranks or relative levels, health
#' mean, index, covariate means, and WLS point regression. For positive integer
#' weights, this point estimator is exactly the estimator obtained by expanding
#' each row by its weight and applying the original equal-weight KE method.
#' Stata `pweight`, `aweight`, and `fweight` regressions share these WLS point
#' coefficients for a common positive weight vector, but their statistical
#' meanings and variance estimators are not interchangeable. This function does
#' not reinterpret a sampling weight as an analytic or frequency weight; survey
#' inference follows the declared design.
#'
#' @section Variance estimation:
#' Kessels and Erreygers (2019) define the point index, direct-regression
#' coefficients, t/F tests, and logworth, but do not specify a covariance
#' estimator or report inference for the aggregate index. Accordingly, every
#' VCE option in this function applies only to the direct-regression
#' coefficients. The aggregate index is always reported as a point estimate;
#' its inferential fields remain `NA`.
#'
#' With ordinary independent observations, `"linearized"` reports the HC1
#' Huber--White sandwich covariance used by Stata
#' `regress, vce(robust)`, conditional on the empirically constructed KE target.
#' This is an explicit package extension, not a covariance estimator prescribed
#' by the article. It does not propagate uncertainty from ranks, health or
#' socioeconomic means, or the correction scale.
#'
#' In survey mode, `"linearized"` analogously reports the Taylor design
#' covariance of the direct-regression coefficients conditional on the
#' constructed target. It accounts for declared final sampling weights,
#' strata, and first-stage PSUs, but does not linearize the construction of
#' ranks, relative levels, means, or correction scaling. This conditional
#' survey covariance is also an explicit package extension.
#'
#' Kessels and Erreygers (2019, equations 12 and 17) provide influence
#' functions for the rank- and level-dependent *indices* in their separate RIF
#' discussion, but not a joint influence function for the coefficients of the
#' direct regression. Those formulas are not used to add index inference to the
#' direct method.
#'
#' Outside survey mode, bootstrap uses ordinary observation-level samples
#' generated by [boot::boot()] and jackknife uses the deterministic delete-one
#' samples from [resample::jackknife()]. The complete KE target, including ranks
#' or relative levels, means, index scaling, and regression, is recomputed inside
#' every replicate, but only the coefficient vector enters the replicate
#' covariance. Ordinary non-MSE covariance follows Stata's prefix convention:
#' bootstrap uses the sample covariance of the valid coefficient replicates and
#' jackknife uses \eqn{(N-1)/N} times their centered cross-product. When a
#' replicate returns a non-finite coefficient vector, the complete replicate is
#' omitted and `N` in these formulas becomes the number of valid replicates,
#' matching the non-MSE missing-replicate convention of Stata's prefixes. The
#' number and identity of omitted replicates remain available in `raw`.
#'
#' Survey replicate methods remain distinct: they use design replicate weights
#' rather than resampling analytic rows. Survey jackknife uses JK1 for
#' unstratified designs and JKn for stratified designs. JK1 is centered on the
#' replicate mean; JKn is centered separately on the mean of the replicates
#' that alter each stratum, matching Stata's default non-MSE pseudovalue
#' convention.
#'
#' Survey bootstrap uses the Rao--Wu rescaled construction implemented by
#' `survey::as.svrepdesign(type = "subbootstrap")`: within each stratum,
#' \eqn{n_h-1} PSUs are sampled with replacement and their weights are
#' multiplied by \eqn{n_h/(n_h-1)}. The coefficient covariance is centered on
#' the replicate mean and uses Stata's `svy bootstrap` finite-replicate scale
#' \eqn{1/B}. The native `survey` generator scale \eqn{1/(B-1)} is retained in
#' `raw$replication$generator_scale` for auditability.
#'
#' Both survey replicate methods use PSUs minus strata degrees of freedom.
#' They require at least two PSUs in every stratum and represent an
#' ultimate-cluster, with-replacement design with one final sampling weight.
#' They are not a first-stage PPS bootstrap and do not represent FPCs,
#' calibration replicate weights, or lower sampling stages. A failed survey
#' replicate stops estimation by default; `relax = TRUE` permits an explicitly
#' exploratory complete-replicate deletion and rescaling.
#'
#' @param data A data frame containing the health outcome, socioeconomic
#' variable, covariates, and optional survey design variables.
#' @param dep_var Character string. Name of the health variable.
#' @param indep_vars Character vector. Names of the explanatory variables to
#' include in the direct regression model. Each element may also be an R formula
#' term such as `"age * education"` or `"education:region"`. Categorical
#' predictors should be coded as factors or character variables; their first
#' factor level is the reference unless the user relevels them before calling
#' the function.
#' @param ses_var Optional character string. Name of the socioeconomic variable
#' used to construct ranks or relative levels. Required unless the relevant
#' precalculated variable is supplied.
#' @param index_type Character string. Type of bivariate inequality index.
#' Options are `"rank"` and `"level"`.
#' @param correction Character string. Index scaling/correction. Options are
#' `"erreygers"`, `"wagstaff"`, `"standard"`, and `"generalized"`.
#' @param precalc_rank_var Optional character string. Name of a precalculated
#' fractional rank variable for `index_type = "rank"`. Values should typically
#' lie between 0 and 1. When supplied, these positions take precedence over
#' `ses_var` and remain fixed inside jackknife and bootstrap replicates; they
#' are repeated or reweighted with the sampled records rather than recomputed.
#' @param precalc_level_var Optional character string. Name of a precalculated
#' relative socioeconomic level variable for `index_type = "level"`, such as
#' income divided by mean income. Values must be finite and nonnegative, with
#' at least one positive value. As with a precalculated rank, the supplied
#' positions remain fixed inside all replication methods.
#' @param dep_min Optional numeric. Theoretical lower bound of the health
#' variable. Required for `"erreygers"` and `"wagstaff"`.
#' @param dep_max Optional numeric. Theoretical upper bound of the health
#' variable. Required for `"erreygers"` and `"wagstaff"`.
#' @param use_svy Logical. If `TRUE`, applies complex survey design settings.
#' @param weight_var Optional character string. Name of the sampling/design
#' weight variable. At the point-estimation stage, positive weights are
#' normalized by their total and used consistently for ranks or levels, means,
#' the index, and the WLS direct regression.
#' @param strata_var Optional character string. Name of the stratification
#' variable.
#' @param psu_var Optional character string. Name of the primary sampling unit
#' variable.
#' @param vce_method Character string. Variance estimation method. Options are
#' `"linearized"`, `"jackknife"`, and `"bootstrap"`.
#' @param boot_reps Integer. Number of bootstrap replications when
#' `vce_method = "bootstrap"`.
#' @param lonely_psu Character string. Taylor-linearization handling of strata
#' with one PSU. Options are `"fail"`, `"remove"`, `"certainty"`, `"adjust"`,
#' and `"average"`. Survey jackknife and bootstrap instead require at least two
#' PSUs in every stratum.
#' @param level Numeric. Confidence level for confidence intervals. Default is
#' `0.95`.
#' @param seed Optional integer. Random seed for reproducibility of replicate
#' resampling. Defaults to `NULL` (non-deterministic); set an explicit value to
#' reproduce results. When an explicit seed is supplied, the prior state of the
#' global random number generator is restored on exit. Ordinary jackknife is
#' deterministic and also leaves the caller's RNG state unchanged. Bootstrap
#' with `seed = NULL` follows normal R behavior and advances the caller's random
#' number stream.
#' @param relax Logical. If `TRUE`, selected validation problems are stored as
#' diagnostics rather than stopping execution.
#' @param quiet Logical. If `TRUE`, suppresses console output and progress bars.
#'
#' @return
#' A named list. The inferential columns (`Term`, `Estimate`, `Std_Error`,
#' `Statistic`, `P_Value`, `Conf_Low`, `Conf_High`) are harmonized with
#' [oby_decomp()], [f_decomp()], and [wvw_decomp()]; method-specific columns are
#' preserved. Consistent with Kessels and Erreygers (2019), no additive
#' contribution or grouped-contribution tables are produced, because in the
#' direct approach the coefficients are marginal effects, not contributions.
#' \describe{
#'   \item{summary_stats}{Sample size, population size, number of strata, number
#'   of PSUs, and design degrees of freedom.}
#'   \item{model_metrics}{Model label and R-squared from the direct regression.}
#'   \item{diagnostics}{Warnings and diagnostic messages generated during data
#'   validation, model estimation, or variance estimation.}
#'   \item{results_overall}{The selected index point estimate and health bounds
#'   (`Min_Health`, `Max_Health`). The harmonized inferential columns
#'   (`Std_Error`, `Statistic`, `P_Value`, confidence interval) are always `NA`:
#'   the article does not report aggregate-index inference, and the package VCE
#'   extensions are intentionally limited to the direct-regression
#'   coefficients.}
#'   \item{results_anova}{Variable-level F tests and logworth values used to
#'   rank the importance of variables or factor blocks.}
#'   \item{results_detailed}{Detailed marginal effects by model term: `Estimate`
#'   (the marginal effect), `Std_Error`, `Statistic`, t-test `P_Value`, and
#'   confidence interval, plus the KE-specific F-test p value (`Prob_F`) and
#'   `Logworth` importance measure. `Omitted` identifies coefficients excluded
#'   because of exact collinearity; these follow Stata's zero-coefficient,
#'   zero-variance display convention.}
#'   \item{raw}{Unrounded audit information, including the analytic sample,
#'   normalized weights, rank or relative-level position, individual KE target,
#'   model matrix, coefficients, covariance matrices, fitted values, residuals,
#'   index identities, configuration, and replicate estimates when applicable.}
#' }
#'
#' @examples
#' data(lunadecomp_example)
#'
#' fit_ke <- ke_decomp(
#'   data = lunadecomp_example,
#'   dep_var = "health_score",
#'   indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
#'   ses_var = "ses",
#'   index_type = "rank",
#'   correction = "erreygers",
#'   dep_min = 0,
#'   dep_max = 1,
#'   vce_method = "linearized",
#'   quiet = TRUE
#' )
#'
#' fit_ke$results_overall
#' fit_ke$results_anova
#' fit_ke$results_detailed
#'
#' @author
#' Adrian Vasquez-Mejia, MD, MSc \cr
#' Oscar J. Mujica, MD, MPH, PHE, FACE \cr
#' Antonio Sanhueza, MPH, MSc, PhD
#'
#' @references
#' Kessels, R., & Erreygers, G. (2019). A direct regression approach to
#' decomposing socioeconomic inequality of health. \emph{Health Economics},
#' 28(7), 884--905. \doi{10.1002/hec.3891}
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
#' @importFrom dplyr %>% select all_of mutate mutate_if row_number case_when arrange desc
#' @importFrom tidyr drop_na
#' @importFrom tibble tibble add_row as_tibble
#' @importFrom rlang .data
#' @importFrom MASS ginv
#' @importFrom survey svydesign svyglm as.svrepdesign withReplicates
#' @importFrom stats as.formula coef cov.wt lm.wfit pf pt qt vcov
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
ke_decomp <- function(
    data, dep_var, indep_vars, ses_var = NULL,
    index_type = "rank", correction = "erreygers",
    precalc_rank_var = NULL, precalc_level_var = NULL,
    dep_min = NULL, dep_max = NULL, use_svy = FALSE, weight_var = NULL,
    strata_var = NULL, psu_var = NULL, vce_method = "linearized",
    boot_reps = 500, lonely_psu = "adjust", level = 0.95,
    seed = NULL, relax = FALSE, quiet = FALSE
) {

  is_scalar_character <- function(value) {
    is.character(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      nzchar(trimws(value))
  }
  is_optional_variable <- function(value) {
    is.null(value) || is_scalar_character(value)
  }
  is_scalar_logical <- function(value) {
    is.logical(value) &&
      length(value) == 1L &&
      !is.na(value)
  }

  if (!is.data.frame(data)) {
    stop("Error: data must be a data.frame or a data.frame subclass.")
  }
  if (nrow(data) < 1L) {
    stop("Error: data must contain at least one observation.")
  }
  if (!is_scalar_character(dep_var)) {
    stop("Error: dep_var must be one non-empty character string.")
  }
  variable_arguments <- list(
    ses_var = ses_var,
    precalc_rank_var = precalc_rank_var,
    precalc_level_var = precalc_level_var,
    weight_var = weight_var,
    strata_var = strata_var,
    psu_var = psu_var
  )
  invalid_variable_arguments <- names(variable_arguments)[
    !vapply(variable_arguments, is_optional_variable, logical(1))
  ]
  if (length(invalid_variable_arguments) > 0L) {
    stop(
      "Error: ",
      paste(invalid_variable_arguments, collapse = ", "),
      " must be NULL or one non-empty character string.",
      call. = FALSE
    )
  }
  for (argument_name in c("use_svy", "relax", "quiet")) {
    if (!is_scalar_logical(get(argument_name))) {
      stop(
        "Error: ", argument_name,
        " must be exactly TRUE or FALSE.",
        call. = FALSE
      )
    }
  }
  if (
    !is.numeric(level) ||
      length(level) != 1L ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1
  ) {
    stop("Error: level must be one finite number strictly between 0 and 1.")
  }
  if (
    !is.null(seed) &&
      (
        !is.numeric(seed) ||
          length(seed) != 1L ||
          !is.finite(seed) ||
          seed < 0 ||
          seed > .Machine$integer.max ||
          seed != as.integer(seed)
      )
  ) {
    stop(
      "Error: seed must be NULL or one integer between 0 and ",
      .Machine$integer.max, ".",
      call. = FALSE
    )
  }
  valid_vce <- c("linearized", "jackknife", "bootstrap")
  valid_corr <- c("erreygers", "wagstaff", "standard", "generalized")
  valid_lonely_psu <- c(
    "fail", "remove", "certainty", "adjust", "average"
  )
  if (
    !is_scalar_character(vce_method) ||
      !(vce_method %in% valid_vce)
  ) {
    stop(
      "Error: vce_method must be 'linearized', 'jackknife', or 'bootstrap'."
    )
  }
  if (
    !is_scalar_character(correction) ||
      !(correction %in% valid_corr)
  ) {
    stop(
      paste0(
        "Error: correction must be 'erreygers', 'wagstaff', ",
        "'standard', or 'generalized'."
      )
    )
  }
  if (
    !is_scalar_character(index_type) ||
      !(index_type %in% c("rank", "level"))
  ) {
    stop("Error: index_type must be 'rank' or 'level'.")
  }
  if (
    !is_scalar_character(lonely_psu) ||
      !(lonely_psu %in% valid_lonely_psu)
  ) {
    stop(
      "Error: lonely_psu must be 'fail', 'remove', 'certainty', ",
      "'adjust', or 'average'."
    )
  }
  if (index_type == "rank" && !is.null(precalc_level_var)) {
    stop(
      "Error: precalc_level_var is only valid when index_type = 'level'."
    )
  }
  if (index_type == "level" && !is.null(precalc_rank_var)) {
    stop(
      "Error: precalc_rank_var is only valid when index_type = 'rank'."
    )
  }

  # C1: restore the caller's RNG state on exit instead of leaving it altered.
  # resample::jackknife() can initialize .Random.seed despite being
  # deterministic, so preserve an absent state for ordinary jackknife as well.
  ordinary_jackknife_requested <-
    identical(vce_method, "jackknife") &&
    !isTRUE(use_svy) &&
    is.null(weight_var)
  if (!is.null(seed) || ordinary_jackknife_requested) {
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

  if (correction == "wagstaff" && index_type != "rank") {
    stop(
      paste0(
        "Error: correction = 'wagstaff' is defined only for ",
        "index_type = 'rank'. No published Wagstaff normalization ",
        "exists for the level-dependent KE index."
      )
    )
  }
  if (
    vce_method == "bootstrap" &&
      (
        !is.numeric(boot_reps) ||
        length(boot_reps) != 1L ||
          !is.finite(boot_reps) ||
          boot_reps < 2 ||
          boot_reps != as.integer(boot_reps)
      )
  ) {
    stop("Error: boot_reps must be a single integer greater than or equal to 2.")
  }

  if (index_type == "rank" && is.null(ses_var) && is.null(precalc_rank_var)) {
    stop("Error: For rank-dependent index, provide either 'ses_var' or 'precalc_rank_var'.")
  }
  if (index_type == "level" && is.null(ses_var) && is.null(precalc_level_var)) {
    stop("Error: For level-dependent index, provide either 'ses_var' or 'precalc_level_var'.")
  }

  active_ses_var <- if (
    (index_type == "rank" && is.null(precalc_rank_var)) ||
      (index_type == "level" && is.null(precalc_level_var))
  ) {
    ses_var
  } else {
    NULL
  }
  active_precalc_var <- if (index_type == "rank") {
    precalc_rank_var
  } else {
    precalc_level_var
  }

  if (
    !is.character(indep_vars) ||
      length(indep_vars) < 1L ||
      anyNA(indep_vars) ||
      any(!nzchar(trimws(indep_vars)))
  ) {
    stop(
      "Error: indep_vars must contain one or more non-empty variable names or formula terms."
    )
  }
  form_base <- tryCatch(
    stats::as.formula(
      paste("~", paste(indep_vars, collapse = " + ")),
      env = parent.frame()
    ),
    error = function(e) {
      stop(
        "Error: indep_vars could not be parsed as model terms: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  model_variables <- all.vars(form_base)
  if (length(model_variables) < 1L) {
    stop("Error: indep_vars does not identify any explanatory variable.")
  }

  old_lonely_psu <- getOption("survey.lonely.psu")
  options(survey.lonely.psu = lonely_psu)
  on.exit(options(survey.lonely.psu = old_lonely_psu), add = TRUE)

  req_vars <- unique(c(
    dep_var, active_ses_var, active_precalc_var,
    model_variables, weight_var, strata_var, psu_var
  ))
  req_vars <- req_vars[!is.na(req_vars) & req_vars != ""]
  missing_variables <- setdiff(req_vars, names(data))
  if (length(missing_variables) > 0L) {
    stop(
      "Error: the following required variables are absent from data: ",
      paste(missing_variables, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  input_n <- nrow(data)
  source_row_name <- ".__ke_source_row__"
  selected_data <- data %>%
    dplyr::select(dplyr::all_of(req_vars))
  complete_case_mask <- stats::complete.cases(selected_data)
  missing_excluded_rows <- which(!complete_case_mask)
  selected_data[[source_row_name]] <- seq_len(input_n)
  df_clean <- selected_data[complete_case_mask, , drop = FALSE] %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    droplevels()

  N_obs <- nrow(df_clean)
  if (!quiet && N_obs < nrow(data)) cat(sprintf("\n[!] Data Prep: %d rows with NAs were dropped. N = %d\n", nrow(data) - N_obs, N_obs))
  if (N_obs == 0) stop("CRITICAL: Zero observations remaining after NA removal.")

  coerce_finite_numeric <- function(variable, role) {
    converted <- suppressWarnings(
      as.numeric(as.character(df_clean[[variable]]))
    )
    invalid <- which(!is.finite(converted))
    if (length(invalid) > 0L) {
      source_rows <- df_clean[[source_row_name]][invalid]
      stop(
        "Error: ", role, " '", variable,
        "' must contain only finite numeric values; invalid source row(s): ",
        paste(utils::head(source_rows, 8L), collapse = ", "),
        if (length(source_rows) > 8L) ", ..." else "",
        ".",
        call. = FALSE
      )
    }
    converted
  }
  df_clean[[dep_var]] <- coerce_finite_numeric(
    dep_var, "dep_var"
  )
  if (!is.null(active_ses_var)) {
    df_clean[[active_ses_var]] <- coerce_finite_numeric(
      active_ses_var, "ses_var"
    )
  }
  if (!is.null(active_precalc_var)) {
    df_clean[[active_precalc_var]] <- coerce_finite_numeric(
      active_precalc_var,
      if (index_type == "rank") {
        "precalc_rank_var"
      } else {
        "precalc_level_var"
      }
    )
  }
  if (use_svy && !is.null(weight_var)) {
    df_clean[[weight_var]] <- coerce_finite_numeric(
      weight_var, "weight_var"
    )
    if (any(df_clean[[weight_var]] <= 0)) {
      stop(
        "Error: weight_var must contain strictly positive sampling weights.",
        call. = FALSE
      )
    }
    if (!is.finite(sum(df_clean[[weight_var]]))) {
      stop(
        "Error: the total sampling weight must be finite.",
        call. = FALSE
      )
    }
  }
  for (variable in model_variables) {
    values <- df_clean[[variable]]
    if (
      (is.numeric(values) || is.integer(values)) &&
        any(!is.finite(values))
    ) {
      stop(
        "Error: model variable '", variable,
        "' must contain only finite values.",
        call. = FALSE
      )
    }
    if (is.factor(values) && nlevels(values) < 2L) {
      stop(
        "Error: factor model variable '", variable,
        "' has fewer than two observed levels after missing-value removal.",
        call. = FALSE
      )
    }
  }
  for (design_variable in c(strata_var, psu_var)) {
    if (
      !is.null(design_variable) &&
        any(!nzchar(trimws(as.character(df_clean[[design_variable]]))))
    ) {
      stop(
        "Error: design variable '", design_variable,
        "' contains an empty label.",
        call. = FALSE
      )
    }
  }

  sort_var <- if (index_type == "rank") {
    if (!is.null(precalc_rank_var)) precalc_rank_var else ses_var
  } else {
    if (!is.null(precalc_level_var)) precalc_level_var else ses_var
  }
  if (!is.null(sort_var) && index_type == "rank") df_clean <- df_clean[order(df_clean[[sort_var]]), ]

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
  diagnostics <- list()

  if (
    use_svy &&
      vce_method %in% c("jackknife", "bootstrap")
  ) {
    unique_psu <- unique(
      df_clean[, c("svy_strata", "svy_psu"), drop = FALSE]
    )
    psu_per_stratum <- table(unique_psu$svy_strata)
    if (any(psu_per_stratum < 2L)) {
      stop(
        "Error: survey replication requires at least two PSUs in every ",
        "stratum. Use Taylor linearization with an explicit lonely_psu ",
        "rule, collapse or redesign singleton strata, or supply an ",
        "appropriate external replicate-weight design.",
        call. = FALSE
      )
    }
    diagnostics$survey_replication_scope <- paste(
      "Survey Replication Scope: ultimate-cluster variance using one",
      "final weight, optional strata, and one PSU stage; no FPC,",
      "lower-stage identifiers, calibration replicate weights, or",
      "first-stage PPS resampling is represented."
    )
    if (!quiet) message(diagnostics$survey_replication_scope)
  }

  # W3: SEs from a variance without masking non-positive terms via abs().
  safe_se <- function(v) {
    neg <- which(v < -1e-8)
    if (length(neg) > 0) {
      diagnostics$neg_var <<- sprintf(
        "Numerical Warning: %d non-positive variance term(s) detected (non-PSD covariance); their SE set to NA.",
        length(neg)
      )
      if (!relax && !quiet) warning(diagnostics$neg_var)
    }
    v[v < 0] <- NA_real_
    sqrt(v)
  }

  diagnostics$overall_inference <- paste(
    "Article Scope: the aggregate KE index is reported as a point",
    "estimate only. Covariance extensions apply to direct-regression",
    "coefficients, so index SE, test, and confidence interval remain NA."
  )

  if (vce_method == "linearized") {
    diagnostics$vce <- paste(
      "Methodological Warning: 'linearized' VCE is the coefficient",
      "covariance conditional on the empirically constructed KE target.",
      "It does not jointly propagate ranks or relative levels, means, or",
      "correction scaling. It is a package extension and is not specified",
      "by Kessels and Erreygers (2019). Use 'jackknife' or 'bootstrap'",
      "when coefficient inference should reconstruct the complete target."
    )
    if (!quiet) warning(diagnostics$vce)
  }

  if (correction %in% c("erreygers", "wagstaff")) {
    if (is.null(dep_min) || is.null(dep_max)) {
      stop(
        paste0(
          "CRITICAL: correction = '", correction,
          "' requires explicit theoretical health bounds in both ",
          "'dep_min' and 'dep_max'. Observed sample extrema are not ",
          "valid substitutes."
        )
      )
    }
    if (
      length(dep_min) != 1L || length(dep_max) != 1L ||
        !is.numeric(dep_min) || !is.numeric(dep_max) ||
        !is.finite(dep_min) || !is.finite(dep_max)
    ) {
      stop(
        "CRITICAL: 'dep_min' and 'dep_max' must be finite numeric scalars."
      )
    }
    if (dep_min >= dep_max) {
      stop("CRITICAL: dep_var bounds are invalid (dep_min >= dep_max).")
    }
    observed_health_range <- range(
      df_clean[[dep_var]],
      na.rm = TRUE
    )
    bound_tolerance <- sqrt(.Machine$double.eps) *
      max(1, abs(dep_min), abs(dep_max))
    if (
      observed_health_range[1] < dep_min - bound_tolerance ||
        observed_health_range[2] > dep_max + bound_tolerance
    ) {
      stop(
        sprintf(
          paste0(
            "CRITICAL: dep_var has observed values outside the ",
            "declared theoretical bounds [%.17g, %.17g]; observed ",
            "range is [%.17g, %.17g]."
          ),
          dep_min, dep_max,
          observed_health_range[1], observed_health_range[2]
        )
      )
    }
  } else {
    dep_min <- NA; dep_max <- NA
  }

  # ==============================================================================
  # DYNAMIC TARGET GENERATOR (M-SCALAR EXPANDED FOR WAGSTAFF)
  # ==============================================================================
  calc_target <- function(
      h_var, s_var, w, idx_type, precalc_pos, min_h, max_h, corr_type,
      return_components = FALSE) {
    sigma <- w / sum(w, na.rm = TRUE)
    mu_h <- sum(sigma * h_var, na.rm = TRUE)
    mu_s <- NA_real_

    if (idx_type == "rank") {
      if (!is.null(precalc_pos)) {
        f_tilde <- 2 * precalc_pos
      } else {
        vals <- rle(as.numeric(s_var))
        agg_sigma <- tapply(sigma, rep(seq_along(vals$lengths), vals$lengths), sum)
        cum_sigma <- cumsum(agg_sigma)
        lag_cum_sigma <- c(0, cum_sigma[-length(cum_sigma)])
        f_tilde_group <- agg_sigma + 2 * lag_cum_sigma
        f_tilde <- rep(f_tilde_group, vals$lengths)
      }
      u_i <- (f_tilde * h_var) - mu_h
      position <- f_tilde / 2
      target_base <- u_i

      # DYNAMIC SCALAR RESOLUTION
      if (corr_type == "erreygers") {
        M_val <- 4 / (max_h - min_h)
      } else if (corr_type == "wagstaff") {
        denom <- (max_h - mu_h) * (mu_h - min_h)
        # Silently returning 0 here would produce an all-zero composite (index = 0)
        # and a meaningless regression; stop instead (consistent with wvw_decomp).
        if (denom <= 1e-12) stop("CRITICAL: Wagstaff correction is unstable because the outcome mean is too close to its theoretical bounds.")
        M_val <- (max_h - min_h) / denom
      } else if (corr_type == "standard") {
        if (mu_h <= 1e-12) stop("CRITICAL: Mean of dep_var is zero, negative, or nearly zero. The 'standard' index requires a ratio-scale outcome with a strictly positive mean.")
        M_val <- 1 / mu_h
      } else if (corr_type == "generalized") {
        M_val <- 1
      }

    } else {
      if (!is.null(precalc_pos)) {
        y_tilde <- precalc_pos
      } else {
        mu_s <- sum(sigma * s_var, na.rm = TRUE)
        # KE (2019) assume a ratio-scale SES variable; relative levels s/mu_s are
        # meaningless (or divide by ~0) when the SES mean is zero or negative.
        if (mu_s <= 1e-12) stop("CRITICAL: mean of ses_var is zero or negative. The level-dependent index requires a ratio-scale SES variable with a strictly positive mean.")
        y_tilde <- s_var / mu_s
      }
      d_i <- (y_tilde * h_var) - mu_h
      position <- y_tilde
      target_base <- d_i

      # DYNAMIC SCALAR RESOLUTION
      if (corr_type == "erreygers") {
        M_val <- 1 / (max_h - min_h)
      } else if (corr_type == "wagstaff") {
        stop(
          paste0(
            "INTERNAL ERROR: Wagstaff normalization reached the ",
            "level-dependent target generator."
          )
        )
      } else if (corr_type == "standard") {
        if (mu_h <= 1e-12) stop("CRITICAL: Mean of dep_var is zero, negative, or nearly zero. The 'standard' index requires a ratio-scale outcome with a strictly positive mean.")
        M_val <- 1 / mu_h
      } else if (corr_type == "generalized") {
        M_val <- 1
      }
    }

    target <- as.numeric(target_base * M_val)
    if (!return_components) return(target)

    list(
      target = target,
      target_base = as.numeric(target_base),
      position = as.numeric(position),
      normalized_weight = as.numeric(sigma),
      health_mean = mu_h,
      socioeconomic_mean = mu_s,
      correction_scale = M_val
    )
  }

  precalc_vec <- if (index_type == "rank" && !is.null(precalc_rank_var)) {
    df_clean[[precalc_rank_var]]
  } else if (index_type == "level" && !is.null(precalc_level_var)) {
    df_clean[[precalc_level_var]]
  } else {
    NULL
  }

  # The rank composite uses f_tilde = 2 * rank and assumes fractional ranks in [0, 1].
  if (index_type == "rank" && !is.null(precalc_vec)) {
    rng_r <- range(precalc_vec, na.rm = TRUE)
    if (rng_r[1] < -1e-8 || rng_r[2] > 1 + 1e-8) {
      diagnostics$precalc_rank <- sprintf("Validation Warning: 'precalc_rank_var' has values outside [0, 1] (range %.4f to %.4f). The composite assumes fractional ranks; rescale before use.", rng_r[1], rng_r[2])
      stop(diagnostics$precalc_rank, call. = FALSE)
    }
  }
  if (index_type == "level" && !is.null(precalc_vec)) {
    if (any(precalc_vec < 0) || !any(precalc_vec > 0)) {
      stop(
        paste0(
          "Error: 'precalc_level_var' must be nonnegative and contain ",
          "at least one positive relative socioeconomic level."
        ),
        call. = FALSE
      )
    }
  }

  target_components <- calc_target(
    h_var = df_clean[[dep_var]],
    s_var = if(!is.null(active_ses_var)) df_clean[[active_ses_var]] else rep(1, N_obs),
    w = df_clean$svy_weight,
    idx_type = index_type,
    precalc_pos = precalc_vec,
    min_h = dep_min, max_h = dep_max, corr_type = correction,
    return_components = TRUE
  )
  df_clean$TARGET_Y <- target_components$target

  Index_Value <- sum((df_clean$svy_weight / sum(df_clean$svy_weight)) * df_clean$TARGET_Y)

  # ==============================================================================
  # DICTIONARY & DESIGN MATRIX (Bug Fix 5b)
  # ==============================================================================
  display_dict <- list(); auto_groups <- list()
  X_full <- stats::model.matrix(form_base, data = df_clean)
  c_names <- colnames(X_full)

  if (ncol(X_full) >= N_obs) {
    stop(sprintf("CRITICAL MODEL VALIDATION: Predictors (%d) exceeds or equals observations (%d). Unidentified model.", ncol(X_full), N_obs))
  }

  assign_attr <- attr(X_full, "assign")
  term_labels <- attr(stats::terms(form_base), "term.labels")

  for (i in seq_along(term_labels)) {
    term <- term_labels[i]
    cols <- which(assign_attr == i)
    if(length(cols) > 0) {
      auto_groups[[term]] <- cols

      term_variables <- all.vars(
        stats::as.formula(paste("~", term))
      )
      if (
        length(term_variables) == 1L &&
          (
            is.factor(df_clean[[term_variables]]) ||
              is.character(df_clean[[term_variables]])
          )
      ) {
        factor_name <- term_variables
        ref_level <- levels(df_clean[[factor_name]])[1]
        for (cn in c_names[cols]) {
          lvl <- substring(cn, nchar(factor_name) + 1)
          display_dict[[cn]] <- sprintf(
            "%s (%s vs. Ref: %s)",
            factor_name, lvl, ref_level
          )
        }
      } else if (length(term_variables) > 1L) {
        for (cn in c_names[cols]) {
          display_dict[[cn]] <- sprintf(
            "%s (Interaction: %s)",
            cn, term
          )
        }
      } else {
        for (cn in c_names[cols]) {
          display_dict[[cn]] <- sprintf("%s (Continuous)", term)
        }
      }
    }
  }
  display_dict[["(Intercept)"]] <- "Baseline Reference (Intercept)"

  if (ncol(X_full) > 2) {
    cor_X <- suppressWarnings(stats::cov.wt(X_full[, -1, drop = FALSE], wt = df_clean$svy_weight, cor = TRUE)$cor)
    if (!is.null(cor_X) && !any(is.na(cor_X))) {
      inv_cor <- tryCatch(solve(cor_X), error = function(e) NULL)
      if (is.null(inv_cor)) {
        diagnostics$collinearity <- paste(
          "Model Note: perfect collinearity detected;",
          "aliased coefficients are reported as omitted."
        )
        if (!quiet) message(diagnostics$collinearity)
      } else if (any(diag(inv_cor) > 10)) {
        diagnostics$vif <- "Multicollinearity Warning (VIF > 10)."
        if (!relax && !quiet) warning(diagnostics$vif)
      }
    }
  }

  # ==============================================================================
  # VARIANCE ESTIMATION ENGINES
  # ==============================================================================
  form_w <- stats::as.formula("~svy_weight")
  form_strata <- if (use_svy && !is.null(strata_var)) stats::as.formula("~svy_strata") else NULL
  form_psu <- if (use_svy && !is.null(psu_var)) stats::as.formula("~svy_psu") else ~1
  des_base <- survey::svydesign(ids = form_psu, strata = form_strata, weights = form_w, data = df_clean, nest = use_svy)
  form_mod <- stats::update(form_base, TARGET_Y ~ .)
  replication_audit <- NULL
  model_object <- NULL
  model_vcov_type <- NULL
  aliased_coefficients <- stats::setNames(
    rep(FALSE, length(c_names)),
    c_names
  )

  if (vce_method == "linearized") {
    if (!use_svy) {
      mod_ordinary <- suppressWarnings(stats::lm.wfit(
        x = X_full,
        y = df_clean$TARGET_Y,
        w = rep(1, N_obs)
      ))
      model_object <- mod_ordinary
      model_vcov_type <- paste(
        "Stata-style HC1 independent-observation sandwich covariance",
        "conditional on the constructed target"
      )
      b_ordinary <- mod_ordinary$coefficients
      aliased_coefficients <- is.na(b_ordinary)
      names(aliased_coefficients) <- c_names
      estimable_columns <- sort(
        mod_ordinary$qr$pivot[seq_len(mod_ordinary$rank)]
      )
      X_estimable <- X_full[, estimable_columns, drop = FALSE]
      bread <- solve(crossprod(X_estimable))
      residual <- as.numeric(mod_ordinary$residuals)
      hc0 <- bread %*%
        crossprod(X_estimable * residual) %*%
        bread
      hc1 <- N_obs / (N_obs - mod_ordinary$rank) * hc0
      b_svy <- b_ordinary
      V_svy <- hc1
      rownames(V_svy) <- colnames(V_svy) <-
        c_names[estimable_columns]
      design_df <- N_obs - mod_ordinary$rank
    } else {
      mod_svy <- suppressWarnings(
        survey::svyglm(form_mod, design = des_base)
      )
      model_object <- mod_svy
      model_vcov_type <-
        "survey design covariance conditional on the constructed target"
      b_svy <- stats::coef(mod_svy)
      V_svy <- stats::vcov(mod_svy)
      # Stata svy inference uses the design degrees of freedom (PSUs minus
      # strata), not svyglm's default model-residual denominator df.
      design_df <- survey::degf(des_base)
      aliased_coefficients[names(b_svy)] <- is.na(b_svy)
    }

    # Aliased (collinear) coefficients: zero them with a diagnostic instead of
    # letting NA propagate through every downstream table.
    if (any(is.na(b_svy))) {
      diagnostics$aliased <- sprintf(
        paste0(
          "Model Note: %d coefficient(s) aliased by exact ",
          "collinearity and reported as omitted with zero estimate ",
          "and variance."
        ),
        sum(is.na(b_svy))
      )
      if (!quiet) message(diagnostics$aliased)
      b_svy[is.na(b_svy)] <- 0
    }
    b_star <- rep(0, length(c_names)); names(b_star) <- c_names
    b_star[names(b_svy)] <- b_svy
    V_full <- matrix(0, nrow = length(c_names), ncol = length(c_names), dimnames = list(c_names, c_names))
    match_idx <- intersect(rownames(V_svy), c_names)
    V_full[match_idx, match_idx] <- V_svy[match_idx, match_idx]
    # The article reports the aggregate index as a point estimate and does not
    # specify its inference. Package VCE extensions are deliberately limited to
    # the coefficients of the direct regression.
    se_index <- NA_real_
  } else {
    ordinary_resampling <- !use_svy
    des_rep <- NULL
    survey_generator_scale <- NA_real_
    survey_effective_scale <- NA_real_
    survey_replicate_type <- NULL
    replicate_center_groups <- NULL
    replicate_factors_audit <- NULL
    if (ordinary_resampling) {
      n_replicates <- if (vce_method == "bootstrap") {
        as.integer(boot_reps)
      } else {
        N_obs
      }
      generator_scale <- if (vce_method == "bootstrap") {
        1 / (n_replicates - 1)
      } else {
        (N_obs - 1) / N_obs
      }
      generator_rscales <- rep(1, n_replicates)
      generator_type <- if (vce_method == "bootstrap") {
        "ordinary multinomial bootstrap by analytic observation"
      } else {
        "ordinary delete-one jackknife"
      }
    } else {
      des_rep <- if (vce_method == "bootstrap") {
        survey::as.svrepdesign(
          des_base,
          type = "subbootstrap",
          replicates = boot_reps
        )
      } else {
        jackknife_type <- if (is.null(form_strata)) {
          "JK1"
        } else {
          "JKn"
        }
        survey::as.svrepdesign(
          des_base,
          type = jackknife_type
        )
      }
      replicate_factors_audit <- as.matrix(des_rep$repweights)
      n_replicates <- length(des_rep$rscales)
      generator_scale <- des_rep$scale
      generator_rscales <- des_rep$rscales
      survey_replicate_type <- if (
        vce_method == "bootstrap"
      ) {
        "Rao-Wu rescaled bootstrap"
      } else {
        jackknife_type
      }
      generator_type <- survey_replicate_type
      survey_generator_scale <- des_rep$scale
      # Stata's svy bootstrap with one bootstrap sample represented by each
      # bsrweight() variable uses bsn(1), hence scale 1/B. The survey package
      # supplies the Rao--Wu replicate weights; only its finite-B scale
      # convention 1/(B-1) is replaced here.
      survey_effective_scale <- if (
        vce_method == "bootstrap"
      ) {
        1 / n_replicates
      } else {
        generator_scale
      }
      if (
        vce_method == "jackknife" &&
          !is.null(form_strata)
      ) {
        candidate_groups <- vapply(
          seq_len(ncol(replicate_factors_audit)),
          function(rep_index) {
            changed <- abs(
              replicate_factors_audit[, rep_index] - 1
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
        # Standard JKn alters one stratum per replicate. Stata's default
        # non-MSE pseudovalue formula is equivalent to centering each such
        # block on its own stratum-specific replicate mean.
        if (!anyNA(candidate_groups)) {
          replicate_center_groups <- candidate_groups
        }
      }
      design_df <- max(1, n_psu - n_strata)
    }
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
    m_base <- suppressWarnings(stats::lm.wfit(x = X_full, y = df_clean$TARGET_Y, w = df_clean$svy_weight))
    model_object <- m_base
    model_vcov_type <- paste0(
      vce_method,
      " replicate covariance of the complete recomputed target"
    )
    b_star <- m_base$coefficients
    aliased_coefficients <- is.na(b_star)
    names(aliased_coefficients) <- c_names
    if (any(is.na(b_star))) {
      diagnostics$aliased <- sprintf(
        paste0(
          "Model Note: %d coefficient(s) aliased by exact ",
          "collinearity and reported as omitted with zero estimate ",
          "and variance."
        ),
        sum(is.na(b_star))
      )
      if (!quiet) message(diagnostics$aliased)
    }
    b_star[is.na(b_star)] <- 0; names(b_star) <- c_names
    design_df <- if (ordinary_resampling) {
      if (vce_method == "bootstrap") Inf else max(1, N_obs - 1)
    } else {
      max(1, n_psu - n_strata)
    }

    calc_rep <- function(w, data) {
      if (!quiet) {
        env_pb$count <- env_pb$count + 1
        utils::setTxtProgressBar(
          pb,
          min(env_pb$count, n_replicates)
        )
      }
      na_ret <- stats::setNames(
        rep(NA_real_, ncol(X_full)),
        c_names
      )
      if (sum(w) == 0) return(na_ret)

      # tryCatch: calc_target can now stop() (e.g. an unstable Wagstaff denominator
      # in a given replicate); a failed replicate must return NA, not abort the run.
      tryCatch({
        precalc_vec_rep <- if (index_type == "rank" && !is.null(precalc_rank_var)) {
          data[[precalc_rank_var]]
        } else if (index_type == "level" && !is.null(precalc_level_var)) {
          data[[precalc_level_var]]
        } else {
          NULL
        }

        tgt_r <- calc_target(
          data[[dep_var]],
          if(!is.null(active_ses_var)) data[[active_ses_var]] else rep(1, nrow(data)),
          w, index_type,
          precalc_vec_rep,
          dep_min, dep_max, correction
        )

        idx_v <- w > 0
        m_r <- suppressWarnings(stats::lm.wfit(x = X_full[idx_v, , drop = FALSE], y = tgt_r[idx_v], w = w[idx_v]))
        cfs <- m_r$coefficients; cfs[is.na(cfs)] <- 0
        stats::setNames(cfs, c_names)
      }, error = function(e) na_ret)
    }

    if (ordinary_resampling && vce_method == "bootstrap") {
      bootstrap_audit <- new.env(parent = emptyenv())
      bootstrap_audit$calls <- 0L
      bootstrap_audit$weights <- matrix(
        0,
        nrow = N_obs,
        ncol = n_replicates
      )
      bootstrap_result <- boot::boot(
        data = seq_len(N_obs),
        statistic = function(data, indices) {
          bootstrap_audit$calls <- bootstrap_audit$calls + 1L
          frequencies <- tabulate(indices, nbins = N_obs)
          if (bootstrap_audit$calls > 1L) {
            replicate_number <- bootstrap_audit$calls - 1L
            if (replicate_number <= n_replicates) {
              bootstrap_audit$weights[, replicate_number] <-
                frequencies
            }
          }
          calc_rep(frequencies, df_clean)
        },
        R = n_replicates,
        sim = "ordinary",
        stype = "i",
        simple = TRUE
      )
      res_matrix <- list(
        theta = bootstrap_result$t0,
        replicates = bootstrap_result$t
      )
      replicate_weights_audit <- bootstrap_audit$weights
    } else if (ordinary_resampling) {
      row_ids <- seq_len(N_obs)
      jackknife_result <- resample::jackknife(
        data = row_ids,
        statistic = function(retained_rows) {
          if (length(retained_rows) == N_obs) {
            jackknife_weights <- rep(1, N_obs)
          } else {
            jackknife_weights <- rep(
              N_obs / (N_obs - 1),
              N_obs
            )
            jackknife_weights[
              setdiff(row_ids, retained_rows)
            ] <- 0
          }
          calc_rep(jackknife_weights, df_clean)
        },
        statisticNames = c_names,
        trace = FALSE
      )
      res_matrix <- list(
        theta = jackknife_result$observed,
        replicates = jackknife_result$replicates
      )
      replicate_weights_audit <- matrix(
        N_obs / (N_obs - 1),
        nrow = N_obs,
        ncol = N_obs
      )
      diag(replicate_weights_audit) <- 0
    } else {
      res_matrix <- suppressWarnings(
        survey::withReplicates(
          des_rep,
          calc_rep,
          return.replicates = TRUE
        )
      )
      replicate_weights_audit <- tryCatch(
        stats::weights(des_rep, type = "analysis"),
        error = function(e) NULL
      )
    }
    colnames(res_matrix$replicates) <- c_names
    if (!quiet) close(pb)
    valid_idx <- stats::complete.cases(res_matrix$replicates)
    fail_rate <- sum(!valid_idx) / length(valid_idx)
    if (
      !ordinary_resampling &&
        any(!valid_idx) &&
        !relax
    ) {
      stop(
        "CRITICAL: ",
        sum(!valid_idx),
        " survey replicate(s) failed. Dropping a survey replicate can ",
        "invalidate the design variance. Inspect the sparse or singular ",
        "replicates, or rerun with relax = TRUE for an explicitly ",
        "exploratory approximation.",
        call. = FALSE
      )
    }
    if (fail_rate > 0.05) {
      diagnostics$resampling <- sprintf("Statistical Alert: %.1f%% of resampling iterations failed. SEs might be less reliable.", fail_rate * 100)
      if (!relax && fail_rate > 0.2) stop(diagnostics$resampling) else if (!quiet) warning(diagnostics$resampling, call. = FALSE)
    }
    if (sum(valid_idx) < 2) stop("CRITICAL: Variance estimation failed.")
    reps <- res_matrix$replicates[valid_idx, , drop=FALSE]
    rscales <- if (
      length(generator_rscales) == length(valid_idx)
    ) {
      generator_rscales[valid_idx]
      } else {
        generator_rscales
      }
    valid_center_groups <- if (
      !is.null(replicate_center_groups)
    ) {
      replicate_center_groups[valid_idx]
    } else {
      NULL
    }
    if (ordinary_resampling) {
      valid_count <- sum(valid_idx)
      effective_scale <- if (vce_method == "bootstrap") {
        1 / (valid_count - 1)
      } else {
        (valid_count - 1) / valid_count
      }
      replicate_center <- colMeans(reps)
      center_method <- "replicate mean"
      design_df <- if (vce_method == "bootstrap") {
        Inf
      } else {
        max(1, valid_count - 1)
      }
      rep_adjust <- effective_scale / generator_scale
    } else {
      rep_adjust <- length(valid_idx) / sum(valid_idx)
      effective_scale <- survey_effective_scale * rep_adjust
      if (is.null(valid_center_groups)) {
        replicate_center <- colMeans(reps)
        center_method <- "replicate mean"
      } else {
        group_centers <- lapply(
          unique(valid_center_groups),
          function(group) {
            colMeans(
              reps[
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
        colnames(replicate_center) <- c_names
        center_method <-
          "stratum-specific replicate mean"
      }
    }
    if (any(!valid_idx)) {
      diagnostics$rep_adjust <- sprintf(
        paste0(
          "Variance Note: %d failed replicate(s) dropped; ",
          "the effective replicate scale is %.8g."
        ),
        sum(!valid_idx),
        effective_scale
      )
    }
    diffs_all <- if (is.matrix(replicate_center)) {
      reps - replicate_center
    } else {
      sweep(reps, 2, replicate_center, "-")
    }
    V_all <- effective_scale *
      crossprod(diffs_all * sqrt(rscales))
    rownames(V_all) <- colnames(V_all) <- c_names
    V_full <- V_all
    rownames(V_full) <- colnames(V_full) <- c_names
    se_index <- NA_real_
    replication_audit <- list(
      method = vce_method,
      design_type = generator_type,
      generator_type = generator_type,
      ordinary_resampling = ordinary_resampling,
      ordinary_bootstrap =
        ordinary_resampling && vce_method == "bootstrap",
      ordinary_jackknife =
        ordinary_resampling && vce_method == "jackknife",
      engine = if (
        ordinary_resampling && vce_method == "bootstrap"
      ) {
        "boot::boot"
      } else if (ordinary_resampling) {
        "resample::jackknife"
      } else {
        "survey::withReplicates"
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
      estimates = res_matrix$replicates,
      valid_estimates = reps,
      valid = valid_idx,
      requested_replicates = length(valid_idx),
      valid_replicates = sum(valid_idx),
      failed_replicates = sum(!valid_idx),
      replicate_factors = replicate_factors_audit,
      replicate_weights = replicate_weights_audit,
      generator_scale = generator_scale,
      effective_scale = effective_scale,
      scale = effective_scale,
      rscales = generator_rscales,
      used_rscales = rscales,
      failed_replicate_adjustment = rep_adjust,
      variance_rescale_after_failures = rep_adjust,
      center = if (is.matrix(replicate_center)) {
        replicate_center
      } else {
        stats::setNames(replicate_center, c_names)
      },
      center_groups = valid_center_groups,
      center_method = center_method,
      mse = FALSE,
      full_sample = b_star,
      vcov = V_all,
      covariance = V_full,
      index_variance = NULL,
      pending_methodological_review = FALSE
    )
  }

  # ==============================================================================
  # INFERENTIAL METRICS & WALD TESTS (LOGWORTH)
  # ==============================================================================
  se_flat <- safe_se(diag(V_full))
  stat_flat <- b_star / ifelse(is.na(se_flat) | se_flat == 0, NA_real_, se_flat)
  prob_t <- 2 * stats::pt(
    abs(stat_flat),
    df = design_df,
    lower.tail = FALSE
  )

  tb_anova <- tibble::tibble(Variable = character(), DF = integer(), F_Statistic = numeric(), Prob_F = numeric(), Logworth = numeric())
  prob_F_flat <- rep(NA, length(b_star)); names(prob_F_flat) <- c_names
  logworth_flat <- rep(NA, length(b_star)); names(logworth_flat) <- c_names

  for (var_name in names(auto_groups)) {
    g_cols <- c_names[auto_groups[[var_name]]]
    # which() guards against NA diagonal entries (a plain logical subset would
    # inject NA elements into valid_cols and crash the V_full indexing below).
    dv <- diag(V_full)[g_cols]
    valid_cols <- g_cols[which(!is.na(dv) & dv > 0)]
    if (length(valid_cols) == 0) next

    beta_g <- b_star[valid_cols]
    V_g <- V_full[valid_cols, valid_cols, drop=FALSE]
    V_inv <- tryCatch(solve(V_g), error = function(e) tryCatch(MASS::ginv(V_g), error = function(e2) NULL))

    if (!is.null(V_inv)) {
      W_stat <- as.numeric(t(beta_g) %*% V_inv %*% beta_g) / length(valid_cols)
      if (W_stat < 0) W_stat <- 0
      log_p_f <- stats::pf(
        W_stat,
        length(valid_cols),
        design_df,
        lower.tail = FALSE,
        log.p = TRUE
      )
      p_f <- exp(log_p_f)
      lw <- -log_p_f / log(10)

      prob_F_flat[valid_cols] <- p_f
      logworth_flat[valid_cols] <- lw
      tb_anova <- rbind(tb_anova, data.frame(Variable = var_name, DF = length(valid_cols), F_Statistic = W_stat, Prob_F = p_f, Logworth = lw))
    }
  }
  tb_anova <- tb_anova %>%
  dplyr::arrange(dplyr::desc(.data$Logworth)) %>%
  tibble::as_tibble()

  # ==============================================================================
  # RESULTS ASSEMBLY & OUTPUT (Bug Fix 5a)
  # ==============================================================================
  q_val <- 1 - ((1 - level) / 2)
  display_terms <- unname(sapply(c_names, function(x) ifelse(x %in% names(display_dict), display_dict[[x]], x)))

  # Shared inferential columns (Term, Estimate, Std_Error, Statistic, P_Value,
  # Conf_Low, Conf_High) align with oby_decomp/f_decomp/wvw_decomp. 'Estimate' is
  # the marginal effect (NOT a contribution, per Kessels & Erreygers 2019). The
  # KE-specific F-test p-value (Prob_F) and Logworth importance are preserved.
  tb_detailed <- tibble::tibble(
    Term = display_terms,
    Estimate = as.numeric(b_star),
    Std_Error = as.numeric(se_flat),
    Statistic = as.numeric(stat_flat),
    P_Value = as.numeric(prob_t),
    Prob_F = as.numeric(prob_F_flat),
    Logworth = as.numeric(logworth_flat),
    Omitted = as.logical(aliased_coefficients),
    Conf_Low = as.numeric(
      b_star - stats::qt(q_val, df = design_df) * se_flat
    ),
    Conf_High = as.numeric(
      b_star + stats::qt(q_val, df = design_df) * se_flat
    )
  )

  Metric_Name <- dplyr::case_when(
    index_type == "rank" & correction == "erreygers" ~ "Erreygers' Rank Index (Bounded Absolute)",
    index_type == "rank" & correction == "wagstaff" ~ "Wagstaff's Rank Index (Bounded Relative)",
    index_type == "rank" & correction == "generalized" ~ "Generalized Concentration Index (Unbounded Absolute)",
    index_type == "rank" & correction == "standard" ~ "Standard Concentration Index (Unbounded Relative)",
    index_type == "level" & correction == "erreygers" ~ "Erreygers' Level Index (Bounded Absolute)",
    index_type == "level" & correction == "generalized" ~ "Generalized Level Index (Unbounded Absolute)",
    index_type == "level" & correction == "standard" ~ "Standard Level Index (Unbounded Relative)"
  )

  # Kessels and Erreygers report the index point but do not report its standard
  # error, test, or confidence interval. Keep the harmonized columns as explicit
  # NA values; all VCE extensions in this function target regression
  # coefficients only.
  ov_stat <- NA_real_
  tb_overall <- tibble::tibble(
    Term = Metric_Name,
    Estimate = as.numeric(Index_Value),
    Std_Error = as.numeric(se_index),
    Statistic = as.numeric(ov_stat),
    P_Value = NA_real_,
    Conf_Low = NA_real_,
    Conf_High = NA_real_,
    Min_Health = dep_min,
    Max_Health = dep_max
  )

  SS_tot <- sum(df_clean$svy_weight * (df_clean$TARGET_Y - Index_Value)^2)
  res_fit <- df_clean$TARGET_Y - as.vector(X_full %*% b_star)
  target_spread <- diff(range(df_clean$TARGET_Y))
  target_spread_tolerance <- 100 * .Machine$double.eps *
    max(1, max(abs(df_clean$TARGET_Y)))

  Model_Label <- if (use_svy) "Survey-weighted direct regression" else "Direct regression (OLS/WLS)"
  fit_metrics <- tibble::tibble(
    Model = Model_Label,
    R_Squared = if (
      is.finite(SS_tot) &&
        target_spread > target_spread_tolerance
    ) {
      1 - (sum(df_clean$svy_weight * res_fit^2) / SS_tot)
    } else {
      NA_real_
    }
  )

  source_rows_sorted <- df_clean[[source_row_name]]
  complete_case_rows <- selected_data[[source_row_name]][complete_case_mask]
  model_matrix_audit <- X_full
  rownames(model_matrix_audit) <- as.character(source_rows_sorted)
  fitted_values <- as.vector(X_full %*% b_star)
  names(fitted_values) <- as.character(source_rows_sorted)
  residual_values <- as.numeric(df_clean$TARGET_Y - fitted_values)
  names(residual_values) <- as.character(source_rows_sorted)
  coefficient_mean_identity <- sum(
    b_star * colSums(
      X_full * target_components$normalized_weight
    )
  )
  overall_vcov <- matrix(
    se_index^2,
    nrow = 1,
    ncol = 1,
    dimnames = list("Index", "Index")
  )
  target_audit <- data.frame(
    source_row = source_rows_sorted,
    sort_value = as.numeric(df_clean[[sort_var]]),
    health = as.numeric(df_clean[[dep_var]]),
    weight = as.numeric(df_clean$svy_weight),
    normalized_weight = target_components$normalized_weight,
    position = target_components$position,
    fractional_rank = if (index_type == "rank") {
      target_components$position
    } else {
      rep(NA_real_, N_obs)
    },
    relative_level = if (index_type == "level") {
      target_components$position
    } else {
      rep(NA_real_, N_obs)
    },
    target_base = target_components$target_base,
    correction_scale = rep(
      target_components$correction_scale,
      N_obs
    ),
    target = target_components$target,
    fitted = fitted_values,
    residual = residual_values,
    stringsAsFactors = FALSE
  )
  analytic_data_audit <- as.data.frame(
    df_clean[, c(source_row_name, req_vars), drop = FALSE]
  )
  names(analytic_data_audit)[
    names(analytic_data_audit) == source_row_name
  ] <- "source_row"

    raw <- list(
    schema_version = "1.0",
    settings = list(
      dep_var = dep_var,
      indep_vars = indep_vars,
      ses_var = ses_var,
      index_type = index_type,
      correction = correction,
      precalc_rank_var = precalc_rank_var,
      precalc_level_var = precalc_level_var,
      dep_min = dep_min,
      dep_max = dep_max,
      use_svy = use_svy,
      weight_var = weight_var,
      strata_var = strata_var,
      psu_var = psu_var,
      vce_method = vce_method,
      boot_reps = boot_reps,
      lonely_psu = lonely_psu,
      level = level,
      seed = seed,
      relax = relax,
      numeric_rounding = "none"
    ),
    sample = list(
      input_n = input_n,
      analytic_n = N_obs,
      analytic_source_rows = source_rows_sorted,
      complete_case_source_rows_before_sort = complete_case_rows,
      excluded_source_rows = list(
        missing_required_values = missing_excluded_rows
      ),
      sort_variable = sort_var,
      analytic_data = analytic_data_audit
    ),
      design = list(
      point_weight_semantics = if (!is.null(weight_var)) {
        paste(
          "sampling/design weight; normalized consistently for",
          "position, target mean, index, covariate means, and WLS"
        )
      } else {
        "equal observation weight"
      },
      point_estimation_equation = paste(
        "weighted least squares on the unweighted individual KE target;",
        "equivalent to integer frequency expansion at the point estimate"
      ),
      weights = stats::setNames(
        as.numeric(df_clean$svy_weight),
        as.character(source_rows_sorted)
      ),
      normalized_weights = stats::setNames(
        target_components$normalized_weight,
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
    position = list(
      source = if (!is.null(active_precalc_var)) {
        paste0("precalculated variable: ", active_precalc_var)
      } else {
        paste0("recomputed from socioeconomic variable: ", active_ses_var)
      },
      precalculated = !is.null(active_precalc_var),
      replicate_semantics = if (!is.null(active_precalc_var)) {
        paste(
          "user-supplied positions remain fixed within all jackknife and",
          "bootstrap replicates; sampled records retain their supplied value"
        )
      } else if (vce_method %in% c("jackknife", "bootstrap")) {
        paste(
          "positions are recomputed from the socioeconomic variable and",
          "the replicate weights in every replicate"
        )
      } else {
        "positions are computed from the socioeconomic variable"
      }
    ),
    target = target_audit,
    model = list(
      formula = form_mod,
      matrix = model_matrix_audit,
      assign = assign_attr,
      term_labels = term_labels,
      term_groups = auto_groups,
      contrasts = attr(X_full, "contrasts"),
      factor_levels = lapply(
        df_clean[intersect(model_variables, names(df_clean))],
        function(variable) {
          if (is.factor(variable)) levels(variable) else NULL
        }
      ),
      coefficients = b_star,
      aliased = aliased_coefficients,
      vcov = V_full,
      vcov_type = model_vcov_type,
      fitted = fitted_values,
      residuals = residual_values,
      object_class = class(model_object)
    ),
    index = list(
      health_mean = target_components$health_mean,
      socioeconomic_mean = target_components$socioeconomic_mean,
      correction_scale = target_components$correction_scale,
      target_mean = Index_Value,
      regression_at_covariate_means = coefficient_mean_identity,
      identity_error = coefficient_mean_identity - Index_Value,
      weighted_residual_mean = sum(
        target_components$normalized_weight * residual_values
      )
    ),
    estimates = list(
      overall = stats::setNames(Index_Value, Metric_Name),
      detailed = stats::setNames(as.numeric(b_star), c_names),
      anova = tb_anova
    ),
    vcov = list(
      model = V_full,
      overall = overall_vcov,
      scope = if (vce_method == "linearized") {
        "conditional_on_empirically_constructed_target"
      } else {
        "coefficient_covariance_with_complete_target_recomputation"
      },
      includes_generated_target_uncertainty =
        vce_method != "linearized",
      model_complete_for_constructed_target = TRUE,
      model_complete_for_full_estimator =
        vce_method != "linearized",
      coefficient_covariance_includes_target_reconstruction =
        vce_method != "linearized",
      overall_complete = FALSE,
      joint_index_coefficient_covariance = NULL,
      recommended_complete_methods =
        if (vce_method == "linearized") {
          c("jackknife", "bootstrap")
        } else {
          character()
        }
    ),
    replication = replication_audit,
    known_limitations = list(
      linearized = if (vce_method == "linearized") {
        if (use_svy) {
          paste(
            "The survey Taylor covariance conditions on the empirically",
            "constructed target; uncertainty from ranks, means, and",
            "correction scaling is not propagated. This conditional design",
            "covariance is a package extension and is not specified by",
            "Kessels and Erreygers (2019)."
          )
        } else {
          paste(
            "The model covariance conditions on the empirically constructed",
            "target; uncertainty from ranks, means, and correction scaling is",
            "not propagated. This HC1 covariance is a package extension and is",
            "not specified by Kessels and Erreygers (2019)."
          )
        }
      } else {
        NULL
      },
      replication = NULL,
      overall_inference = paste(
        "Kessels and Erreygers (2019) report the aggregate index as a",
        "point estimate without standard error, test, or confidence interval.",
        "Package covariance extensions intentionally target only the",
        "direct-regression coefficients."
      ),
      survey_design = if (use_svy) {
        paste(
          "The built-in design represents one final weight, optional strata,",
          "and one PSU stage. Lower stages, FPC, calibration, and external",
          "replicate weights are not represented."
        )
      } else {
        NULL
      }
    )
  )

  if (!quiet) {
    cat("\n", rep("=", 85), "\n DIRECT REGRESSION DECOMPOSITION\n", rep("=", 85), "\n", sep="")
    cat(sprintf(" Number of obs    = %-15d Index Base Type  = %s\n", N_obs, toupper(index_type)))
    cat(sprintf(" Pop. size        = %-15.2f Metric Corrected = %s\n", pop_size, toupper(correction)))
    cat(sprintf(" Design df        = %-15s VCE Engine       = %s\n", ifelse(is.infinite(design_df), "Inf", design_df), toupper(vce_method)))
    cat(sprintf(" Dep Variable     = %-15s Index Value      = %.5f\n", dep_var, Index_Value))
    cat(rep("-", 85), "\n", sep="")

    if (length(diagnostics) > 0) { cat("\n [DIAGNOSTICS & WARNINGS]\n"); for (warn in diagnostics) cat(" *", warn, "\n"); cat(rep("-", 85), "\n", sep="") }
    tb_overall_print <- tb_overall
    tb_overall_print$Estimate <- round(tb_overall_print$Estimate, 6)
    tb_overall_print$Std_Error <- round(tb_overall_print$Std_Error, 6)
    tb_overall_print$Statistic <- round(tb_overall_print$Statistic, 2)
    tb_overall_print$P_Value <- round(tb_overall_print$P_Value, 4)
    tb_overall_print$Conf_Low <- round(tb_overall_print$Conf_Low, 5)
    tb_overall_print$Conf_High <- round(tb_overall_print$Conf_High, 5)
    fit_metrics_print <- fit_metrics
    fit_metrics_print$R_Squared <- round(
      fit_metrics_print$R_Squared,
      4
    )
    tb_detailed_print <- tb_detailed
    tb_detailed_print$Estimate <- round(tb_detailed_print$Estimate, 6)
    tb_detailed_print$Std_Error <- round(tb_detailed_print$Std_Error, 6)
    tb_detailed_print$Statistic <- round(tb_detailed_print$Statistic, 2)
    tb_detailed_print$P_Value <- round(tb_detailed_print$P_Value, 4)
    tb_detailed_print$Prob_F <- round(tb_detailed_print$Prob_F, 4)
    tb_detailed_print$Logworth <- round(tb_detailed_print$Logworth, 3)
    tb_detailed_print$Conf_Low <- round(tb_detailed_print$Conf_Low, 5)
    tb_detailed_print$Conf_High <- round(tb_detailed_print$Conf_High, 5)

    cat("\n [OVERALL INEQUALITY]\n"); print(as.data.frame(tb_overall_print), row.names = FALSE)
    cat("\n [MODEL FIT METRICS]\n"); print(as.data.frame(fit_metrics_print), row.names = FALSE)
    cat("\n [ANOVA: IMPORTANCE RANKING BY VARIABLE (LOGWORTH)]\n"); print(as.data.frame(tb_anova %>% dplyr::mutate_if(is.numeric, ~round(., 4))), row.names = FALSE)
    cat("\n [DETAILED MARGINAL EFFECTS]\n"); print(as.data.frame(tb_detailed_print[, c("Term", "Estimate", "Std_Error", "Statistic", "P_Value", "Prob_F", "Logworth")]), row.names = FALSE)
    cat(rep("-", 85), "\n\n", sep="")
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
    results_anova = tibble::as_tibble(tb_anova),
    results_detailed = tb_detailed,
    raw = raw
  )))
}
