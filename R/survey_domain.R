# Validate the optional survey-domain mask.
.ld_validate_analysis_subset <- function(data, analysis_subset, use_svy) {
  if (is.null(analysis_subset)) return(NULL)
  if (!isTRUE(use_svy)) {
    stop("Error: analysis_subset requires use_svy = TRUE.", call. = FALSE)
  }
  if (!is.logical(analysis_subset) || !is.null(dim(analysis_subset)) ||
      length(analysis_subset) != nrow(data) || anyNA(analysis_subset)) {
    stop("Error: analysis_subset must be a logical vector of length nrow(data), without missing values.", call. = FALSE)
  }
  if (!any(analysis_subset)) {
    stop("Error: analysis_subset selects no observations.", call. = FALSE)
  }
  unname(analysis_subset)
}

# Build the parent design independently of model missingness. Preserve the
# caller's analytical row order, including SES sorting.
.ld_survey_domain_frame <- function(data, analytic_rows, weight_var,
                                    strata_var, psu_var) {
  n <- nrow(data)
  design_vars <- unique(c(weight_var, strata_var, psu_var))
  valid_design <- if (length(design_vars)) {
    stats::complete.cases(data[, design_vars, drop = FALSE])
  } else rep(TRUE, n)
  w <- if (is.null(weight_var)) rep(1, n) else {
    values <- data[[weight_var]]
    if (is.numeric(values)) as.numeric(values) else
      suppressWarnings(as.numeric(as.character(values)))
  }
  if (any(valid_design & (!is.finite(w) | w < 0))) {
    stop("Error: survey-domain weights must be finite and non-negative, including outside the analytical domain.", call. = FALSE)
  }
  valid_design <- valid_design & !is.na(w) & w > 0
  strata <- if (is.null(strata_var)) rep("1", n) else as.character(data[[strata_var]])
  psu <- if (is.null(psu_var)) as.character(seq_len(n)) else as.character(data[[psu_var]])
  if (!length(analytic_rows) || anyNA(analytic_rows) ||
      anyDuplicated(analytic_rows) ||
      any(!analytic_rows %in% which(valid_design))) {
    stop("Error: analytical rows are not a non-empty subset of the valid survey design.", call. = FALSE)
  }
  eligible_rows <- which(valid_design)
  participating_strata <- unique(strata[analytic_rows])
  design_rows <- which(valid_design & strata %in% participating_strata)
  frame <- data.frame(
    svy_weight = w[design_rows],
    svy_strata = strata[design_rows],
    svy_psu = psu[design_rows],
    row.names = as.character(design_rows),
    stringsAsFactors = FALSE
  )
  psus <- unique(frame[c("svy_strata", "svy_psu")])
  n_strata <- length(unique(frame$svy_strata))
  list(
    frame = frame,
    analytic_positions = match(analytic_rows, design_rows),
    metadata = list(
      selection = "explicit analysis_subset and method-specific complete cases",
      eligible_design_rows = eligible_rows,
      design_rows = design_rows,
      analytic_rows = as.integer(analytic_rows),
      outside_domain_rows = setdiff(design_rows, analytic_rows),
      excluded_design_rows = setdiff(seq_len(n), eligible_rows),
      excluded_strata = setdiff(unique(strata[eligible_rows]), participating_strata),
      design_psus = nrow(psus),
      design_strata = n_strata,
      design_df = nrow(psus) - n_strata,
      domain_psus = nrow(unique(data.frame(
        stratum = strata[analytic_rows], psu = psu[analytic_rows]
      ))),
      replicate_weights_row_order = "analytical source rows, after method sorting"
    )
  )
}
