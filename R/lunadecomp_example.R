#' Simulated Example Data for Decomposition Methods
#'
#' @description
#' A simulated individual-level dataset for examples, tests, and vignettes in
#' the `lunadecomp` package. The dataset includes a two-group comparison
#' variable, socioeconomic ranking variables, survey design variables, continuous
#' and binary outcomes, and common health inequality covariates.
#'
#' @format A data frame with 1,000 rows and 17 variables:
#' \describe{
#'   \item{id}{Individual identifier.}
#'   \item{group}{Binary comparison group coded as 0 and 1.}
#'   \item{ses}{Continuous socioeconomic status variable.}
#'   \item{weight}{Sampling weight.}
#'   \item{strata}{Simulated survey strata.}
#'   \item{psu}{Simulated primary sampling unit.}
#'   \item{age}{Age in years.}
#'   \item{female}{Sex indicator, factor with levels `male` and `female`.}
#'   \item{rural}{Binary rural residence indicator.}
#'   \item{insurance}{Binary health insurance indicator.}
#'   \item{sanitation}{Binary improved sanitation indicator.}
#'   \item{education}{Educational attainment: primary, secondary, or higher.}
#'   \item{region}{Geographic region.}
#'   \item{rank}{Fractional socioeconomic rank based on `ses`.}
#'   \item{y_continuous}{Simulated continuous health-related outcome.}
#'   \item{y_binary}{Simulated binary health-related outcome.}
#'   \item{health_score}{Simulated bounded health score ranging from 0 to 1.}
#' }
#'
#' @source Simulated data generated for package examples and unit tests.
#'
#' @examples
#' data(lunadecomp_example)
#' head(lunadecomp_example)
#' table(lunadecomp_example$group)
#'
"lunadecomp_example"