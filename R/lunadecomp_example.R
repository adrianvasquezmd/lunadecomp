#' Simulated Example Data for Decomposition Methods
#'
#' @description
#' A simulated individual-level dataset for examples and tests in
#' the `lunadecomp` package. The dataset includes a two-group comparison
#' variable, socioeconomic ranking variables, survey design variables, continuous
#' and binary outcomes, and common health inequality covariates.
#'
#' @format A data frame with 1,000 rows and 17 variables:
#' \describe{
#'   \item{id}{Individual identifier.}
#'   \item{group}{Binary comparison group. The examples treat 0 as favored and
#'   1 as disadvantaged.}
#'   \item{ses}{Continuous socioeconomic status score, with higher values
#'   indicating higher status. Values may be negative.}
#'   \item{weight}{Positive simulated sampling weight.}
#'   \item{strata}{Simulated survey strata.}
#'   \item{psu}{Simulated primary sampling unit.}
#'   \item{age}{Age in years.}
#'   \item{female}{Sex indicator, factor with levels `male` and `female`.}
#'   \item{rural}{Residence, factor with levels `urban` and `rural`.}
#'   \item{insurance}{Health insurance, factor with levels `uninsured` and
#'   `insured`.}
#'   \item{sanitation}{Sanitation, factor with levels `unimproved` and
#'   `improved`.}
#'   \item{education}{Educational attainment, factor with levels `primary`,
#'   `secondary`, and `higher`.}
#'   \item{region}{Geographic region, factor with levels `coast`, `highlands`,
#'   `jungle`, and `capital`.}
#'   \item{rank}{Unweighted fractional position computed as
#'   `rank(ses) / (nrow(lunadecomp_example) + 1)`. This supplied position differs
#'   from ranks recomputed by a decomposition using its analytical sample and
#'   weights.}
#'   \item{y_continuous}{Simulated continuous health-related outcome.}
#'   \item{y_binary}{Simulated binary health-related outcome.}
#'   \item{health_score}{Simulated bounded health score ranging from 0 to 1.}
#' }
#'
#' @details
#' The outcomes, covariates, and survey design are simulated; the data do not
#' represent a population survey. For domain examples, pass the full dataset
#' and select analytical rows with `analysis_subset`.
#'
#' @source Simulated data generated for package examples and unit tests.
#' @seealso [lunadecomp-package]
#'
#' @examples
#' data(lunadecomp_example)
#' head(lunadecomp_example)
#' table(lunadecomp_example$group)
#'
"lunadecomp_example"
