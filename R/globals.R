# CRAN checks may not see non-standard evaluation variables used in
# model formulas and dplyr expressions. These declarations are only for
# static code checking and do not affect package calculations.

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".lunadecomp_weight"
  ))
}
