# Data generation script for lunadecomp_example
# This script creates a simulated dataset used in examples and tests.

set.seed(20260531)

n <- 1000

# Internal numeric indicators used only to generate outcomes.
# The exported dataset stores the binary covariates as factors so that
# package examples exercise categorical-variable handling.
group_num <- rbinom(n, size = 1, prob = 0.45)
female_num <- rbinom(n, size = 1, prob = 0.52)
rural_num <- rbinom(n, size = 1, prob = 0.35)
insurance_num <- rbinom(n, size = 1, prob = 0.70)
sanitation_num <- rbinom(n, size = 1, prob = 0.68)

lunadecomp_example <- data.frame(
  id = seq_len(n),

  # Two-group comparison variable. Keep numeric 0/1 because decomposition
  # functions expect the comparison group to be coded as 0 and 1.
  group = group_num,

  # Socioeconomic ranking variable
  ses = rnorm(n, mean = 0, sd = 1),

  # Survey design variables
  weight = runif(n, min = 0.5, max = 2.5),
  strata = sample(paste0("S", 1:8), n, replace = TRUE),
  psu = sample(paste0("PSU", 1:120), n, replace = TRUE),

  # Covariates
  age = round(rnorm(n, mean = 42, sd = 13)),
  female = factor(
    ifelse(female_num == 1, "female", "male"),
    levels = c("male", "female")
  ),
  rural = factor(
    ifelse(rural_num == 1, "rural", "urban"),
    levels = c("urban", "rural")
  ),
  insurance = factor(
    ifelse(insurance_num == 1, "insured", "uninsured"),
    levels = c("uninsured", "insured")
  ),
  sanitation = factor(
    ifelse(sanitation_num == 1, "improved", "unimproved"),
    levels = c("unimproved", "improved")
  ),
  education = factor(
    sample(c("primary", "secondary", "higher"), n, replace = TRUE,
           prob = c(0.35, 0.45, 0.20)),
    levels = c("primary", "secondary", "higher")
  ),
  region = factor(
    sample(c("coast", "highlands", "jungle", "capital"), n, replace = TRUE),
    levels = c("coast", "highlands", "jungle", "capital")
  )
)

# Fractional socioeconomic rank, useful for concentration-index methods
lunadecomp_example$rank <- rank(
  lunadecomp_example$ses,
  ties.method = "average"
) / (n + 1)

# Continuous outcome
lunadecomp_example$y_continuous <-
  50 +
  3.0 * lunadecomp_example$ses -
  4.0 * group_num +
  0.20 * lunadecomp_example$age -
  2.5 * rural_num +
  3.5 * insurance_num +
  2.0 * (lunadecomp_example$education == "secondary") +
  5.0 * (lunadecomp_example$education == "higher") +
  rnorm(n, mean = 0, sd = 8)

# Binary outcome probability
eta <-
  -0.8 +
  0.50 * lunadecomp_example$ses -
  0.45 * group_num +
  0.015 * lunadecomp_example$age -
  0.50 * rural_num +
  0.45 * insurance_num +
  0.35 * sanitation_num +
  0.25 * (lunadecomp_example$education == "secondary") +
  0.55 * (lunadecomp_example$education == "higher")

p <- stats::plogis(eta)

lunadecomp_example$y_binary <- rbinom(n, size = 1, prob = p)

# Bounded health outcome, useful for concentration-index decompositions
lunadecomp_example$health_score <- pmin(
  pmax(
    0,
    0.50 +
      0.12 * lunadecomp_example$ses -
      0.08 * group_num -
      0.10 * rural_num +
      0.08 * insurance_num +
      rnorm(n, mean = 0, sd = 0.15)
  ),
  1
)

usethis::use_data(
  lunadecomp_example,
  overwrite = TRUE,
  compress = "xz"
)
