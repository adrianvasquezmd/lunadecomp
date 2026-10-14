# lunadecomp

<!-- badges: start -->
[![R-CMD-check](https://github.com/adrianvasquezmd/lunadecomp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adrianvasquezmd/lunadecomp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`lunadecomp` provides four methods for decomposing two-group outcome gaps
or socioeconomic health inequality indices. This documentation describes
version **0.1.1**.

| Function | Target | Supported outcome model |
|---|---|---|
| `oby_decomp()` | Oaxaca-Blinder / Yun explained and unexplained group gaps | OLS, linear probability, logit, probit |
| `f_decomp()` | Fairlie binary-outcome gap and detailed explained effects | Logit, probit |
| `wvw_decomp()` | Determinant contributions to a concentration-type index | OLS, logit, probit |
| `ke_decomp()` | Direct-regression marginal effects and variable importance | Regression of a rank- or level-dependent composite target |

## Installation

Install version 0.1.1 from GitHub:

```r
install.packages("remotes")
remotes::install_github("adrianvasquezmd/lunadecomp@v0.1.1")
```

Load the package and its simulated example data:

```r
library(lunadecomp)
data(lunadecomp_example)
d <- lunadecomp_example
predictors <- c("age", "rural", "insurance", "sanitation", "education")
```

The dataset has 1,000 observations, 17 variables, continuous and binary
outcomes, a bounded health score, categorical predictors, and simulated
weights, strata, and PSUs. It is intended for demonstrating the methods.
Use `?lunadecomp` for the package overview and `?lunadecomp_example` for
variable definitions.

## Oaxaca-Blinder and Yun decomposition

Declare the comparison order explicitly: the first group is the favored
group and the second is the disadvantaged group. The reported gap is
`E(Y | favored) - E(Y | disadvantaged)`, estimated from mean group-model
predictions.

```r
fit_ob <- oby_decomp(
  data = d,
  dep_var = "y_continuous",
  group_var = "group",
  group_levels = c(0, 1),
  indep_vars = predictors,
  model_type = "ols",
  ref_method = "pooled",
  normalize = c("insurance", "education"),
  vce_method = "linearized",
  quiet = TRUE
)

fit_ob$results_overall
fit_ob$results_detailed_explained
fit_ob$results_detailed_unexplained
```

Use `model_type = "lpm"`, `"logit"`, or `"probit"` for a binary outcome
coded 0/1. Numeric outcomes stored as text or factors are converted by their
values, not their factor codes.

`group_levels` can select two levels of a variable containing additional
groups. Other groups are excluded from the models; in survey mode their
valid design information is retained in participating strata. Alternatively,
`favored_group` identifies the favored level in data with exactly two groups.

The default `normalize = NULL` retains treatment coding. Named factors or
`normalize = "all"` use equal-category normalization and include the omitted
level in detailed results. Factor-by-continuous interactions are supported:

```r
fit_ob_interaction <- oby_decomp(
  data = d,
  dep_var = "y_continuous",
  group_var = "group",
  group_levels = c(0, 1),
  indep_vars = c("education * age", "rural"),
  normalize = "education",
  groupings = list(
    Education = c("education", "education:age"),
    Age = "age",
    Residence = "rural"
  ),
  quiet = TRUE
)
fit_ob_interaction$results_grouped_explained
```

For nonlinear models, `pooled_anchor = "favored"`, `"disadvantaged"`, or
`"centered"` changes the pooled counterfactual and can change the components.
It does not change the favored-minus-disadvantaged orientation of the gap.

## Fairlie nonlinear decomposition

```r
fit_fairlie <- f_decomp(
  data = d,
  dep_var = "y_binary",
  group_var = "group",
  group_levels = c(0, 1),
  indep_vars = predictors,
  model_type = "logit",
  ref_method = "pooled",
  vce_method = "linearized",
  reps = 25,
  seed = 2026,
  quiet = TRUE
)

fit_fairlie$results_overall
fit_fairlie$results_detailed_explained
fit_fairlie$raw$matching$precision
```

Fairlie reports observed group means and their difference in `results_overall`.
`results_detailed_unexplained` is `NULL`: Fairlie does not allocate the
unexplained component by predictor. Factor levels are exchanged jointly.
Formula interactions are rejected; a precomputed interaction is appropriate
only when all its main effects and interaction columns are exchanged in one
joint `groupings` block.

`reps` controls matching and randomized block-order iterations; `boot_reps`
controls sampling bootstrap replicates. The examples use 25 to keep runtime
short. The defaults are 1,000 and 500, respectively. Assess stability as these
counts increase; a small `reps` can affect both contributions and sampling
standard errors. `raw$matching$precision` reports algorithmic Monte Carlo
error separately from sampling uncertainty.

Detailed effects are not rescaled to the full-sample explained total. Their
sum can differ from that total under sampled matching. The linearized
detailed covariance has unestimated off-diagonal entries stored as zeros;
`raw$detailed_vcov_complete` is `FALSE`. Do not use that matrix for joint
tests. Jackknife and bootstrap return complete covariance matrices for the
replicated estimator.

Both gap methods support pooled, Neumark, Cotton, and Reimers references.
The favored-group reference is `"group1"` in `oby_decomp()` and `"group0"`
in `f_decomp()`; the disadvantaged-group reference is `"group2"` and
`"group1"`, respectively.

## Wagstaff, van Doorslaer, and Watanabe decomposition

```r
fit_wvw <- wvw_decomp(
  data = d,
  dep_var = "health_score",
  indep_vars = predictors,
  ses_var = "ses",
  model_type = "ols",
  correction = "erreygers",
  dep_min = 0,
  dep_max = 1,
  vce_method = "linearized",
  quiet = TRUE
)

fit_wvw$results_overall
fit_wvw$results_detailed
fit_wvw$results_grouped
```

In the detailed table, `Estimate` is the determinant's contribution. Other
columns include `Marginal Effects`, `Mean_X`, `Elasticity`, `Concentration
Index`, `% Contribution Explained`, and `% Contribution Total`.

Available corrections are `"standard"`, `"generalized"`, `"erreygers"`, and
`"wagstaff"`. Bounded corrections require theoretical `dep_min` and `dep_max`
for OLS; logit/probit use bounds 0 and 1. Observed sample extrema are not
substituted for theoretical bounds. `is_shortfall = TRUE` records an adverse
outcome orientation without recoding the outcome.

WVW accepts predictor column names, not formula transformations or
interactions. For nonlinear models, factors and binary indicators receive
discrete probability changes; continuous predictors receive average
derivatives of the inverse link. The overall
table includes a `Non-linear Approximation Error` row so the components
reconcile with the total index.

## Kessels-Erreygers direct regression

```r
fit_ke <- ke_decomp(
  data = d,
  dep_var = "health_score",
  indep_vars = predictors,
  ses_var = "ses",
  index_type = "rank",
  correction = "erreygers",
  dep_min = 0,
  dep_max = 1,
  vce_method = "linearized",
  quiet = TRUE
)

fit_ke$results_overall
fit_ke$results_anova
fit_ke$results_detailed
```

`Estimate` in the detailed table is a direct-regression marginal effect.
These effects are not additive contributions to the aggregate index.
`results_anova` reports variable-level F tests and logworth. Formula terms
and interactions, such as `indep_vars = c("age * education", "rural")`, are
supported.

The aggregate index is point-only under every variance method: its standard
error, test, and confidence interval are `NA`. Linearized coefficient
inference is conditional on the constructed target and does not propagate
uncertainty in its ranks, means, or correction scale. Replicate methods
reconstruct the target and refit the regression.

For a level-dependent index, use a meaningful socioeconomic level, such as
income, with `index_type = "level"`. Erreygers and Wagstaff corrections
require explicit health bounds; Wagstaff is available only for rank-dependent
indices. The default KE correction is `"erreygers"`.

## Survey designs and variance estimation

Declare designs explicitly with `use_svy = TRUE`, `weight_var`, `strata_var`,
and `psu_var`. A missing weight variable means unit weights; without a PSU
variable, each design-source row is its own PSU. PSU labels are nested within
strata.

A weight variable automatically activates survey mode in all four functions.
For Fairlie, WVW, and KE, a stratum or PSU variable also activates it. For
OBY, set `use_svy = TRUE` when declaring strata or PSUs without a weight.
Weights represent sampling weights, not frequency counts.

| Function | `linearized` | `jackknife` / `bootstrap` |
|---|---|---|
| OBY | Robust delta method or survey Taylor covariance, with the approximations documented in its help | Refits models and recomputes the decomposition |
| Fairlie | Available without declared strata or PSUs; detailed cross-covariances are incomplete | Required with declared strata or PSUs; re-runs matching within each sampling replicate |
| WVW | Joint covariance of contributions, total, residual, and nonlinear approximation error | Recomputes the model and decomposition |
| KE | Coefficient covariance conditional on the constructed target | Reconstructs the target and estimates coefficient covariance; no aggregate-index inference |

Ordinary replication resamples observations. Survey replication uses
PSU-based design weights, JK1/JKn jackknife, or Rao-Wu rescaled bootstrap.
Survey replicate inference uses `PSUs - strata` degrees of freedom. Linearized
weight-only conventions differ by method; consult the function help.

### Survey domains and a common analytical sample

Pass the full parent survey sample and select the domain with
`analysis_subset`. Filtering `data` first removes the excluded records from
the design. A mask must be logical, have length `nrow(data)`, contain no
missing values, and select at least one record.

This example selects participants aged 40 or older and a common complete-case
sample across the outcomes and predictors:

```r
analysis_vars <- c("y_continuous", "y_binary", "health_score", "group",
                   "ses", predictors)
analysis_domain <- complete.cases(d[, analysis_vars]) &
  !is.na(d$age) & d$age >= 40

fit_wvw_domain <- wvw_decomp(
  data = d,
  dep_var = "health_score",
  indep_vars = predictors,
  ses_var = "ses",
  correction = "erreygers",
  dep_min = 0,
  dep_max = 1,
  use_svy = TRUE,
  weight_var = "weight",
  strata_var = "strata",
  psu_var = "psu",
  analysis_subset = analysis_domain,
  vce_method = "linearized",
  quiet = TRUE
)

fit_wvw_domain$summary_stats
fit_wvw_domain$raw$survey_domain[c("design_psus", "domain_psus", "design_df")]
```

All four functions accept this mask. Model missingness and group selection
can further restrict it. Valid positive-weight design rows in participating
strata are retained even when their outcomes or covariates are missing;
PSUs with no analytical observations still contribute design information.
Entirely outside-domain strata are excluded. Empirical ranks and model
estimates use only analytical records. `analysis_subset = NULL` retains each
method's existing complete-case policy.

For Fairlie with this design, use a replicate variance estimator:

```r
fit_fairlie_survey <- f_decomp(
  data = d,
  dep_var = "y_binary",
  group_var = "group",
  group_levels = c(0, 1),
  indep_vars = predictors,
  use_svy = TRUE,
  weight_var = "weight",
  strata_var = "strata",
  psu_var = "psu",
  analysis_subset = analysis_domain,
  vce_method = "bootstrap",
  reps = 25,
  boot_reps = 25,
  seed = 2026,
  quiet = TRUE
)

fit_fairlie_survey$results_overall
fit_fairlie_survey$raw$matching$precision
fit_fairlie_survey$raw$replication$failed_replicates
```

### Design assumptions and diagnostics

`lonely_psu` controls Taylor covariance and supported jackknife conversions.
It does not control survey-bootstrap singleton handling. Bootstrap requires
at least two PSUs per stratum by default. `bootstrap_singleton = "certainty"`
holds a singleton PSU fixed and assigns zero first-stage variance; use it only
when that assumption is justified. WVW and KE survey jackknife require at
least two PSUs per stratum.

A failed survey replicate stops estimation by default. `relax = TRUE` allows
an explicitly exploratory adjustment that drops failed replicates; it does
not restore the full design. Inspect `diagnostics` and `raw$replication`.

The built-in design interface supports one final weight, optional strata,
and one PSU stage under an ultimate-cluster/with-replacement approximation.
It does not accept finite-population corrections, lower-stage identifiers,
first-stage PPS resampling, recalibration, or external replicate weights.

### Supplied ranks and reproducibility

For WVW and KE, `precalc_rank_var` supplies finite fractional ranks in `[0, 1]`.
KE also accepts `precalc_level_var`. Supplied positions take precedence over
`ses_var` and remain fixed for each sampled record in all replicates.
SES-derived positions are recomputed with the replicate weights.

An explicit `seed` makes stochastic calculations reproducible and restores
the caller's random-number state. With `seed = NULL`, stochastic calculations
advance the normal R random stream. Ordinary OBY, WVW, and KE jackknife is
deterministic; Fairlie jackknife still includes stochastic matching.

## Returned results

Every function returns an invisible named list with `summary_stats`,
`model_metrics`, `diagnostics`, and `raw`. Numeric results retain full
precision; display rounding does not alter the returned values.

| Function | Main result components |
|---|---|
| OBY | `results_overall`, `results_detailed_explained`, `results_detailed_unexplained`, `results_grouped_explained`, `results_grouped_unexplained`, `models` |
| Fairlie | `results_overall`, `results_detailed_explained`, `models_coefficients`; `results_detailed_unexplained` is `NULL` |
| WVW | `results_overall`, `results_detailed`, `results_grouped` |
| KE | `results_overall`, `results_detailed`, `results_anova` |

`level` sets the confidence level and must be one finite number strictly
between 0 and 1. Contribution percentages use the explained, unexplained,
or total component named in the column; they are not causal effect estimates.
OBY grouped tables are `NULL` unless `groupings` is supplied; WVW groups
contributions by predictor automatically when no custom grouping is supplied.

Choose the method for the target quantity, declare the groups or socioeconomic
ordering, select a supported variance estimator, and examine the sample,
model, matching, and replication diagnostics before interpreting results.
Report those choices together with the outcome bounds and any domain mask.

## Implementation credits

| Function | Method and implementation sources |
|---|---|
| `oby_decomp()` | Translation and adaptation to R of decomposition calculations in Ben Jann's [`oaxaca`](https://github.com/benjann/oaxaca) 4.1.1, based on Oaxaca, Blinder, Yun, and Jann. The R interface, diagnostics, and survey extensions are developed in `lunadecomp`. |
| `f_decomp()` | Fairlie's method, drawing on and adapting Ben Jann's [`fairlie`](https://github.com/benjann/fairlie) 1.0.7, including weighted matching and detailed variance calculations. Complex survey support and additional reference options are package extensions. |
| `wvw_decomp()` | Wagstaff, van Doorslaer, and Watanabe's method, guided by Chapter 13 and the `AHE13.do` example in O'Donnell et al.'s World Bank guide. Selected calculations were compared with `rineq` 0.3.0 and Stata `conindex` 1.7. |
| `ke_decomp()` | Direct implementation of Kessels and Erreygers (2019). The documented covariance estimators, survey support, and replication procedures are package extensions. |

Comparisons with other software establish agreement for the configurations
validated; they do not establish equivalence of every model or variance
estimator. Source revisions and the MIT notices for Ben Jann's code are
recorded in [COPYRIGHTS](inst/COPYRIGHTS), which is included in the installed
package.

Survey estimation uses `survey` (Lumley, 2004). Replicate-weight methods
follow Rao and Wu (1988) and Rao, Wu, and Yue (1992). The singleton-certainty
bootstrap uses `svrep` (Schneider, 2025) and the construction of Beaumont
and Émond (2022).

## Methodological references

- Oaxaca, R. (1973). Male-female wage differentials in urban labor markets. *International Economic Review*, 14(3), 693–709. [doi:10.2307/2525981](https://doi.org/10.2307/2525981).
- Blinder, A. S. (1973). Wage discrimination: Reduced form and structural estimates. *The Journal of Human Resources*, 8(4), 436–455. [doi:10.2307/144855](https://doi.org/10.2307/144855).
- Neumark, D. (1988). Employers' discriminatory behavior and the estimation of wage discrimination. *The Journal of Human Resources*, 23(3), 279–295. [doi:10.2307/145830](https://doi.org/10.2307/145830).
- Oaxaca, R. L., & Ransom, M. R. (1994). On discrimination and the decomposition of wage differentials. *Journal of Econometrics*, 61(1), 5–21. [doi:10.1016/0304-4076(94)90074-4](https://doi.org/10.1016/0304-4076(94)90074-4).
- Yun, M.-S. (2004). Decomposing differences in the first moment. *Economics Letters*, 82(2), 275–280. [doi:10.1016/j.econlet.2003.09.008](https://doi.org/10.1016/j.econlet.2003.09.008).
- Jann, B. (2008). The Blinder-Oaxaca decomposition for linear regression models. *The Stata Journal*, 8(4), 453–479. [doi:10.1177/1536867X0800800401](https://doi.org/10.1177/1536867X0800800401).
- Fairlie, R. W. (2005). An extension of the Blinder-Oaxaca decomposition technique to logit and probit models. *Journal of Economic and Social Measurement*, 30(4), 305–316. [doi:10.3233/JEM-2005-0259](https://doi.org/10.3233/JEM-2005-0259).
- Wagstaff, A., van Doorslaer, E., & Watanabe, N. (2003). On decomposing the causes of health sector inequalities with an application to malnutrition inequalities in Vietnam. *Journal of Econometrics*, 112(1), 207–223. [doi:10.1016/S0304-4076(02)00161-6](https://doi.org/10.1016/S0304-4076(02)00161-6).
- Wagstaff, A. (2005). The bounds of the concentration index when the variable of interest is binary, with an application to immunization inequality. *Health Economics*, 14(4), 429–432. [doi:10.1002/hec.953](https://doi.org/10.1002/hec.953).
- Erreygers, G. (2009). Correcting the concentration index. *Journal of Health Economics*, 28(2), 504–515. [doi:10.1016/j.jhealeco.2008.02.003](https://doi.org/10.1016/j.jhealeco.2008.02.003).
- Kjellsson, G., & Gerdtham, U.-G. (2013). On correcting the concentration index for binary variables. *Journal of Health Economics*, 32(4), 659–670. [doi:10.1016/j.jhealeco.2012.10.012](https://doi.org/10.1016/j.jhealeco.2012.10.012).
- Kessels, R., & Erreygers, G. (2019). A direct regression approach to decomposing socioeconomic inequality of health. *Health Economics*, 28(7), 884–905. [doi:10.1002/hec.3891](https://doi.org/10.1002/hec.3891).
- Yun, M.-S. (2005). A simple solution to the identification problem in detailed wage decompositions. *Economic Inquiry*, 43(4), 766–772. [doi:10.1093/ei/cbi053](https://doi.org/10.1093/ei/cbi053).
- Reimers, C. W. (1983). Labor market discrimination against Hispanic and Black men. *The Review of Economics and Statistics*, 65(4), 570–579. [doi:10.2307/1935925](https://doi.org/10.2307/1935925).
- Cotton, J. (1988). On the decomposition of wage differentials. *The Review of Economics and Statistics*, 70(2), 236–243. [doi:10.2307/1928307](https://doi.org/10.2307/1928307).
- O'Donnell, O., van Doorslaer, E., Wagstaff, A., & Lindelow, M. (2008). *Analyzing Health Equity Using Household Survey Data: A Guide to Techniques and Their Implementation*. World Bank, Chapter 13, 159–164. [doi:10.1596/978-0-8213-6933-3](https://doi.org/10.1596/978-0-8213-6933-3).
- O'Donnell, O., O'Neill, S., Van Ourti, T., & Walsh, B. (2016). conindex: Estimation of concentration indices. *The Stata Journal*, 16(1), 112–138. [Source](https://www.stata-journal.com/article.html?article=st0427).
- Lumley, T. (2004). Analysis of complex survey samples. *Journal of Statistical Software*, 9(1), 1–19. [doi:10.18637/jss.v009.i01](https://doi.org/10.18637/jss.v009.i01).
- Rao, J. N. K., & Wu, C. F. J. (1988). Resampling inference with complex survey data. *Journal of the American Statistical Association*, 83(401), 231–241. [doi:10.1080/01621459.1988.10478591](https://doi.org/10.1080/01621459.1988.10478591).
- Rao, J. N. K., Wu, C. F. J., & Yue, K. (1992). Some recent work on resampling methods for complex surveys. *Survey Methodology*, 18(2), 209–217. [Source](https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X199200214486).
- Beaumont, J.-F., & Émond, N. (2022). A bootstrap variance estimation method for multistage sampling and two-phase sampling when Poisson sampling is used at the second phase. *Stats*, 5(2), 339–357. [doi:10.3390/stats5020019](https://doi.org/10.3390/stats5020019).

## Software references

- Jann, B. (2023). *oaxaca: Stata module to compute the Blinder-Oaxaca decomposition*. Version 4.1.1 (24 April 2023). [Source](https://github.com/benjann/oaxaca).
- Jann, B. (2008). *fairlie: Stata module to compute nonlinear decomposition of binary outcome differentials*. Version 1.0.7 (16 June 2008); reference snapshot dated 27 February 2023. [Source](https://github.com/benjann/fairlie).
- Devleesschauwer, B., Willimès, S., Van Malderen, C., Konings, P., Speybroeck, N., & Meili, K. (2025). *rineq: Concentration Index and Decomposition for Health Inequalities*. R package version 0.3.0. [doi:10.32614/CRAN.package.rineq](https://doi.org/10.32614/CRAN.package.rineq).
- Schneider, B. (2025). *svrep: Tools for Creating, Updating, and Analyzing Survey Replicate Weights*. R package version 0.9.1. [doi:10.32614/CRAN.package.svrep](https://doi.org/10.32614/CRAN.package.svrep).


## Citation

If you use `lunadecomp` in academic work, technical reports, or policy analyses, please cite both:

1. The `lunadecomp` R package.
2. The methodological reference corresponding to the decomposition method used.

You can obtain the package citation in R with:

```r
citation("lunadecomp")
```

Recommended methodological citations depend on the function used:

| Function | Methodological framework | Recommended citation |
|---|---|---|
| `oby_decomp()` | Oaxaca-Blinder and Yun decomposition | Oaxaca (1973), Blinder (1973), Yun (2004, 2005), Jann (2008); Jann, `oaxaca` 4.1.1 |
| `f_decomp()` | Fairlie nonlinear decomposition | Fairlie (2005); Jann, `fairlie` 1.0.7 |
| `wvw_decomp()` | Wagstaff-type concentration index decomposition | Wagstaff, van Doorslaer, and Watanabe (2003); O'Donnell et al. (2008) |
| `ke_decomp()` | Kessels-Erreygers direct regression decomposition | Kessels and Erreygers (2019) |

The full methodological references are provided in the documentation of each
function. Cite Reimers (1983), Cotton (1988), or Neumark (1988) when using
the corresponding reference option, and the survey or resampling references
when those procedures are used.

## License

This project is licensed under the [GNU Affero General Public License v3.0
or later](LICENSE.md) (`AGPL-3.0-or-later`; `License: AGPL (>= 3)` in R).

Copyright and third-party notices are provided in [COPYRIGHTS](inst/COPYRIGHTS).
Adapted portions of Ben Jann's `oaxaca` and `fairlie` retain their MIT notices.

## Authors

- Adrián Vásquez-Mejía, MD, MSc — lead developer, principal author, and maintainer.
- Oscar J. Mujica, MD, MPH, PHE, FACE — conceptual contribution and scientific review.
- Antonio Sanhueza, MPH, MSc, PhD — conceptual contribution and scientific review.
