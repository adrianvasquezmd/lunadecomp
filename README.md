# lunadecomp

<!-- badges: start -->
[![R-CMD-check](https://github.com/adrianvasquezmd/lunadecomp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adrianvasquezmd/lunadecomp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`lunadecomp` is an R package for applied decomposition analysis in health inequality research. It provides a unified, reproducible, and survey-aware implementation of regression-based decomposition methods used to quantify and explain socioeconomic differences in health outcomes.

The package is designed for epidemiology, public health, health economics, and health inequality monitoring workflows where analysts need to decompose either:

1. **Two-group outcome gaps**, such as differences in health outcomes between population groups; or
2. **Socioeconomic health inequality indices**, such as concentration-type indices and corrected bounded inequality measures.

`lunadecomp` emphasizes transparent statistical workflows, interpretable outputs, support for complex survey designs, and CRAN-compatible documentation, examples, and tests.

## Implemented methods

| Function | Methodological framework | Main estimand | Typical use case |
|---|---|---|---|
| `oby_decomp()` | Oaxaca-Blinder and Yun decomposition | Two-group outcome gap | Decompose mean or probability differences between two groups into explained and unexplained components. |
| `f_decomp()` | Fairlie nonlinear decomposition | Binary outcome gap | Decompose binary outcome differences using logit or probit counterfactual predictions. |
| `wvw_decomp()` | Wagstaff, van Doorslaer, and Watanabe decomposition | Socioeconomic concentration-type index | Decompose health inequality into determinant-specific contributions based on marginal effects, elasticities, and concentration indices. |
| `ke_decomp()` | Kessels-Erreygers direct regression approach | Rank- or level-dependent bivariate inequality index | Estimate direct-regression marginal effects and variable importance for composite socioeconomic-health inequality targets. |

## Installation

You can install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("adrianvasquezmd/lunadecomp")
```

Then load the package:

```r
library(lunadecomp)
```

## Example dataset

The package includes a simulated dataset for examples, tests, and vignettes:

```r
data(lunadecomp_example)
str(lunadecomp_example)
```

The dataset contains individual-level observations with:

- a two-group comparison variable, `group`;
- a socioeconomic variable, `ses`;
- survey design variables, `weight`, `strata`, and `psu`;
- common covariates such as age, rural residence, insurance, sanitation, and education;
- continuous, binary, and bounded health-related outcomes.

This dataset is simulated and should be used only for demonstration, unit tests, and teaching examples.

## Quick start

### Oaxaca-Blinder / Yun decomposition

Use `oby_decomp()` to decompose a two-group outcome gap into explained and unexplained components.

```r
fit_ob <- oby_decomp(
  data = lunadecomp_example,
  dep_var = "y_continuous",
  group_var = "group",
  indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
  model_type = "ols",
  ref_method = "pooled",
  vce_method = "linearized",
  quiet = TRUE
)

fit_ob$results_overall
fit_ob$results_detailed_explained
fit_ob$results_detailed_unexplained
```

For binary outcomes, the function also supports:

```r
model_type = "lpm"    # Linear probability model
model_type = "logit"  # Yun-type nonlinear decomposition
model_type = "probit" # Yun-type nonlinear decomposition
```

Declare the comparison order explicitly with the favored level first and the
disadvantaged level second:

```r
group_levels = c("highest_quintile", "lowest_quintile")
```

`group_var` may contain additional levels. `oby_decomp()` retains the two
selected values, filters all other levels, and reports
`Difference = E(Y | favored) - E(Y | disadvantaged)`. Reversing the vector
reverses the comparison. The older `favored_group` argument remains available
when the data already contain exactly two groups.

The detailed output includes term-level estimates, uncertainty intervals, and percentage contributions to the explained, unexplained, and total gaps.

For categorical predictors, detailed results can be normalized so they do not
depend on the omitted category. Normalize selected factors by name:

```r
normalize = c("insurance", "education")
```

or normalize every factor included in `indep_vars`:

```r
normalize = "all"
```

Normalization uses equal-category deviation contrasts, retains one detailed
row for every level (including the originally omitted level), and propagates
the full covariance matrix. `normalize = NULL`, the default, retains ordinary
treatment coding.

The same argument automatically normalizes factor-by-continuous interactions:

```r
fit_interaction <- oby_decomp(
  data = analysis_data,
  dep_var = "outcome",
  group_var = "group",
  indep_vars = c("education * experience", "age"),
  normalize = "education"
)
```

Here the main effects of `education` are centered around the intercept and the
`education:experience` coefficients are centered independently around the
coefficient of `experience`. The specification must be hierarchical and use
one factor with one numeric continuous variable per normalized interaction.

### Fairlie nonlinear decomposition

Use `f_decomp()` for binary outcome gaps estimated using logit or probit models.

```r
fit_fairlie <- f_decomp(
  data = lunadecomp_example,
  dep_var = "y_binary",
  group_var = "group",
  indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
  model_type = "logit",
  ref_method = "pooled",
  vce_method = "linearized",
  reps = 25,
  quiet = TRUE
)

fit_fairlie$results_summary
fit_fairlie$results_detailed
```

The detailed Fairlie table reports both:

- `% Contribution Explained`: contribution divided by the total explained component;
- `% Contribution Total`: contribution divided by the observed total gap.

For a survey bootstrap in `oby_decomp()`, `f_decomp()`, `wvw_decomp()`, or
`ke_decomp()`, `lonely_psu` does not control Rao--Wu replicate construction.
If a stratum contains one PSU and that PSU can defensibly be treated as
selected with certainty, use `bootstrap_singleton = "certainty"`. The
singleton is then held fixed with replicate factor 1 and contributes zero
first-stage variance. The default, `bootstrap_singleton = "fail"`, requires at
least two PSUs per stratum and stops instead of making that assumption
silently.

### Wagstaff-type decomposition

Use `wvw_decomp()` to decompose a socioeconomic concentration-type index into determinant-specific contributions.

```r
fit_wvw <- wvw_decomp(
  data = lunadecomp_example,
  dep_var = "health_score",
  indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
  ses_var = "ses",
  model_type = "ols",
  correction = "standard",
  vce_method = "linearized",
  quiet = TRUE
)

fit_wvw$results_overall
fit_wvw$results_detailed
fit_wvw$results_grouped
```

The detailed table reports:

- `Marginal Effects`;
- `Elasticity`;
- `Concentration Index`;
- `Contribution`;
- `% Contribution Explained`;
- `% Contribution Total`.

Available correction methods include:

```r
correction = "standard"     # Relative concentration index
correction = "generalized"  # Generalized or absolute concentration index
correction = "erreygers"    # Erreygers correction for bounded outcomes
correction = "wagstaff"     # Wagstaff normalization for bounded outcomes
```

For bounded corrections, provide theoretical outcome bounds whenever possible:

```r
correction = "erreygers"
dep_min = 0
dep_max = 1
```

`wvw_decomp()` uses one VCE interface:

```r
vce_method = "linearized"
vce_method = "jackknife"
vce_method = "bootstrap"
```

Without survey design variables these select the observation-level
Huber--White/Taylor influence, delete-one jackknife, and ordinary multinomial
bootstrap. With strata or PSU, they select survey Taylor, JK1/JKn, and the
Rao--Wu rescaled bootstrap. JKn uses stratum-specific replicate centers and
survey inference uses `number of PSUs - number of strata` degrees of freedom.

The built-in WVW survey contract represents one final weight, optional strata,
and one PSU stage under an ultimate-cluster/with-replacement approximation.
FPC, lower stages, first-stage PPS resampling, recalibration, and external
replicate weights require a future external survey-design interface and are
not approximated silently.

### Kessels-Erreygers direct regression decomposition

Use `ke_decomp()` for the direct regression approach to rank-dependent or level-dependent socioeconomic health inequality indices.

```r
fit_ke <- ke_decomp(
  data = lunadecomp_example,
  dep_var = "health_score",
  indep_vars = c("age", "rural", "insurance", "sanitation", "education"),
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

Unlike the Wagstaff decomposition, the Kessels-Erreygers direct regression approach should not be interpreted as a percentage-contribution decomposition of the observed index. The reported coefficients are marginal effects on the constructed individual-level composite target. Variable importance is summarized using F tests and logworth values.

Factors use their first R factor level as the reference. Formula terms and
interactions can be supplied directly:

```r
fit_ke_interaction <- ke_decomp(
  data = survey_data,
  dep_var = "health",
  indep_vars = c("age", "education", "age:education"),
  ses_var = "wealth_index",
  index_type = "rank",
  correction = "generalized",
  vce_method = "linearized",
  quiet = TRUE
)
```

For ordinary independent observations, `vce_method = "linearized"` uses the
HC1 Huber–White covariance equivalent to Stata `regress, vce(robust)`.
This coefficient inference is conditional on the constructed KE target and is
an explicit package extension: Kessels and Erreygers (2019) do not specify
their covariance estimator. `jackknife` and `bootstrap` reconstruct ranks or
relative levels, means, correction scale, target, and regression in every
replicate, but their covariance also applies only to the regression
coefficients.

For ordinary resampling, bootstrap draws analytic rows with replacement through
`boot::boot()` and uses the sample covariance of the valid coefficient
replicates. Jackknife uses the deterministic delete-one samples from
`resample::jackknife()` and the non-MSE factor `(N - 1) / N`. Both are centered
at the replicate mean, matching Stata's ordinary prefix conventions. If a
replicate is undefined, the complete coefficient vector is omitted and the
scale is recomputed from the valid count. An explicit `seed` makes bootstrap
reproducible without changing the caller's RNG state; `seed = NULL` advances
the normal R random stream.

When `precalc_rank_var` or `precalc_level_var` is supplied, that position is
treated as external record-level information and remains fixed inside
jackknife and bootstrap replicates. The sampled records retain their supplied
position while health means, correction scales, targets, and regressions are
recomputed. When position is instead constructed from `ses_var`, it is
recomputed from the replicate weights. A supplied precalculated position takes
precedence over `ses_var`.

Required numeric inputs must be finite. Sampling weights must be strictly
positive, and numeric weights stored as text are parsed by value rather than
as factor codes. Missing required values are removed jointly and their source
rows are available in `raw$sample`; invalid infinite values are rejected.

In survey mode, `linearized` uses the Taylor design covariance of the
direct-regression coefficients conditional on the constructed KE target.
Supplying a weight, stratum, or PSU variable automatically activates survey
mode. Designs without an explicit weight use unit weights, and designs without
an explicit PSU treat each analytic row as a PSU. Coefficient t tests,
confidence intervals, and F tests use design degrees of freedom
`nested PSUs - strata`, matching Stata `svy`. This covariance accounts for the
declared design but does not propagate the sampling variability of the
empirical ranks, relative levels, means, or correction scale.

In survey replicate mode, jackknife uses JK1 without strata and JKn with
strata. JKn centers each block on the mean of the replicates that alter the
same stratum. Bootstrap uses Rao--Wu rescaled PSU sampling
(`survey::as.svrepdesign(type = "subbootstrap")`) and the Stata `svy
bootstrap` non-MSE scale `1 / B`; the generator's native `1 / (B - 1)` scale
is retained in `raw` for auditability. Both methods reconstruct ranks or
levels, means, correction scale, target, and coefficients in every replicate
and use `nested PSUs - strata` degrees of freedom. At least two PSUs are
required per stratum, and a failed survey replicate stops estimation unless
`relax = TRUE` is explicitly requested for an exploratory approximation.

The article reports the aggregate index as a point estimate without a standard
error, test, or confidence interval. `ke_decomp()` follows that scope under
every VCE method: the index `Estimate` is reported and its harmonized
inferential fields remain `NA`. The package does not calculate a replicate or
survey variance for the aggregate index.

For `correction = "erreygers"` or `"wagstaff"`, `dep_min` and `dep_max`
are required theoretical outcome bounds; observed sample extrema are not used
as substitutes. Wagstaff normalization is available only with
`index_type = "rank"`. Although corrected aggregate indices retain their
affine-invariance or mirror properties, the KE direct-regression coefficients
need not remain invariant after translating or complementing the outcome.

The completed ten-phase KE validation and final capability matrix distinguish
article-defined calculations, Stata-composed comparisons, index
triangulation, algebraic checks, simulation diagnostics, and unsupported
survey features. PPS, FPC, lower stages, calibration, external replicate
weights, and aggregate-index inference remain explicitly outside the current
contract.

## Conceptual distinction between methods

### Two-group gap decomposition

`oby_decomp()` and `f_decomp()` are used when the target is a difference between two groups:

```text
Observed gap = outcome in group 0 - outcome in group 1
```

These methods partition the observed gap into explained and unexplained components. The explained component reflects differences in observed covariates. The unexplained component reflects coefficient differences and residual structure. It should not be interpreted mechanically as discrimination, inequity, or causal effect without additional substantive assumptions.

### Socioeconomic inequality index decomposition

`wvw_decomp()` and `ke_decomp()` are used when the target is a socioeconomic inequality index rather than a two-group gap.

`wvw_decomp()` follows the conventional concentration-index decomposition logic, where determinant-specific contributions are based on marginal effects, covariate means, and covariate concentration indices.

`ke_decomp()` follows the direct regression logic of Kessels and Erreygers, where the inequality index is reformulated as the mean of an individual-level composite variable, and that composite variable is regressed directly on explanatory variables.

## Survey design support

The decomposition functions expose optional survey design variables where
supported by their documented method:

```r
weight_var = "weight"
strata_var = "strata"
psu_var = "psu"
use_svy = TRUE
```

For `wvw_decomp()`, supplying any of `weight_var`, `strata_var`, or `psu_var`
activates survey-aware behavior automatically. Reused PSU labels are nested
within strata.

For `oby_decomp()`, `wvw_decomp()`, and `ke_decomp()`, the public
variance-estimator names are:

```r
vce_method = "linearized"        # Huber-White/Taylor as documented
vce_method = "jackknife"
vce_method = "bootstrap"
```

With a declared survey design, their bootstrap uses Rao--Wu rescaled
PSU-within-stratum replicates. Each stratum must contain at least two PSUs.
Outside survey mode, `oby_decomp()` generates deterministic delete-one
jackknife samples with `resample::jackknife()` and retains the Stata-compatible
JK1 scale and inference.

Other decomposition functions retain their documented method-specific VCE
interfaces and may reject combinations for which a complete design-based VCE
has not been established. Always report the exact design and VCE rather than
treating the three labels as interchangeable across functions.

## Output structure

Most functions return an invisible named list. Common components include:

| Component | Description |
|---|---|
| `summary_stats` | Sample size, population size, degrees of freedom, and design information. |
| `model_metrics` | Model fit statistics when applicable. |
| `diagnostics` | Validation and numerical diagnostics. |
| `results_overall` or `results_summary` | Main decomposition summary. |
| `results_detailed` | Term-level decomposition or marginal-effect table. |
| `results_grouped` | Grouped or domain-level results when applicable. |
| `models` or `models_coefficients` | Fitted model objects or coefficient tables when applicable. |

## Recommended workflow

A typical analysis workflow is:

1. Clean and label the input dataset.
2. Select and order the comparison with
   `group_levels = c(favored, disadvantaged)`; binary outcomes still require
   0/1 coding.
3. Encode categorical predictors as factors.
4. Choose the decomposition method based on the estimand.
5. Start with the function's linearized VCE for exploratory runs.
6. Use `jackknife` or `bootstrap` for final inference when feasible.
7. Review diagnostics before interpreting results.
8. Report the method, reference structure, variance estimator, survey design, and any bounded-index correction.

## Methodological references

- Oaxaca, R. (1973). Male-female wage differentials in urban labor markets. *International Economic Review*, 14(3), 693. doi:10.2307/2525981.
- Blinder, A. S. (1973). Wage discrimination: Reduced form and structural estimates. *The Journal of Human Resources*, 8(4), 436. doi:10.2307/144855.
- Neumark, D. (1988). Employers' discriminatory behavior and the estimation of wage discrimination. *The Journal of Human Resources*, 23(3), 279. doi:10.2307/145830.
- Oaxaca, R. L., & Ransom, M. R. (1994). On discrimination and the decomposition of wage differentials. *Journal of Econometrics*, 61(1), 5-21. doi:10.1016/0304-4076(94)90074-4.
- Yun, M.-S. (2004). Decomposing differences in the first moment. *Economics Letters*, 82(2), 275-280. doi:10.1016/j.econlet.2003.09.008.
- Jann, B. (2008). The Blinder-Oaxaca decomposition for linear regression models. *The Stata Journal*, 8(4), 453-479. doi:10.1177/1536867X0800800401.
- Fairlie, R. W. (2005). An extension of the Blinder-Oaxaca decomposition technique to logit and probit models. *Journal of Economic and Social Measurement*, 30(4), 305-316. doi:10.3233/JEM-2005-0259.
- Wagstaff, A., van Doorslaer, E., & Watanabe, N. (2003). On decomposing the causes of health sector inequalities with an application to malnutrition inequalities in Vietnam. *Journal of Econometrics*, 112(1), 207-223. doi:10.1016/S0304-4076(02)00161-6.
- Wagstaff, A. (2005). The bounds of the concentration index when the variable of interest is binary, with an application to immunization inequality. *Health Economics*, 14(4), 429-432. doi:10.1002/hec.953.
- Erreygers, G. (2009). Correcting the concentration index. *Journal of Health Economics*, 28(2), 504-515. doi:10.1016/j.jhealeco.2008.02.003.
- Kjellsson, G., & Gerdtham, U.-G. (2013). On correcting the concentration index for binary variables. *Journal of Health Economics*, 32(4), 659-670. doi:10.1016/j.jhealeco.2012.10.012.
- Kessels, R., & Erreygers, G. (2019). A direct regression approach to decomposing socioeconomic inequality of health. *Health Economics*, 28(7), 884-905. doi:10.1002/hec.3891.

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
| `oby_decomp()` | Oaxaca-Blinder and Yun decomposition | Oaxaca (1973), Blinder (1973), Yun (2004, 2005), Jann (2008) |
| `f_decomp()` | Fairlie nonlinear decomposition | Fairlie (2005) |
| `wvw_decomp()` | Wagstaff-type concentration index decomposition | Wagstaff, van Doorslaer, and Watanabe (2003) |
| `ke_decomp()` | Kessels-Erreygers direct regression decomposition | Kessels and Erreygers (2019) |

The full methodological references are provided in the documentation of each function.

## Development status

`lunadecomp` is under active development. The current version is intended for methodological development, reproducible research workflows, and applied health inequality analyses. Interfaces and output names may evolve before a stable release.

The package currently passes:

```text
R CMD check: 0 errors, 0 warnings, 0 notes
Unit tests: 1,228 passing expectations in 163 test blocks
WVW validation: 16,074 accepted cross-reference checks
```

## License

This project is licensed under the GNU Affero General Public License v3.0 or later.

## Authors

- Adrián Vásquez-Mejía, MD, MSc
- Oscar J. Mujica, MD, MPH, PHE, FACE
- Antonio Sanhueza, MPH, MSc, PhD
