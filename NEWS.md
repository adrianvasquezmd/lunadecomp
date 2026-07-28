# lunadecomp 0.1.0

* Initial release.
* Added `oby_decomp()` for Oaxaca-Blinder and Yun two-group outcome gap decomposition (OLS, LPM, logit, and probit models).
* Added `f_decomp()` for Fairlie nonlinear decomposition of binary outcome gaps.
* Added `wvw_decomp()` for Wagstaff, van Doorslaer, and Watanabe concentration-index decomposition, with standard, generalized, Erreygers, and Wagstaff corrections.
* `wvw_decomp()` now returns unrounded public numeric tables and a versioned
  `raw` audit block containing analytic-sample rows, ranks, normalized design
  weights, model matrices and covariance, determinant components, estimates,
  complete supported covariance matrices, influence functions, and replication
  internals. Design boundaries and conditional benchmarks are labeled
  explicitly.
* Validated the unweighted, no-tie OLS point decomposition in `wvw_decomp()`
  against Stata `regress`, composed Stata covariance primitives, `conindex`
  1.7, and `rineq` 0.3.0. Coefficients, determinant indices, elasticities,
  contributions, total, explained component, residual, algebraic closure, row
  order, determinant/outcome scale, and affine SES-rank invariance are covered
  by permanent tests.
* Validated positive unequal-weight OLS point decompositions against Stata
  `pweight`, `aweight`, `fweight`, `conindex` 1.7, explicit population
  covariances, and direct `rineq` primitives. The production convention is
  probability-normalized and invariant to global weight scaling; Stata's
  historical `fweight` covariance and the distinct finite-sample correction
  in weighted `rineq::contribution.lm()` are retained only as documented
  comparison conventions.
* Completed the eleven-phase `wvw_decomp()` validation program. The frozen
  evidence contains 16,074 accepted cross-reference checks spanning the
  historical World Bank example, Stata 19.5, `conindex`, `rineq`, numerical
  influence derivatives, shared replicate weights, and algebraic identities.
* Completed joint linearized inference for WVW contributions, total,
  residual, and nonlinear approximation error. The public covariance propagates
  fitted coefficients, average marginal effects, determinant concentration
  covariances, outcome mean, correction scale, and empirical fractional ranks.
* Aligned WVW ordinary replication with Stata's non-MSE conventions:
  observation delete-one jackknife with `(N-1)/N`, multinomial bootstrap with
  `1/(B-1)`, replicate-mean centering, and valid-replicate handling.
* Aligned WVW survey inference under the documented ultimate-cluster contract:
  Taylor covariance, JK1/JKn with stratum-specific JKn centers, Rao--Wu
  rescaled bootstrap with final `1/B`, nested PSU labels, `PSU - strata`
  degrees of freedom, and Stata-equivalent lonely-PSU rules.
* Survey-design extensions not represented by the WVW interface—FPC,
  lower-stage identifiers, first-stage PPS resampling, recalibration, and
  external replicate weights—are now documented as unsupported rather than
  approximated silently.
* Added a final capability matrix, phase-level evidence audit, and frozen
  artifact hashes. Local bibliographic material, validation artifacts, and
  reference source trees are excluded from the distributable R package while
  remaining available in the development workspace.
* Added `ke_decomp()` for the Kessels-Erreygers direct regression decomposition of socioeconomic health inequality indices.
* Validated `ke_decomp()` ordinary bootstrap and delete-one jackknife against
  independent reconstruction and Stata 19.5 for rank- and level-dependent
  targets under generalized, standard, Erreygers, and rank-only Wagstaff
  corrections. Ordinary replication uses `boot::boot()` and
  `resample::jackknife()`, replicate-mean centering, Stata non-MSE scaling, and
  Stata-compatible omission and rescaling of failed coefficient replicates.
  The aggregate index remains point-only under both methods.
* Validated `ke_decomp()` weight-only and weight--stratum--PSU Taylor
  coefficient covariance against an independent PSU-score linearization and
  Stata 19.5 `svy: regress`. Any supplied design variable now activates survey
  mode, unit weights are used when no weight is declared, and survey
  coefficient inference uses design degrees of freedom (`PSUs - strata`)
  instead of `svyglm` model-residual degrees. The linearized covariance remains
  explicitly conditional on the constructed KE target.
* Validated `ke_decomp()` survey jackknife and bootstrap by reconstructing all
  1,368 coefficient replicates independently and in Stata 19.5 with identical
  replicate weights. Survey jackknife now selects JK1 without strata and JKn
  with strata; JKn uses stratum-specific non-MSE centers. Survey bootstrap now
  uses Rao--Wu `subbootstrap` weights and Stata's final `1/B` scale while
  retaining the native generator scale for auditability. Both methods use
  `PSUs - strata` degrees of freedom, reject singleton strata, and stop rather
  than silently discard failed survey replicates.
* Completed the ten-phase `ke_decomp()` validation program. Final edge-case
  checks now enforce scalar API arguments, finite numeric inputs, strictly
  positive sampling weights, explicit degenerate-factor errors, and stable
  constant-target diagnostics. Numeric weights stored as text are parsed by
  value instead of factor code.
* Defined precalculated KE ranks and relative levels as external positions that
  remain fixed inside ordinary and survey jackknife/bootstrap replicates.
  Eight replicate scenarios agree with independent reconstruction, while
  SES-derived positions continue to be recomputed from every replicate's
  weights. The final local capability matrix labels article-defined,
  Stata-composed, triangulated, algebraic, simulation-only, and unsupported
  features.
* Added the simulated `lunadecomp_example` dataset for examples, tests, and vignettes.
* Support for complex survey designs (weights, strata, PSU) and linearized, jackknife, and bootstrap variance estimation across all core functions.
* Aligned `oby_decomp()` linearized point estimates and robust delta-method
  covariance with Ben Jann's `oaxaca, vce(robust)` for the six shared two-fold
  reference structures.
* `oby_decomp()` now preserves unrounded numeric output and exposes raw
  estimates, covariance matrices, coefficients, covariate means, and analytic
  sample indices.
* The `oby_decomp()` VCE interface now consists of `linearized`, `jackknife`,
  and `bootstrap`. Linearized estimation uses the Huber--White sandwich outside
  survey mode and Taylor design linearization when survey mode is active.
* Documented and validated binary and multi-level categorical predictors,
  including grouped contributions and covariance-aware grouped standard errors.
* Added Yun's equal-category normalization for additive factor predictors in
  `oby_decomp()`. Users may select factors individually with
  `normalize = c("factor1", "factor2")` or normalize every modeled factor with
  `normalize = "all"`. Coefficients, covariate means, detailed contributions,
  grouped results, and their full covariance matrices are transformed
  consistently, including within jackknife and bootstrap replicates.
* Extended the same `normalize` interface to factor-by-continuous interactions.
  `indep_vars` now accepts R formula terms such as
  `"education * experience"`; selecting `normalize = "education"`
  automatically normalizes both its main effect and its interaction block,
  using the continuous coefficient as the interaction anchor.
* Validated LPM, Logit, and Probit two-fold decompositions against Ben Jann's
  `oaxaca` 4.1.1 across all six shared reference structures, including
  normalized factors and factor-by-continuous interactions.
* Aligned nonlinear linearized VCEs with Stata's maximum-likelihood sandwich
  conventions: tight binary-model convergence, the `N/(N-1)` finite-sample
  correction for group Logit/Probit equations, and the observed-information
  Hessian for Probit.
* Added `group_levels = c(favored, disadvantaged)` to `oby_decomp()`. It
  explicitly selects and orders two comparison levels, filters other levels of
  a multilevel group variable, records the full selection mapping, and makes
  the orientation of the reported gap unambiguous. The earlier
  `favored_group` interface remains supported for two-level data.
* Aligned `oby_decomp()` linearized inference with Stata: normal-asymptotic
  tests outside survey mode and design-based t tests with `PSU - strata`
  degrees of freedom in survey mode. Two-sided p-values now use a numerically
  stable direct upper-tail calculation, preserving probabilities far below
  machine epsilon instead of rounding them to zero.
* Aligned jackknife and bootstrap replication with Stata conventions:
  delete-one/JKn centering, fixed full-sample Cotton weights, ordinary
  observation bootstrap, survey PSU-within-stratum bootstrap weights,
  `bsn(1)` scale `1/B`, method-appropriate degrees of freedom, and unrounded
  replicate diagnostics in `raw$replication`.
* Began direct validation of `f_decomp()` against Ben Jann's `fairlie` 1.0.7.
  The deterministic unweighted core now matches Logit and Probit for
  `group0`/`reference(0)` and `group1`/`reference(1)`, including coefficients,
  individual predictions, total explained effects, and detailed
  contributions.
* `f_decomp()` now returns unrounded public tables and a structured `raw`
  component containing estimates, available covariance matrices, model
  coefficients, predictions, analytic-sample rows, decomposition blocks,
  matching diagnostics, and estimation settings.
* Detailed Fairlie contributions are no longer rescaled to the full-sample
  explained total. The matched-sample sum and its difference from the
  full-sample total are retained for auditability.
* `f_decomp()` now estimates only the models required by the selected
  reference, preventing separation or convergence failures in unused group
  models.
* Added `group_levels = c(favored, disadvantaged)` to `f_decomp()`, including
  explicit selection of two levels from a multilevel group variable, internal
  0/1 recoding, filtering metadata, and an auditable group mapping.
* Added pooled counterfactual anchors to `f_decomp()`:
  `pooled_anchor = "favored"`, `"disadvantaged"`, or `"centered"`. Logit and
  Probit results for all three anchors and both group orientations were
  validated against `fairlie, pooled(varlist)`. Reimers and Cotton remain
  available as package extensions without a native `fairlie` comparator.
* Aligned unequal-group matching in `f_decomp()` with Ben Jann's RNG logic:
  the complete smaller group is no longer needlessly permuted, while a ranked
  simple random sample without replacement is drawn from the larger group.
  Matching metadata now records the method and sampled group.
* Added `raw$matching$detailed_monte_carlo_se`, which measures finite-
  replication algorithmic error separately from statistical sampling
  uncertainty. Unequal-group Logit and Probit matching was validated against
  shared-draw harnesses and exact combinatorial expectations.
* Validated `f_decomp()` fixed and randomized block order against `fairlie`
  using binary indicators and a three-level factor. Automatic R factor
  expansion is exactly equivalent to manually generated and jointly grouped
  dummy variables.
* Hardened `groupings`: unmatched specifications, overlapping blocks, unnamed
  entries, and splitting the columns of a multi-level factor now stop with
  explicit errors. This prevents duplicated swaps and invalid categorical
  counterfactual profiles.
* Validated `f_decomp()` probability weights against `fairlie [pweight=...]`
  for Logit and Probit, including weighted means and models, pooled and
  group-specific references, globally rescaled weights, group-imbalanced
  weights, and zero-weight rows.
* Weighted Fairlie matching now has an explicit audited contract: both ranked
  groups are sampled independently with replacement, proportional to their
  weights, using `floor((N0 + N1) / 2)` matches per replication. Weight
  metadata and the separate detailed Monte Carlo error are returned under
  `raw$weighting` and `raw$matching`.
* `f_decomp()` now rejects negative, non-finite, and nonnumeric weights with
  clear errors, excludes and records zero-weight observations, and verifies
  that both selected groups retain positive-weight observations. A global
  positive rescaling of weights leaves all point estimates invariant.
* Aligned `f_decomp()` linearized detailed inference with
  `fairlie, vce(robust)`. Outside a complex survey, `linearized` now always
  uses the Huber--White maximum-likelihood sandwich, Stata's `N/(N-1)`
  correction for each Logit/Probit equation, the observed-information Probit
  Hessian, and normal-asymptotic inference.
* Fairlie block gradients and their per-matching delta variances are now
  reproduced directly. The detailed covariance remains diagonal, matching
  Ben Jann's documented convention, and is explicitly marked incomplete to
  prevent unsupported joint interpretation.
* Completed the ten-phase `f_decomp()` validation program. Ordinary
  jackknife/bootstrap and survey JKn/Rao--Wu replication now return a complete
  covariance for the replicated estimator, while incomplete complex-survey
  Taylor linearization is rejected explicitly.
* Validated the Reimers and Cotton `f_decomp()` extensions through exact
  coefficient, covariance, orientation, and equal-weight identities. Cotton's
  full-sample weighted group share is now fixed across sampling replicates, so
  linearized and replicate inference target the same conditional reference.
* Clarified the interaction contract for Fairlie decomposition: formula
  interactions are rejected, while precomputed interactions are valid only
  when the complete interacting system is exchanged as one joint block.
* Replaced the misleading “propensity overlap” wording with a descriptive
  predicted-outcome-probability overlap diagnostic. Its ranges and weighted
  and unweighted off-overlap proportions are returned without trimming any
  observations.
