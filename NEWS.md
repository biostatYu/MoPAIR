# MoPAIR 0.2.0

* Renamed the package and public API to MoPAIR: `fit_mopair()`,
  `fit_mopair_grid()`, `plot_mopair()`, and `mopair_example_data()`.
* Made the efficient logistic score test the default for binary outcome paths.
  `binomial_test = "wald"` retains the previous logistic Wald option.
* Added right-censored Cox outcome paths, configurable tie handling, and
  member-level hazard ratios with nominal 95% confidence limits.
* Preserved uncorrected DACT, equal-weight ACAT, and module-level BH selection.
* Aligned named fixed-module labels by mediator identifier and added explicit
  family-wide selection flags for multiple exposures.
* Recorded complete-case counts, member-specific sample sizes, unavailable fits,
  test type, coefficient scale, and score-based effect proxies.
* Added independent regression, Rao-score, and Cox agreement tests and examples
  for all three outcome families.

The new binary default changes outcome-path evidence relative to version 0.1.0.
Use the explicit Wald option when reproducing analyses that used that test.
