# ivctools 0.2.1

## Bug fixes and usability

* `compare_iv_cv(se = "bootstrap")` with `n_boot < 1000` silently returned
  `NA` standard errors because the minimum-valid-replicates floor was
  hard-coded at 1000; the floor now adapts to `n_boot`, and a genuine
  bootstrap failure produces an explicit warning instead of a silent `NA`.
* `ivc_cli_test()` now refuses a continuous treatment with a clear error
  (previously it produced one all-`NA` row per unique treatment value).
* `print.ivc_cli()` distinguishes three pair-level outcomes -- "violation
  detected", "inconclusive (some levels not testable)", and "no evidence of
  violation" -- instead of conflating rejection with untestability.

## Consistency

* `estimate_iv()` gains the same `cluster` argument as `estimate_cv()`.
* `compare_iv_cv()` gains the same `weak_loading_warn` check as
  `estimate_cv()`, so the main diagnostic entry point now warns when the
  compass loading is too weak for a stable correction. `ivc_power()`
  suppresses this per-replicate warning inside its simulation loop.
* `print.ivc()` reports when standard errors are cluster-robust and shows
  the number of clusters.

## Documentation

* Small-cluster caution (roughly G < 40; no finite-sample df correction)
  added to the `cluster` documentation of all three estimators.
* `ivc_cli_test()` documents the post-selection caveat of screening many
  candidate pairs, and notes that the tetrad SE assumes within-level
  independence (behavior under clustered data not yet established).

# ivctools 0.2.0

## New features

* `estimate_cv()` and `compare_iv_cv()` gain a `cluster` argument. With
  `se = TRUE` / `se = "mom"` the sandwich "meat" is built from cluster-summed
  moment contributions (cluster-robust SE); with `se = "bootstrap"` whole
  clusters are resampled. Implements the multilevel extension of
  Lee & Kim (2026), *JEEV*, 39(2), 267-291. Point estimates are unchanged.
* New `ivc_cli_test()`: conditional local independence (tetrad) test for
  candidate compass variables, with delta-method SEs and Bonferroni
  adjustment across treatment levels; screens all pairs when more than two
  candidates are supplied. Implements Lee & Kim (2026), *JEEV*, 39(1), 55-75.
* `estimate_cv()` gains `weak_loading_warn`: warns when the compass loading
  denominator cov(P, Z | A) is statistically indistinguishable from zero,
  in which case the ratio `delta_hat` and its SE are unreliable.

## Bug fixes

* `.ivc_mom_se()`: the moment functions for `tau_comp` and `tau_IV` omitted
  the intercept (they were the moments of no-intercept estimators). The error
  was invisible when the treatment and instrument had mean zero, but it
  distorted standard errors for binary treatments. Moments are now weighted
  and centered; the IV SE now matches the textbook HC0 formula exactly.

# ivctools 0.1.0

* Initial release.
* Core estimators: `estimate_iv()`, `estimate_cv()`, `compare_iv_cv()`.
* Simulation and power analysis: `ivc_simulate()`, `ivc_power()`.
* Visualization: `plot()` method / `ivc_plot()`.
