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
