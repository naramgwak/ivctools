# ivctools

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/naramgwak/ivctools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/naramgwak/ivctools/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**ivctools** implements the **IV-Compass (IVC)** diagnostic for assessing the
exogeneity of a *single* instrument in pretest–posttest designs. Standard
overidentification tests (Sargan, Hansen J) require more than one instrument;
the IVC diagnostic instead exploits a baseline measure (the pretest) that is
already present in most educational and longitudinal data.

The idea: if the candidate instrument *Z* is valid, the conventional IV estimate
and a *compass-corrected* estimate built from the pretest should agree. Their
contrast, Δ = τ̂_IV − τ̂_comp, has expectation zero under a valid instrument; a
confidence interval for Δ that excludes zero is evidence against exogeneity.

## Installation

```r
# install.packages("remotes")
remotes::install_github("naramgwak/ivctools")
```

## Quick start

```r
library(ivctools)

# Simulated pretest–posttest data (kappa = 0.3 => mild exogeneity violation)
d <- ivc_simulate(n = 1000, gamma = 1, kappa = 0.3, seed = 1)

# The IVC diagnostic
fit <- compare_iv_cv(d, outcome = "Y", pretest = "P",
                     treat = "A", instrument = "Z", se = "mom")
fit
#> IV-Compass (IVC) exogeneity diagnostic
#>   tau_IV   = ...
#>   tau_comp = ...   (delta_hat = ...)
#>   Delta    = ...   (SE = ...)
#>   95% CI = [...]   p = ...
#>   Conclusion: interval excludes 0 -> evidence against exogeneity

plot(fit)
```

## Core functions

| Function | Purpose |
|---|---|
| `estimate_iv()`    | Conventional IV (Wald) estimate τ̂_IV |
| `estimate_cv()`    | Compass-corrected estimate τ̂_comp (and loading ratio δ̂) |
| `compare_iv_cv()`  | The IVC diagnostic: Δ with bootstrap or analytic MoM inference |
| `ivc_cli_test()`   | Conditional local independence (tetrad) test for candidate compass variables |
| `ivc_power()`      | Simulation-based power / Type-I-error analysis |
| `ivc_simulate()`   | Data-generating process for the single-factor model |
| `plot()` / `ivc_plot()` | Visualize the two estimates and the contrast |

Covariates and survey weights are supported in the estimators via
Frisch–Waugh–Lovell partialling (`covariates = `, `weights = `); with neither,
the estimators reduce exactly to the covariate-free forms.

For data with a natural grouping (e.g. students nested in schools), pass
`cluster = ` to `estimate_cv()` or `compare_iv_cv()` for cluster-robust
standard errors (analytic or bootstrap); point estimates are unaffected.
`estimate_cv()` also warns by default when the compass loading is too weak
for the correction (and its SE) to be trustworthy.

## Interpretation

The diagnostic converts a difficult credibility question into an equality
restriction that can be examined empirically. **Non-rejection is not proof of
validity.** The diagnostic is most informative when the instrument is relevant
and the pretest–posttest structure yields a stable compass correction; when
conditional relevance collapses, the contrast becomes noisy and intervals widen.

## Acknowledgment

Developed under the National Research Foundation of Korea, Humanities and Social
Sciences Academic Research Professor program (Type B, Growth Research):
"Developing an R-Based Analytical Toolkit for Detecting Instrumental Variable
Exogeneity Violations Using Pretest Scores."

## License

MIT © Naram Gwak
