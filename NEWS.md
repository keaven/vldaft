# vldaft 0.2.0

## Features

* Added generic Poisson-mixture cure-model wrappers for the `gamlss2`
  family layer: `Cure_Weibull()`, `Cure_Logistic()`, `Cure_Normal()`,
  `Cure_Cauchy()`, `Cure_Gamma()`, and `Cure_Exponential()`.
* Added a cure-model vignette (`vignettes/cure-mixture.Rmd`) describing
  the likelihood, first/second partials, simulation setup, and the new
  family constructors.
* Added a motivational vignette (`vignettes/weibull-motivation.Rmd`)
  with interactive `gglite` visualizations showing how location, scale,
  the VLDAFT coupling parameter, and cure rate shape Weibull survival
  and hazard functions.
* Added a companion slide deck (`slides/vldaft-weibull-demo.qmd`) for
  presenting the Weibull AFT and cure-model methodology to a statistical
  audience.

## Internal changes

* Refactored `R/family_aft.R` so the baseline AFT density/survival/score
  helpers can be reused by both the plain AFT families and the new
  cure-model wrappers.
* Added self-contained simulation helpers and test coverage for the
  cure-family likelihood, scores, left-truncation handling, and
  intercept-only parameter recovery.

# vldaft 0.1.0

First public release.

## Features

* `vldaft()` fits the Anderson (1991) varying-location-and-dispersion
  accelerated-failure-time (VLDAFT) model with Weibull, logistic, normal,
  Cauchy, or gamma error distributions, right-censored and/or left-truncated
  data, and an optional polynomial coupling of `log(sigma)` to the fitted
  location `mu`.
* Dual computational back-ends selectable via `backend = c("c", "rust")`;
  the Rust back-end uses `extendr-api` and keeps `f64` precision end-to-end.
* S3 methods: `print`, `summary`, `coef`, `vcov`, `logLik` (with optional
  Jacobian correction), and `predict` (linear predictor, log-sigma, sigma,
  predicted median, predicted survival).
* gamlss2 family wrappers: `AFT_Weibull()`, `AFT_Logistic()`, `AFT_Normal()`,
  `AFT_Cauchy()`, `AFT_Gamma()`.

## Internal changes compared with the pre-release prototype

* The C data path now uses double precision throughout (the prototype stored
  the design matrix in `float`). This removes the ~1e-7 gap between the C
  and Rust back-ends.
* Replaced the legacy `MACHEPS = 3e-39` with `DBL_EPSILON`, which is what the
  modified Cholesky step actually needs.
* Removed the unused backward-elimination / influence-function / `fopen`
  machinery from the C sources. Only the R-exposed engine remains.
* Added a compile-time parameter-count check (`MAXCOV = 30`) with a proper
  R-level error before any stack-allocated buffer is touched.
* Added `src/Makevars.win` so the package builds on Windows.
* Added a standard R-CMD-check GitHub Actions workflow covering Linux,
  macOS, and Windows.
