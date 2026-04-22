# vldaft

<!-- badges: start -->
[![R-CMD-check](https://github.com/keaven/vldaft/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/keaven/vldaft/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Varying Location and Dispersion Accelerated Failure Time (VLDAFT)
regression, implementing the model from Anderson (1991). Both the
location (μ) and the dispersion (σ) of the log-survival time can depend
on covariates, and the dispersion can depend on the fitted location
through a polynomial coupling — enabling non-proportional hazards
under the Weibull and other parametric distributions. The C code that
supported the original 1991 manuscript was updated to the standards
required for an R package. Either a C or a Rust back end can be used
for the computations.

## Installation

``` r
# install.packages("pak")
pak::pak("keaven/vldaft")
```

### System requirements

A Rust toolchain (`rustc` + `cargo`) is required. Install from
<https://rustup.rs>.

## Usage

The Framingham Heart Study data used in Anderson (1991) are
confidential and are **not** shipped with this package. They live in
the companion private repository
[`keaven/vldaft.data`](https://github.com/keaven/vldaft.data),
which also carries the regression tests that depend on them.

The example below uses the freely available `lung` dataset from
the `survival` package so that anyone can run it out of the box.

``` r
library(vldaft)
library(survival)

# Constant-dispersion (standard Weibull AFT)
fit1 <- vldaft(
  Surv(time, status) ~ age + sex + ph.ecog,
  data  = na.omit(lung),
  dist  = "weibull",
  theta = 0
)

# Dispersion coupled to location via a linear theta polynomial
fit2 <- vldaft(
  Surv(time, status) ~ age + sex + ph.ecog,
  data  = na.omit(lung),
  dist  = "weibull",
  theta = 1
)

summary(fit2)

# Likelihood-ratio test for "is dispersion proportional to location?"
2 * (fit2$loglik - fit1$loglik)
```

## Features

- **Five error distributions**: Weibull, logistic, normal, Cauchy, gamma.
- **Censoring**: right censoring, left censoring, left truncation.
- **Scale–location coupling**: polynomial `theta` parameters linking
  `log(sigma)` to `mu`.
- **Scale covariates**: separate covariates for the scale (dispersion)
  model via formula syntax `Surv(...) ~ loc_vars | scale_vars`.
- **Dual back ends**: C and Rust implementations, selectable via
  `backend = c("c", "rust")`.
- **`gamlss2` families**: `AFT_Weibull()`, `AFT_Logistic()`,
  `AFT_Normal()`, `AFT_Cauchy()`, `AFT_Gamma()`.

## Reference

Anderson, K. M. (1991). A nonproportional hazards Weibull accelerated
failure time regression model. *Biometrics*, 47(1), 281–288.
