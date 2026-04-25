# vldaft

<!-- badges: start -->
[![R-CMD-check](https://github.com/keaven/vldaft/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/keaven/vldaft/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Varying Location and Dispersion Accelerated Failure Time (VLDAFT)
regression, implementing the model from Anderson (1991) and supplementing
it with an optional Poisson mixture cure model for survival distributions
with an apparent plateau. Both the
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
<https://rustup.rs>. Source installs run `cargo build`; unless the
crate cache is already populated, Cargo needs network access to
download Rust dependencies.

## Usage

The Framingham Heart Study data used in Anderson (1991) are
confidential and are **not** shipped with this package. The example below
uses the freely available `lung` dataset from the `survival` package so
that anyone can run it out of the box.

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
lr_theta <- 2 * (fit2$loglik - fit1$loglik)
lr_theta
pchisq(lr_theta, df = 1, lower.tail = FALSE)

# Intercept-only Poisson mixture cure model
fit3 <- vldaft(
  Surv(time, status) ~ age + sex + ph.ecog,
  data  = na.omit(lung),
  dist  = "weibull",
  theta = 0,
  cure  = TRUE,
  control = list(maxiter = 1000L, acc = 1e-8)
)

summary(fit3)

# Likelihood-ratio test for adding the cure fraction
# Cure likelihoods are evaluated on the time scale, so compare fit1 there too.
lr_cure <- 2 * (fit3$loglik - as.numeric(logLik(fit1, jacobian = TRUE)))
lr_cure
pchisq(lr_cure, df = 1, lower.tail = FALSE)
```

## Features

- **Five error distributions**: Weibull, logistic, normal, Cauchy, gamma.
- **Poisson-mixture cure families**: generic cure wrappers over the AFT
  baselines, including `Cure_Weibull()` and `Cure_Exponential()`.
- **Censoring**: right censoring, left censoring, interval censoring,
  and left truncation.
- **Scale–location coupling**: polynomial `theta` parameters linking
  `log(sigma)` to `mu`.
- **Scale covariates**: separate covariates for the scale (dispersion)
  model via formula syntax `Surv(...) ~ loc_vars | scale_vars`.
- **Dual back ends**: C and Rust implementations, selectable via
  `backend = c("c", "rust")`, including interval-censored fits.
- **`gamlss2` families**: `AFT_Weibull()`, `AFT_Logistic()`,
  `AFT_Normal()`, `AFT_Cauchy()`, `AFT_Gamma()`, plus
  `Cure_Weibull()`, `Cure_Logistic()`, `Cure_Normal()`,
  `Cure_Cauchy()`, `Cure_Gamma()`, and `Cure_Exponential()`.

## Reference

Anderson, K. M. (1991). A nonproportional hazards Weibull accelerated
failure time regression model. *Biometrics*, 47(1), 281–288.
