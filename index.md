# vldaft

Varying Location and Dispersion Accelerated Failure Time (VLDAFT)
regression, implementing the model from Anderson (1991). Both the
location (μ) and dispersion (σ) of the log-survival time can depend on
covariates, and the dispersion can depend on the fitted location through
theta polynomial coupling — enabling nonproportional hazards under the
Weibull and other parametric distributions. C programming from
supporting the original 1991 manuscript was updated to standards
required for R. Either C or Rust backends can be used for computattions.

## Installation

``` r
# install.packages("pak")
pak::pak("keaven/vldaft")
```

### System requirements

A Rust toolchain (`rustc` + `cargo`) is required. Install from
<https://rustup.rs>.

## Usage

``` r
library(vldaft)
library(survival)
data(MLT50)

# Standard Weibull AFT (constant scale)
fit1 <- vldaft(Surv(YTCHD, CHD) ~ ln.Age. + ln.SBP. + ln.Chol. + MRW + Smoke,
               data = MLT50, dist = "weibull", theta = 0)

# Scale coupled to location via theta polynomial
fit2 <- vldaft(Surv(YTCHD, CHD) ~ ln.Age. + ln.SBP. + ln.Chol. + MRW + Smoke,
               data = MLT50, dist = "weibull", theta = 1)

summary(fit2)
```

## Features

- **Five error distributions**: Weibull, logistic, normal, Cauchy, gamma
- **Censoring**: right censoring, left censoring, left truncation
- **Scale–location coupling**: polynomial theta parameters linking σ to
  μ
- **Scale covariates**: separate covariates for the scale (dispersion)
  model via `theta_vars`
- **Dual backends**: C and Rust implementations selectable via
  `backend = c("c", "rust")`
- **gamlss2 families**:
  [`AFT_Weibull()`](https://keaven.github.io/vldaft/reference/AFT_families.md),
  [`AFT_Logistic()`](https://keaven.github.io/vldaft/reference/AFT_families.md),
  [`AFT_Normal()`](https://keaven.github.io/vldaft/reference/AFT_families.md),
  [`AFT_Cauchy()`](https://keaven.github.io/vldaft/reference/AFT_families.md),
  [`AFT_Gamma()`](https://keaven.github.io/vldaft/reference/AFT_families.md)
  for use with the gamlss2 framework
- **Datasets**: Framingham Heart Study subsets from Anderson (1991)

## Reference

Anderson, K. M. (1991). A nonproportional hazards Weibull accelerated
failure time regression model. *Biometrics*, 47(1), 281–288.
