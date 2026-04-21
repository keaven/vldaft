# AFT Family Constructors for gamlss2

These functions create gamlss2 family objects for accelerated failure
time models with five error distributions. Each family handles
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) responses with
right censoring, left censoring, and left truncation.

## Usage

``` r
AFT_Weibull(theta = 0L)

AFT_Logistic(theta = 0L)

AFT_Normal(theta = 0L)

AFT_Cauchy(theta = 0L)

AFT_Gamma(theta = 0L, nu = 1)
```

## Arguments

- theta:

  Integer, order of the theta polynomial linking log(sigma) to mu.
  Default 0 (no coupling).

## Value

A `gamlss2.family` object.

## Details

The model is: \$\$\log(T) = \mu + \sigma\_{eff} \cdot \varepsilon\$\$
where \\\log(\sigma\_{eff}) = \log(\sigma) + \sum_k \theta_k (\mu -
\bar\mu)^k\\

Parameters:

- mu:

  Location parameter (identity link)

- sigma:

  Scale parameter (log link)

- theta1, theta2, ...:

  Coupling parameters (identity link), if theta \> 0
