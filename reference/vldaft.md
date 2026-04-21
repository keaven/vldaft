# Fit an Accelerated Failure Time Model with Non-Proportional Hazards

Fits the Anderson (1991) AFT regression model where both location (mu)
and scale (sigma) are functions of covariates, and sigma can depend on
mu through a theta polynomial.

## Usage

``` r
vldaft(
  formula,
  data,
  dist = c("weibull", "logistic", "normal", "cauchy", "gamma"),
  theta = 0L,
  theta_vars = NULL,
  nu = 1,
  init = NULL,
  control = list(),
  adjust = TRUE,
  backend = c("c", "rust")
)
```

## Arguments

- formula:

  A formula of the form `Surv(time, status) ~ loc_vars | scale_vars`.
  Left of `|` specifies location (mu) covariates, right specifies scale
  (gamma) covariates. If no `|` is present, all covariates enter the
  location model and scale has intercept only.

- data:

  A data frame.

- dist:

  Character string specifying the error distribution: `"weibull"`
  (default), `"logistic"`, `"normal"`, `"cauchy"`, or `"gamma"`.

- theta:

  Integer, order of the theta polynomial linking log(sigma) to mu. 0 =
  no coupling (default), 1 = linear in mu, 2 = quadratic, etc.

- theta_vars:

  A one-sided formula specifying which location covariates form mu\*
  (the mu that feeds into the theta polynomial for sigma). Default
  `NULL` means all location covariates.

- nu:

  Numeric, shape parameter for the gamma distribution (only used when
  `dist = "gamma"`).

- init:

  Numeric vector of initial parameter values, or `NULL` for zeros.

- control:

  A list of control parameters. See Details.

- adjust:

  Logical, whether to center covariates (default `TRUE`).

- backend:

  Character, which computation backend to use: `"c"` (default) or
  `"rust"`.

## Value

An object of class `"vldaft"`.

## Details

The model assumes \$\$\log(T) = \mu + \sigma \cdot \varepsilon\$\$ where
the location is a linear function of covariates, \$\$\mu = \beta_0 +
\beta_1 x_1 + \cdots + \beta_p x_p\$\$ and the log-dispersion can depend
on covariates and on the location: \$\$\log(\sigma) = \gamma_0 +
\gamma_1 z_1 + \cdots + \gamma_q z_q + \theta_1 \mu^{\*} + \theta_2
\mu^{\*2} + \cdots\$\$ Here \\\mu^{\*} = \mu - \bar{\mu}\\ is the
centered location. When \\\theta_1 = 0\\, this reduces to the standard
linear location AFT model with proportional hazards under a Weibull
distribution. The error distribution \\F\\ can be Weibull (extreme
value), logistic, normal, Cauchy, or gamma.

Control parameters:

- acc:

  Convergence accuracy (default 0.0001)

- maxiter:

  Maximum Newton-Raphson iterations (default 50)

- maxhalv:

  Maximum step halvings (default 20)

## References

Anderson, K.M. (1991). A nonproportional hazards Weibull accelerated
failure time regression model. *Biometrics*, 47, 281-288.
