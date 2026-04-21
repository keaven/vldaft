# Methods for vldaft Objects

Print, summary, coefficient extraction, variance-covariance matrix, and
log-likelihood for `vldaft` model fits.

## Usage

``` r
# S3 method for class 'vldaft'
print(x, ...)

# S3 method for class 'vldaft'
summary(object, ...)

# S3 method for class 'summary.vldaft'
print(x, ...)

# S3 method for class 'vldaft'
coef(object, ...)

# S3 method for class 'vldaft'
vcov(object, ...)

# S3 method for class 'vldaft'
logLik(object, ...)
```

## Arguments

- x, object:

  A `vldaft` object.

- ...:

  Additional arguments (ignored).

## Value

`print` returns `x` invisibly. `summary` returns a `summary.vldaft`
object. `coef` returns the coefficient vector. `vcov` returns the
variance-covariance matrix. `logLik` returns a `logLik` object with `df`
and `nobs` attributes.
