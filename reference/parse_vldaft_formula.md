# Parse vldaft formula

Splits a formula like Surv(time, status) ~ x1 + x2 \| x3 + x4 into
location and scale parts.

## Usage

``` r
parse_vldaft_formula(formula, data)
```

## Arguments

- formula:

  The formula.

- data:

  The data frame.

## Value

A list with response, loc_formula, scale_formula.
