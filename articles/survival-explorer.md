# Survival Distribution Explorer

    ## <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@xiee/utils@1.14.35/css/default.min.css">
    ## <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@xiee/utils@1.14.35/css/tabsets.min.css">
    ## <script src="https://cdn.jsdelivr.net/npm/@xiee/utils@1.14.35/js/tabsets.min.js" defer></script>
    ## <script src="https://cdn.jsdelivr.net/npm/chart.js@4.5.1" defer></script>

## Overview

This interactive vignette lets you explore the five survival
distributions in `vldaft` on the accelerated failure time (AFT) log-time
scale. You can compare families matched on the 25th and 75th
percentiles, or compare multiple parameter sets within a single family.
Toggle between CDF, PDF, and hazard rate.

**AFT parameterization:** $w = \left( \log t - \mu \right)/\sigma$. All
curves below are computed from the standardized error distribution on
$w$, then mapped to time $t$.

## Explorer

- Compare families
- Compare parameters

Distributions

Weibull

Logistic

Normal

Cauchy

Gamma

Match quantiles t25 t75

All selected distributions match these percentiles.

Axis limits x max y max

Reset limits

Display CDF PDF Hazard rate

Gamma shape nu

Shape = exp(nu)

Family Weibull Logistic Normal Cauchy Gamma

Display CDF PDF Hazard rate

Rows

Add row

Remove row

Axis limits x max y max

Reset limits

| Label | mu  | sigma | nu  |
|-------|-----|-------|-----|
