# Slides

Presentation material for the VLDAFT methodology talk.

## Contents

- `vldaft-methodology.qmd` — 20-slide Quarto / revealjs deck covering
  the model, the C-to-Rust port, the `gamlss2` extension, and a small
  simulated-data demo of the dispersion-coupling concept and its use
  for risk prediction and RMST-based health-economic comparisons.

## Rendering

The deck is a [Quarto](https://quarto.org) document with embedded R
code chunks. To render:

```sh
quarto render slides/vldaft-methodology.qmd
```

This produces `slides/vldaft-methodology.html`. To preview live with
auto-reload:

```sh
quarto preview slides/vldaft-methodology.qmd
```

## R requirements

The conceptual figures (`fig-theta`, `fig-survival`) only need
`ggplot2`. The model-fitting and prediction chunks (marked with
`#| eval: false` for portability) call into the installed `vldaft`
package; install it once with

```r
pak::pak("keaven/vldaft")
```

and then either flip the chunks to `eval: true` or just discuss the
expected output as written. The numerical results in the talk text
were generated from the same simulated-data recipe used in
`tests/testthat/helper-simdata.R` so they are reproducible.

## Why a `slides/` folder, not a vignette?

The methodology vignette in `vignettes/anderson1991.Rmd` reproduces the
1991 Framingham analysis from sanitized result objects shipped under
`inst/extdata/`. The slides are aimed at a different audience (a
methodology talk rather than a how-to), live outside the package
build (`.Rbuildignore`), and intentionally pull from simulated data so
they require nothing beyond the public package.
