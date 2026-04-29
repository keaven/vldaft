---
name: vldaft-model-development
description: Use when turning a statistical survival-model prompt into vldaft package code, especially new AFT distributions, cure models, censoring or truncation behavior, parameter/link definitions, likelihood/score implementations, C/Rust backend changes, gamlss2 family support, simulation checks, backend parity tests, and model-facing documentation.
---

# vldaft Model Development

## Overview

Use this skill to convert a model idea into a validated `vldaft`
implementation. Start by pinning down the statistical contract, then route the
change through the R API, C and Rust back ends, tests, and docs as needed.

## Required First Steps

Read `inst/site/llms.txt` before editing. Inspect `git status --short` and
preserve unrelated local changes.

For model prompts, produce a short implementation plan before coding unless the
user explicitly asks for analysis only. Include:

- Target model and user-facing API.
- Parameter blocks, links, constraints, and coefficient naming.
- Supported censoring/truncation cases.
- Required implementation surfaces: R only, C, Rust, cure evaluator,
  `gamlss2`, methods, docs, tests.
- Validation strategy.

If the prompt leaves the likelihood, parameterization, or censoring semantics
ambiguous, ask the smallest clarifying question before implementation.

## Model Contract Checklist

Before touching compiled code, write down the contract in local notes or the
turn summary:

- Response scale: time scale, log-time scale, or mixture over a baseline.
- Linear predictors: `mu`, `sigma`, `nu`/cure, theta polynomial terms, or new
  parameter blocks.
- Links and bounds: identity, log, logit, fixed values, shape constraints.
- Likelihood contributions for right, left, interval, and counting-process
  inputs when supported.
- Left-truncation adjustment.
- Prediction outputs and whether they use raw or internally centered
  covariates.
- Log-likelihood scale and Jacobian convention.

Do not bury a new convention in code only; expose it in docs and tests.

## Implementation Routing

Use these routing rules:

- R formula/API changes: update `R/vldaft.R`, methods in `R/methods.R`, tests,
  roxygen/Rd, and README or vignette examples when user-facing.
- Standard non-cure likelihood changes: update both `src/loscreg8.c` and
  `src/rust/src/model.rs`, plus the C/Rust wrappers if call arguments change.
- `.Call()` contract changes: update `src/vldaft_fit.c`, `src/init.c`,
  `src/rust/src/lib.rs`, and the R call site together.
- Cure-model changes: update the R wrapper and `src/vldaft_cure_eval.c`; add
  or adjust `gamlss2` family behavior in `R/family_aft.R` when relevant.
- New parameter blocks: update coefficient naming, summary/print methods,
  `coef()`, `vcov()`, `predict()`, and documentation.
- Compile-time limits: keep C `MAXCOV` and Rust limits in sync.

Keep implementation/tests separate from documentation-only commits when the
diff is large.

## Numerical Validation

Prefer checks that prove the model, not just the plumbing:

- Simulate data from the target model and verify recovery or likelihood
  improvement over a nested model.
- Compare analytic scores with finite differences where practical.
- For standard AFT back-end changes, test C/Rust parity to tight tolerance.
- For `gamlss2` support, compare log-likelihoods and linear predictors against
  standalone `vldaft()` on the same simulated data.
- Exercise censoring cases touched by the change: right, left, interval, and
  left truncation as relevant.

## Validation Ladder

Run the narrowest meaningful checks first, then broaden:

```r
testthat::test_file("tests/testthat/test-vldaft-basic.R")
testthat::test_file("tests/testthat/test-vldaft-cure.R")
devtools::test()
```

For compiled or exported API changes, also run:

```sh
R CMD INSTALL .
```

For docs, vignettes, or pkgdown navigation changes, run:

```r
pkgdown::build_site(devel = TRUE, preview = FALSE)
```

Report any skipped validation and why.
