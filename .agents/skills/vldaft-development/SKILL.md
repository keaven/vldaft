---
name: vldaft-development
description: Draft workflow for general vldaft R/C/Rust package development, backend parity, compiled interface changes, tests, package checks, and maintenance tasks. Use for ordinary implementation and verification work; use vldaft-model-development instead when the task starts from a new statistical model prompt.
---

# vldaft Development

## Draft Status

This is a draft repo-local skill. Prefer the durable project guide in
`inst/site/llms.txt` when there is any conflict, and refine this workflow as
real maintenance tasks reveal better patterns.

## First Steps

Read `inst/site/llms.txt` before editing. Inspect `git status --short` and
preserve unrelated local changes.

Classify the task before changing files:

- R interface or S3 behavior.
- C backend behavior.
- Rust backend behavior.
- `.Call()` registration or argument contract.
- Tests only.
- Documentation only.

If the request starts from a new model idea, switch to
`$vldaft-model-development`.

## Implementation Rules

- Follow existing package style and naming.
- Keep C and Rust implementations in sync for standard non-cure AFT behavior.
- Mirror C partials-engine changes in `src/rust/src/model.rs` when applicable.
- Update `.Call()` wrappers and `src/init.c` together when arguments change.
- Keep coefficient names stable: `gamma:`, `eta:`, `cure:`, and `theta1`,
  `theta2`, ...
- Preserve the log-likelihood scale conventions documented in
  `inst/site/llms.txt`.
- Do not mix unrelated refactors into feature or bug-fix patches.

## File Routing

- R wrapper/API: `R/vldaft.R`.
- Methods and prediction: `R/methods.R`.
- `gamlss2` families: `R/family_aft.R`.
- C standard backend: `src/loscreg8.c`, `src/vldaft_fit.c`, `src/dataio8.c`.
- C cure evaluator: `src/vldaft_cure_eval.c`.
- Rust backend: `src/rust/src/lib.rs`, `src/rust/src/model.rs`.
- Registration: `src/init.c`.
- Compile-time parameter caps: `src/rj8def.h` and Rust limits together.
- Tests: `tests/testthat/`.

## Validation Ladder

Start narrow, then broaden:

```r
testthat::test_file("tests/testthat/test-vldaft-basic.R")
testthat::test_file("tests/testthat/test-vldaft-cure.R")
devtools::test()
```

For compiled-code or exported-interface changes, run:

```sh
R CMD INSTALL .
```

For docs, vignettes, README, or pkgdown navigation changes, run:

```r
pkgdown::build_site(devel = TRUE, preview = FALSE)
```

Report exactly which checks ran and any residual risk.
