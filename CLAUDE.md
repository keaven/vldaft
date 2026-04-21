# vldaft — Project Instructions

## Package overview

vldaft implements the Anderson (1991) Varying Location and Dispersion
AFT model. The core computation is Newton-Raphson maximization of the
log-likelihood with modified Cholesky decomposition, implemented in both
C and Rust backends.

## Project structure

- `R/vldaft.R` — main fitting function
  [`vldaft()`](https://keaven.github.io/vldaft/reference/vldaft.md),
  formula parsing, dispatch to C or Rust
- `R/methods.R` — S3 methods: print, summary, coef, vcov, logLik
- `R/family_aft.R` — gamlss2 family constructors (AFT_Weibull, etc.)
- `src/*.c` — C backend (Newton-Raphson, distributions, matrix algebra,
  data I/O)
- `src/rust/src/` — Rust backend via extendr (same algorithms, f64
  precision)
- `src/init.c` — R routine registration for both C and Rust entry points
- `data/` — Framingham Heart Study datasets (.rda)
- `vignettes/anderson1991.Rmd` — reproduces all models from the 1991
  paper
- `tests/testthat/` — 92 unit tests covering paper results, backend
  parity, methods, edge cases

## Key conventions

- Log-likelihood is on the log-time scale (no Jacobian). To compare with
  `survreg`: `vldaft_ll + sum(d * log(t)) ≈ survreg_ll`
- Coefficient names use prefixes: `gamma:` (location), `eta:` (scale
  covariates), `theta:` (coupling)
- The `backend` argument selects `"c"` or `"rust"` via
  [`.Call()`](https://rdrr.io/r/base/CallExternal.html) dispatch
- C code uses float storage and MAXCOV=30, MAXSS=210 limits; Rust uses
  f64 throughout
- Packed symmetric (lower-triangle) matrix storage for the
  Hessian/covariance

## Build system

- `src/Makevars` compiles C sources and runs `cargo build --release` for
  the Rust static library
- Rust source is in `src/rust/` with `Cargo.toml`; uses extendr-api 0.7
  and libm 0.2
- System requirements: R ≥ 3.5, Rust toolchain (rustc + cargo)

## Testing

- Run tests: `devtools::test()` or `testthat::test_local()`
- All 92 tests should pass; C and Rust backends match to ~1e-7 precision
- Paper log-likelihoods: Model 1 = -2510.96, Model 2 = -2504.40, Model
  3b = -2496.65

## When modifying code

- If changing the C backend, mirror the change in the corresponding Rust
  module
- If changing the R interface arguments, update both the C call
  (`vldaft_fit`) and Rust call (`wrap__vldaft_fit_rust`) in `R/vldaft.R`
  and registration in `src/init.c`
- After changes, rebuild with `R CMD INSTALL .` and run
  `devtools::test()`
