//! extendr entry point for the vldaft Rust backend.

use extendr_api::prelude::*;

mod distributions;
mod matrix;
mod model;
mod newton_raphson;

use model::{Distribution, ModelState};

/// Fit the AFT regression model (Rust backend).
/// The interface mirrors the C `vldaft_fit` entry point one-for-one.
#[extendr]
fn vldaft_fit_rust(
    data: RMatrix<f64>,
    time_col: i32,
    time2_col: i32,
    event_col: i32,
    event_vals: Integers,
    rgtcen_val: i32,
    start_col: i32,
    t0_col: i32,
    loc_cols: Integers,
    scale_cols: Integers,
    ntheta: i32,
    mlo1: i32,
    dist: i32,
    nu: f64,
    init: Nullable<Doubles>,
    acc: f64,
    maxiter: i32,
    maxhalv: i32,
    adjust: i32,
) -> List {
    let nobs = data.nrows();
    let ncol = data.ncols();
    let rdata = data.as_real_slice().unwrap();

    let distribution = match dist {
        1 => Distribution::Weibull,
        2 => Distribution::Logistic,
        3 => Distribution::Normal,
        4 => Distribution::Cauchy,
        5 => Distribution::Gamma { _nu: nu, enu: nu.exp() },
        _ => panic!("Unknown distribution: {}", dist),
    };

    let ev_vals: Vec<i32> = event_vals.iter().map(|v| v.inner()).collect();
    let loc_c:   Vec<i32> = loc_cols.iter().map(|v| v.inner()).collect();
    let scale_c: Vec<i32> = scale_cols.iter().map(|v| v.inner()).collect();

    let mut state = ModelState::from_r_data(
        rdata, nobs, ncol,
        &loc_c, &scale_c,
        ntheta as usize, mlo1,
        time_col as usize, time2_col, event_col as usize,
        &ev_vals, rgtcen_val,
        start_col, t0_col,
        distribution,
        adjust != 0,
    );

    let npar = state.npar;
    if npar > 30 {
        panic!(
            "vldaft (rust backend): too many parameters ({} > 30). \
             Reduce the model or rebuild with a larger compile-time cap.",
            npar
        );
    }

    // MAXCOV = 30 matches the C backend; keep the two in sync.
    let maxcov = 30usize;
    let maxss  = maxcov * (maxcov + 1) / 2;

    let mut beta = vec![0.0f64; maxcov];
    if let NotNull(init_vals) = init {
        let ninit = init_vals.len().min(npar);
        for i in 0..ninit { beta[i] = init_vals[i].inner(); }
    }

    let mut cov  = vec![0.0f64; maxss];
    let mut wk   = vec![0.0f64; maxcov];
    let mut delb = vec![0.0f64; maxcov];
    let mut db   = vec![0.0f64; maxcov];
    let mut ddb  = vec![0.0f64; maxss];
    let mut finit = 0.0f64;
    let mut fmax  = 0.0f64;
    let mut score_stat = 0.0f64;

    let iter = {
        let mut par = |m: usize,
                       beta_in: &[f64],
                       db_out: &mut [f64],
                       ddb_out: &mut [f64],
                       ll: &mut f64| -> i32 {
            state.ac_tpar8(m, beta_in, db_out, ddb_out, ll)
        };
        newton_raphson::nrsolmod(
            &mut par,
            npar,
            &mut finit,
            &mut fmax,
            &mut score_stat,
            &mut beta,
            &mut cov,
            acc,
            maxiter as usize,
            &mut wk,
            &mut delb,
            &mut db,
            &mut ddb,
            maxhalv as usize,
        )
    };

    // Expand packed lower-triangular cov into a symmetric dense matrix.
    let beta_r: Vec<f64> = beta[..npar].to_vec();
    let mut cov_full = vec![0.0f64; npar * npar];
    let mut k = 0usize;
    for i in 0..npar {
        for j in 0..=i {
            cov_full[i + j * npar] = cov[k];
            cov_full[j + i * npar] = cov[k];
            k += 1;
        }
    }
    let cov_matrix = RMatrix::new_matrix(npar, npar, |r, c| cov_full[r + c * npar]);

    list!(
        coefficients = beta_r,
        vcov         = cov_matrix,
        loglik       = fmax,
        loglik_init  = finit,
        score        = score_stat,
        iter         = iter,
        npar         = npar as i32,
        nobs         = nobs as i32
    )
}

extendr_module! {
    mod vldaft;
    fn vldaft_fit_rust;
}
