// newton_raphson.rs — Port of nrmod.c

use crate::matrix::{cholmod, chol, triinv, trimult, lbak, ubak, symmult, mmult1};

const MAXCOV: usize = 30;
const MAXSS: usize = 210;

/// Newton-Raphson optimization with modified Cholesky decomposition.
///
/// `par`: function computing partials (m, beta, db, ddb, ll) -> result code
/// Returns number of iterations on success, negative on failure.
pub fn nrsolmod<F>(
    par: &mut F,
    m: usize,
    finit: &mut f64,
    fmax: &mut f64,
    score: &mut f64,
    beta: &mut [f64],
    cov: &mut [f64],
    acc: f64,
    maxiter: usize,
    wk: &mut [f64],
    delb: &mut [f64],
    db: &mut [f64],
    ddb: &mut [f64],
    maxhalv: usize,
) -> i32
where
    F: FnMut(usize, &[f64], &mut [f64], &mut [f64], &mut f64) -> i32,
{
    let mut steep = [0.0f64; MAXCOV];
    let mut steepdel = [0.0f64; MAXCOV];

    // Compute initial partials
    let mut i = par(m, beta, db, ddb, finit);
    if i != 0 {
        return i;
    }
    *fmax = *finit;

    // Compute steepest ascent variables
    for k in 0..m {
        steep[k] = beta[k] + db[k];
        steepdel[k] = db[k];
    }

    // First NR step
    let mut conv = nr_itmod(m, beta, db, ddb, delb, cov, acc, wk);

    // Compute score statistic
    if conv == 2 {
        *score = -1.0;
    } else {
        symmult(cov, db, wk, m, 1);
        mmult1(wk, db, &mut [0.0; 1], 1, m, 1);
        // Actually we need to write to score
        let mut score_buf = [0.0f64; 1];
        mmult1(wk, db, &mut score_buf, 1, m, 1);
        *score = score_buf[0];
    }

    // Complete iterations
    let mut iter = 1usize;
    while iter < maxiter && (conv == 0 || conv == 2) {
        let fold = *fmax;

        // Compute partials
        conv = par(m, beta, db, ddb, fmax);
        if conv == -2 { return -2; }
        let ftem = fold - *fmax;

        // Step halving if objective decreases or partials error
        for _j in 0..maxhalv {
            if fold <= *fmax && conv >= 0 { break; }

            for k in 0..m {
                delb[k] /= 2.0;
                beta[k] -= delb[k];
            }
            conv = par(m, beta, db, ddb, fmax);
            if conv == -2 { return -2; }

            // Try steepest ascent
            if fold > *fmax || conv < 0 {
                conv = par(m, &steep[..m], db, ddb, fmax);
                if conv == -2 { return -2; }
                if fold > *fmax || conv < 0 {
                    for k in 0..m {
                        steepdel[k] /= 2.0;
                        steep[k] -= steepdel[k];
                    }
                } else {
                    for k in 0..m {
                        beta[k] = steep[k];
                    }
                }
            }
        }

        if fold > *fmax {
            wk[0] = ftem;
            return -6;
        }
        if conv < 0 {
            *fmax = 0.0;
            let mut j = 0;
            for k in 0..m {
                db[k] = 0.0;
                for _ii in 0..=k {
                    ddb[j] = 0.0;
                    j += 1;
                }
            }
            return -7;
        }

        // Steepest ascent variables
        for k in 0..m {
            steep[k] = beta[k] + db[k];
            steepdel[k] = db[k];
        }

        // Next NR step
        conv = nr_itmod(m, beta, db, ddb, delb, cov, acc, wk);
        iter += 1;
    }

    if conv == 0 || conv == 2 {
        return -4;
    }

    // Final partials, inverse
    if par(m, beta, db, ddb, fmax) != 0 {
        return -3;
    }
    if chol(ddb, cov, m) != 0 {
        return -5;
    }
    let size = m * (m + 1) / 2;
    let mut temp = vec![0.0f64; size];
    temp[..size].copy_from_slice(&cov[..size]);
    triinv(&temp, cov, wk, m);
    temp[..size].copy_from_slice(&cov[..size]);
    trimult(&temp, cov, m);

    iter as i32
}

/// Single Newton-Raphson iteration with modified Cholesky
fn nr_itmod(
    n: usize,
    b: &mut [f64],
    db: &mut [f64],
    ddb: &mut [f64],
    delb: &mut [f64],
    ddbinv: &mut [f64],
    acc: f64,
    wk: &mut [f64],
) -> i32 {
    let mut conv: i32;

    // Modified Cholesky
    if cholmod(ddb, ddbinv, delb, n) == 1 {
        conv = 3;
    } else {
        conv = 1;
    }

    // Solve for change in b
    lbak(ddbinv, wk, db, 1, n);
    for j in 0..n {
        wk[j] /= delb[j];
    }
    ubak(ddbinv, delb, wk, 1, n);

    // Find inverse of ddb
    let size = n * (n + 1) / 2;
    let mut temp = vec![0.0f64; size];
    temp[..size].copy_from_slice(&ddbinv[..size]);
    triinv(&temp, ddbinv, wk, n);
    temp[..size].copy_from_slice(&ddbinv[..size]);
    trimult(&temp, ddbinv, n);

    // Update b and check convergence
    let mut j: i64 = -1;
    for k in 0..n {
        j += (k as i64) + 1;
        b[k] += delb[k];
        if conv != 0 && conv != 2 {
            let a = (ddbinv[j as usize]).sqrt();
            let a = delb[k] / a;
            if a > acc || a < -acc {
                conv -= 1;
            }
        }
    }
    conv
}
