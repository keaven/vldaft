//! Model state for the Anderson (1991) AFT engine.
//!
//! Mirrors dataio8.c and loscreg8.c::ac_tpar8 but avoids the C globals by
//! encapsulating everything in `ModelState`.

use crate::distributions;

/// Which error distribution to use. Matches the C `dist_choice` switch.
pub enum Distribution {
    Weibull,
    Logistic,
    Normal,
    Cauchy,
    Gamma { _nu: f64, enu: f64 },
}

impl Distribution {
    pub fn call(&self, event: i32, x: f64, f: &mut [f64; 3]) -> Result<(), i32> {
        match self {
            Distribution::Weibull        => distributions::weibull(event, x, f),
            Distribution::Logistic       => distributions::logistic(event, x, f),
            Distribution::Normal         => distributions::normal_dist(event, x, f),
            Distribution::Cauchy         => distributions::cauchy_dist(event, x, f),
            Distribution::Gamma { enu, .. } => {
                distributions::gamma_dist(event, x, f, *enu)
            }
        }
    }
}

/// All per-fit state, replacing the C globals.
pub struct ModelState {
    pub m: usize,
    pub n: usize,
    pub scale: usize,
    pub mlo: usize,
    pub mlo1: usize,
    pub ntheta: usize,
    pub npar: usize,

    pub tcol: usize,
    pub t2col: i32,
    pub cencol: usize,
    pub startcol: i32,
    pub t0col: usize,

    pub pscale: Vec<usize>,
    pub ploc: Vec<usize>,

    pub reoff: Vec<i32>,
    pub rgtcen: i32,

    pub x: Vec<f64>,        // row-major copy of the R matrix
    pub xbar: Vec<f64>,     // column means
    pub madj: bool,

    pub dist: Distribution,

    /// Scratch record used by `ac_tpar8` so we don't reallocate per row.
    scratch: Vec<f64>,
}

impl ModelState {
    /// Build from an R column-major numeric matrix.
    pub fn from_r_data(
        rdata: &[f64],
        nobs: usize,
        ncol: usize,
        loc_cols: &[i32],
        scale_cols: &[i32],
        ntheta: usize,
        mlo1_override: i32,
        time_col: usize,
        time2_col: i32,
        event_col: usize,
        event_vals: &[i32],
        rgtcen_val: i32,
        start_col: i32,
        t0_col: i32,
        dist: Distribution,
        adjust: bool,
    ) -> Self {
        let nlo = loc_cols.len();
        let nsc = scale_cols.len();

        // Column-major -> row-major double copy and column means.
        let mut x = vec![0.0f64; nobs * ncol];
        let mut xbar = vec![0.0f64; ncol];
        for i in 0..nobs {
            for j in 0..ncol {
                let val = rdata[j * nobs + i];
                x[i * ncol + j] = val;
                xbar[j] += val;
            }
        }
        for j in 0..ncol { xbar[j] /= nobs as f64; }

        let spmlo = nsc + nlo;
        let mlo1 = if ntheta == 0 {
            nlo
        } else if mlo1_override >= 0 {
            mlo1_override as usize
        } else {
            1
        };

        let ploc: Vec<usize>   = loc_cols.iter().map(|&c| c as usize).collect();
        let pscale: Vec<usize> = scale_cols.iter().map(|&c| c as usize).collect();

        ModelState {
            m: ncol,
            n: nobs,
            scale: nsc,
            mlo: nlo,
            mlo1,
            ntheta,
            npar: nsc + nlo + ntheta,
            tcol: time_col,
            t2col: time2_col,
            cencol: event_col,
            startcol: start_col,
            t0col: if t0_col >= 0 { t0_col as usize } else { 0 },
            pscale,
            ploc,
            reoff: event_vals.to_vec(),
            rgtcen: rgtcen_val,
            x,
            xbar,
            madj: adjust,
            dist,
            scratch: vec![0.0f64; spmlo + 5],
        }
    }

    /// Fill `self.scratch` with the record for observation `obs`. Layout
    /// matches getrec() in dataio8.c.
    fn fill_record(&mut self, obs: usize) {
        let scale = self.scale;
        let mlo = self.mlo;
        let spml = scale + mlo;
        let start = obs * self.m;
        let end = start + self.m;
        let xp = &self.x[start..end];
        let y = &mut self.scratch;

        if self.madj {
            if scale > 0 { y[0] = xp[self.pscale[0]]; }
            for i in 1..scale {
                y[i] = xp[self.pscale[i]] - self.xbar[self.pscale[i]];
            }
            if mlo > 0 { y[scale] = xp[self.ploc[0]]; }
            for i in 1..mlo {
                y[scale + i] = xp[self.ploc[i]] - self.xbar[self.ploc[i]];
            }
        } else {
            for i in 0..scale { y[i] = xp[self.pscale[i]]; }
            for i in 0..mlo   { y[i + scale] = xp[self.ploc[i]]; }
        }

        y[spml] = xp[self.tcol].ln();
        y[spml + 1] = 0.0;
        for ev in &self.reoff {
            if xp[self.cencol] as i32 == *ev { y[spml + 1] = 1.0; }
        }
        if xp[self.cencol] as i32 == self.rgtcen { y[spml + 1] = -1.0; }
        if xp[self.cencol] as i32 == 2 { y[spml + 1] = 2.0; }

        if self.startcol >= 0 && xp[self.t0col] > 0.0 {
            y[spml + 2] = xp[self.startcol as usize];
            if y[spml + 2] == 1.0 { y[spml + 3] = xp[self.t0col].ln(); }
        } else {
            y[spml + 2] = 0.0;
        }
        y[spml + 4] = if self.t2col >= 0 && y[spml + 1] == 2.0 {
            xp[self.t2col as usize].ln()
        } else {
            0.0
        };
    }

    /// Compute log-likelihood, gradient, and packed Hessian at `beta`.
    /// Port of ac_tpar8() in loscreg8.c.
    pub fn ac_tpar8(
        &mut self,
        m: usize,
        beta: &[f64],
        db: &mut [f64],
        ddb: &mut [f64],
        ll1: &mut f64,
    ) -> i32 {
        let scale = self.scale;
        let mlo = self.mlo;
        let ntheta = self.ntheta;
        let spmlo = scale + mlo;
        let spmlo1 = if spmlo == m { spmlo } else { scale + self.mlo1 };
        let n_obs = self.n;

        *ll1 = 0.0;
        let mut j = 0usize;
        for i in 0..m {
            db[i] = 0.0;
            for _ii in 0..=i { ddb[j] = 0.0; j += 1; }
        }

        for obs in 0..n_obs {
            self.fill_record(obs);
            let x = &self.scratch;

            let mut mustar = 0.0f64;
            for i in spmlo1..spmlo { mustar += beta[i] * x[i]; }
            let mut mu = mustar;
            for i in scale..spmlo1 { mu += beta[i] * x[i]; }

            let event: i32 = x[spmlo + 1] as i32;
            let start: i32 = if x[spmlo + 2] == 1.0 { 1 } else { 0 };

            let mut log_sigma = 0.0f64;
            for i in 0..scale { log_sigma += beta[i] * x[i]; }
            let mut dpmu = 0.0f64;
            let mut ddpmu = 0.0f64;
            if ntheta > 0 {
                let mut pow_m = 1.0f64;
                dpmu = beta[spmlo];
                log_sigma += beta[spmlo] * mustar;
                for i in 1..ntheta {
                    let prev_pow = pow_m;
                    pow_m *= mustar;
                    dpmu += beta[spmlo + i] * (i + 1) as f64 * pow_m;
                    ddpmu += beta[spmlo + i] * ((i + 1) * i) as f64 * prev_pow;
                    log_sigma += beta[spmlo + i] * pow_m * mustar;
                }
            }
            if !(-39.0..=39.0).contains(&log_sigma) { return -1; }
            let sigma = log_sigma.exp();

            let mut mu1 = vec![0.0f64; m];
            let mut v1 = vec![0.0f64; m];
            let mut v2 = vec![vec![0.0f64; m]; m];
            for i in 0..scale { v1[i] = x[i]; }
            for i in scale..spmlo {
                mu1[i] = x[i];
                if ntheta > 0 && i >= spmlo1 {
                    v1[i] = dpmu * x[i];
                }
            }
            for k in 0..ntheta {
                v1[spmlo + k] = mustar.powi((k + 1) as i32);
            }
            if ntheta > 0 && spmlo > spmlo1 {
                for i in spmlo1..spmlo {
                    for k in spmlo1..spmlo {
                        v2[i][k] = ddpmu * x[i] * x[k];
                    }
                }
                for i in spmlo1..spmlo {
                    let mut pow_m = 1.0f64;
                    for k in 0..ntheta {
                        let val = (k + 1) as f64 * pow_m * x[i];
                        v2[i][spmlo + k] = val;
                        v2[spmlo + k][i] = val;
                        pow_m *= mustar;
                    }
                }
            }

            let residual_partials = |uval: f64| -> (Vec<f64>, Vec<Vec<f64>>) {
                let mut u1 = vec![0.0f64; m];
                let mut u2 = vec![vec![0.0f64; m]; m];
                for i in 0..m {
                    u1[i] = -mu1[i] / sigma - uval * v1[i];
                    for k in 0..m {
                        u2[i][k] = (mu1[i] * v1[k] + mu1[k] * v1[i]) / sigma
                            + uval * v1[i] * v1[k] - uval * v2[i][k];
                    }
                }
                (u1, u2)
            };

            let u_val = (x[spmlo] - mu) / sigma;
            let (u1, u2) = residual_partials(u_val);
            let mut f = [0.0f64; 3];
            let mut f_hi = [0.0f64; 3];
            let mut f0 = [0.0f64; 3];
            let q1: f64;
            let q2: f64;
            let q11: f64;
            let q22: f64;
            let q12: f64;
            let explicit_v: f64;
            let mut u_hi1 = vec![0.0f64; m];
            let mut u_hi2 = vec![vec![0.0f64; m]; m];

            if event == 2 {
                if x[spmlo + 4] <= x[spmlo] { return -1; }
                let u_hi = (x[spmlo + 4] - mu) / sigma;
                let hi_partials = residual_partials(u_hi);
                u_hi1 = hi_partials.0;
                u_hi2 = hi_partials.1;
                if self.dist.call(0, u_val, &mut f).is_err()
                    || self.dist.call(0, u_hi, &mut f_hi).is_err() {
                    return -1;
                }
                let ratio = (f_hi[0] - f[0]).exp();
                if ratio >= 1.0 { return -1; }
                *ll1 += f[0] + (1.0 - ratio).ln();
                let a = f[0].exp();
                let b = f_hi[0].exp();
                let den = a - b;
                if den <= 0.0 { return -1; }
                q1 = a * f[1] / den;
                q2 = -b * f_hi[1] / den;
                q11 = a * (f[2] + f[1] * f[1]) / den - q1 * q1;
                q22 = -b * (f_hi[2] + f_hi[1] * f_hi[1]) / den - q2 * q2;
                q12 = a * f[1] * b * f_hi[1] / (den * den);
                explicit_v = 0.0;
            } else {
                if self.dist.call(event, u_val, &mut f).is_err() {
                    return -1;
                }
                explicit_v = if event == 1 { -1.0 } else { 0.0 };
                *ll1 += f[0] + explicit_v * log_sigma;
                q1 = f[1];
                q2 = 0.0;
                q11 = f[2];
                q22 = 0.0;
                q12 = 0.0;
            }

            let mut u_start1 = vec![0.0f64; m];
            let mut u_start2 = vec![vec![0.0f64; m]; m];
            if start == 1 {
                let u0 = (x[spmlo + 3] - mu) / sigma;
                let start_partials = residual_partials(u0);
                u_start1 = start_partials.0;
                u_start2 = start_partials.1;
                if self.dist.call(0, u0, &mut f0).is_err() { return -1; }
                *ll1 -= f0[0];
            }

            j = 0;
            for i in 0..m {
                let mut grad = q1 * u1[i] + explicit_v * v1[i];
                if event == 2 {
                    grad += q2 * u_hi1[i];
                }
                if start == 1 {
                    grad -= f0[1] * u_start1[i];
                }
                db[i] += grad;

                for k in 0..=i {
                    let mut second = q11 * u1[i] * u1[k] + q1 * u2[i][k]
                        + explicit_v * v2[i][k];
                    if event == 2 {
                        second += q22 * u_hi1[i] * u_hi1[k] + q2 * u_hi2[i][k]
                            + q12 * (u1[i] * u_hi1[k] + u_hi1[i] * u1[k]);
                    }
                    if start == 1 {
                        second -= f0[2] * u_start1[i] * u_start1[k]
                            + f0[1] * u_start2[i][k];
                    }
                    ddb[j] -= second;
                    j += 1;
                }
            }
        }
        0
    }
}
