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
    Gamma { nu: f64, enu: f64 },
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
            scratch: vec![0.0f64; spmlo + 4],
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

        if self.startcol >= 0 && xp[self.t0col] > 0.0 {
            y[spml + 2] = xp[self.startcol as usize];
            if y[spml + 2] == 1.0 { y[spml + 3] = xp[self.t0col].ln(); }
        } else {
            y[spml + 2] = 0.0;
        }
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

            let mut mu = 0.0f64;
            for i in spmlo1..spmlo { mu += beta[i] * x[i]; }
            let mustar = mu;
            for i in scale..spmlo1 { mu += beta[i] * x[i]; }

            let event: i32 = if x[spmlo + 1] == 1.0 { 1 } else { 0 };
            let start: i32 = if x[spmlo + 2] == 1.0 { 1 } else { 0 };

            // Build log(sigma) in u; then convert to sigma.
            let mut u = 0.0f64;
            for i in 0..scale { u += beta[i] * x[i]; }
            let mut dpmu = 0.0f64;
            let mut ddpmu = 0.0f64;
            let mut tem_power = 1.0f64;
            if ntheta > 0 { dpmu = beta[spmlo]; u += dpmu * mustar; }
            for i in 1..ntheta {
                let tem2 = tem_power * ((i + 1) * i) as f64;
                tem_power *= mustar;
                dpmu  += beta[spmlo + i] * tem_power * (i + 1) as f64;
                ddpmu += beta[spmlo + i] * tem2;
                u     += tem_power * mustar * beta[spmlo + i];
            }
            if !(-39.0..=39.0).contains(&u) { return -1; }
            let sigma = u.exp();
            if event == 1 { *ll1 -= u; }

            // Standardized residual(s).
            let u_val = (x[spmlo] - mu) / sigma;
            let mut f = [0.0f64; 3];
            if self.dist.call(x[spmlo + 1] as i32, u_val, &mut f).is_err() {
                return -1;
            }
            let uf2pf1 = u_val * f[2] + f[1];
            let mut fdif = [f[0], f[1], f[2]];
            let mut u0 = 0.0f64;
            let mut u0f2pf1 = 0.0f64;
            let mut f0 = [0.0f64; 3];
            if start == 1 {
                u0 = (x[spmlo + 3] - mu) / sigma;
                if self.dist.call(0, u0, &mut f0).is_err() { return -1; }
                for i in 0..3 { fdif[i] -= f0[i]; }
                u0f2pf1 = u0 * f0[2] + f0[1];
            }
            *ll1 += fdif[0];

            // Partials for gamma (scale) parameters.
            j = 0;
            let mut tem = f[1] * u_val - f0[1] * u0;
            if event == 1 { tem += 1.0; }
            let tem2_s = u_val * uf2pf1 - u0 * u0f2pf1;
            for i in 0..scale {
                db[i] -= tem * x[i];
                let tem3 = x[i] * tem2_s;
                for ii in 0..=i {
                    ddb[j] -= x[ii] * tem3;
                    j += 1;
                }
            }

            // Partials for etas not coupled into sigma.
            let tem_a  = fdif[1] / sigma;
            let tem2_a = (uf2pf1 - u0f2pf1) / sigma;
            let tem3_a = fdif[2] / (sigma * sigma);
            for i in scale..spmlo1 {
                db[i] -= tem_a * x[i];
                let tem4 = tem2_a * x[i];
                for ii in 0..scale { ddb[j] -= x[ii] * tem4; j += 1; }
                let tem4 = tem3_a * x[i];
                for ii in scale..=i { ddb[j] -= tem4 * x[ii]; j += 1; }
            }

            // Partials for etas coupled into sigma via theta.
            if spmlo > spmlo1 {
                let odsput = 1.0 / sigma + dpmu * u_val;
                let mut tem2 = f[2] / (sigma * sigma)
                    + uf2pf1 * dpmu * (2.0 / sigma + u_val * dpmu)
                    - (f[1] * u_val + event as f64) * ddpmu;
                let mut tem3 = odsput * uf2pf1;
                let mut tem  = odsput * f[1];
                let mut odspu0t = 0.0f64;
                if start == 1 {
                    odspu0t = 1.0 / sigma + dpmu * u0;
                    tem2 -= f0[2] / (sigma * sigma)
                        + u0f2pf1 * dpmu * (2.0 / sigma + u0 * dpmu)
                        - f0[1] * u0 * ddpmu;
                    tem3 -= odspu0t * u0f2pf1;
                    tem  -= odspu0t * f0[1];
                }
                tem += dpmu * event as f64;
                for i in spmlo1..spmlo {
                    db[i] -= tem * x[i];
                    let tem4 = tem3 * x[i];
                    for ii in 0..scale { ddb[j] -= x[ii] * tem4; j += 1; }
                    let tem4 = x[i] * (f[2] * odsput - f0[2] * odspu0t + fdif[1] * dpmu) / sigma;
                    for ii in scale..spmlo1 { ddb[j] -= x[ii] * tem4; j += 1; }
                    let tem4 = x[i] * tem2;
                    for ii in spmlo1..=i { ddb[j] -= x[ii] * tem4; j += 1; }
                }

                // Partials for theta.
                let mut tem_ms = mustar;
                let tem2_t = f[1] * u_val - f0[1] * u0;
                let mut tem5 = -tem2_t - event as f64;
                for i in 0..ntheta {
                    db[spmlo + i] -= tem_ms * (tem2_t + event as f64);
                    let tem3_t = tem_ms * uf2pf1;
                    let tem03  = tem_ms * u0f2pf1;
                    let mut tem4 = u_val * tem3_t - u0 * tem03;
                    for ii in 0..scale { ddb[j] -= x[ii] * tem4; j += 1; }
                    tem4 = (tem3_t - tem03) / sigma;
                    for ii in scale..spmlo1 { ddb[j] -= x[ii] * tem4; j += 1; }
                    tem4 = odsput * tem3_t - odspu0t * tem03 + tem5;
                    for ii in spmlo1..spmlo { ddb[j] -= x[ii] * tem4; j += 1; }
                    let mut tem3_theta = tem_ms * (uf2pf1 * u_val - u0f2pf1 * u0);
                    for _ii in 0..=i {
                        tem3_theta *= mustar;
                        ddb[j] -= tem3_theta;
                        j += 1;
                    }
                    if i + 1 < ntheta {
                        tem_ms *= mustar;
                        tem5 *= mustar * (i + 2) as f64;
                    }
                }
            }
        }
        0
    }
}
