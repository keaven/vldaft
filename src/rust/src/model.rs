// model.rs — Port of dataio8.c + ac_tpar8 from loscreg8.c
// Encapsulates all C global state into ModelState struct

use crate::distributions;

const MAXCOV: usize = 30;

/// Distribution function type
pub type DistFn = fn(i32, f64, &mut [f64; 3]) -> Result<(), i32>;

/// Distribution function type for gamma (needs extra params)
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
            Distribution::Weibull => distributions::weibull(event, x, f),
            Distribution::Logistic => distributions::logistic(event, x, f),
            Distribution::Normal => distributions::normal_dist(event, x, f),
            Distribution::Cauchy => distributions::cauchy_dist(event, x, f),
            Distribution::Gamma { nu, enu } => distributions::gamma_dist(event, x, f, *nu, *enu),
        }
    }
}

/// Model state — replaces all C globals
pub struct ModelState {
    pub m: usize,          // number of columns in data
    pub n: usize,          // number of observations
    pub scale: usize,      // number of scale (gamma) parameters
    pub mlo: usize,        // number of location (eta) parameters
    pub mlo1: usize,       // number of location params NOT in theta coupling
    pub ntheta: usize,     // number of theta parameters
    pub npar: usize,       // total parameters = scale + mlo + ntheta
    pub spmlo: usize,      // scale + mlo
    pub spmlo1: usize,     // scale + mlo1

    pub tcol: usize,       // time column index
    pub cencol: usize,     // event/censoring column index
    pub startcol: i32,     // start indicator column (-1 if none)
    pub t0col: usize,      // start time column

    pub pscale: Vec<usize>, // scale covariate column indices
    pub ploc: Vec<usize>,   // location covariate column indices

    pub reoff: Vec<i32>,   // event indicator values
    pub rgtcen: i32,       // right censoring value

    pub x: Vec<f64>,       // data array (row-major, f64 not float)
    pub xbar: Vec<f64>,    // column means
    pub madj: bool,        // whether to center covariates

    pub dist: Distribution,

    // Iterator state
    xindex: usize,
    y: Vec<f64>,           // current observation record
}

impl ModelState {
    /// Load data from R matrix (column-major doubles) and set up model
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

        // Copy R column-major to row-major f64 and compute means
        let mut x = vec![0.0f64; nobs * ncol];
        let mut xbar = vec![0.0f64; ncol];
        for i in 0..nobs {
            for j in 0..ncol {
                let val = rdata[j * nobs + i];
                x[i * ncol + j] = val;
                xbar[j] += val;
            }
        }
        for j in 0..ncol {
            xbar[j] /= nobs as f64;
        }

        let spmlo = nsc + nlo;
        let mlo1 = if ntheta == 0 {
            nlo
        } else if mlo1_override >= 0 {
            mlo1_override as usize
        } else {
            1  // default: 1 location var not coupled
        };
        let spmlo1 = nsc + mlo1;

        let ploc: Vec<usize> = loc_cols.iter().map(|&c| c as usize).collect();
        let pscale: Vec<usize> = scale_cols.iter().map(|&c| c as usize).collect();

        ModelState {
            m: ncol,
            n: nobs,
            scale: nsc,
            mlo: nlo,
            mlo1,
            ntheta,
            npar: nsc + nlo + ntheta,
            spmlo,
            spmlo1,
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
            xindex: 0,
            y: vec![0.0f64; MAXCOV + 4],
        }
    }

    /// Reset data iterator
    pub fn zrewind(&mut self) {
        self.xindex = 0;
    }

    /// Get next observation record (like C getrec)
    /// Returns None when all records processed
    pub fn getrec(&mut self) -> Option<&[f64]> {
        if self.xindex >= self.n * self.m {
            return None;
        }
        let xp = &self.x[self.xindex..self.xindex + self.m];

        if self.madj {
            if self.scale > 0 {
                self.y[0] = xp[self.pscale[0]];
            }
            for i in 1..self.scale {
                self.y[i] = xp[self.pscale[i]] - self.xbar[self.pscale[i]];
            }
            if self.mlo > 0 {
                self.y[self.scale] = xp[self.ploc[0]];
            }
            for i in 1..self.mlo {
                self.y[self.scale + i] = xp[self.ploc[i]] - self.xbar[self.ploc[i]];
            }
        } else {
            for i in 0..self.scale {
                self.y[i] = xp[self.pscale[i]];
            }
            for i in 0..self.mlo {
                self.y[i + self.scale] = xp[self.ploc[i]];
            }
        }

        let spml = self.scale + self.mlo;
        self.y[spml] = xp[self.tcol].ln();
        self.y[spml + 1] = 0.0;
        for ev in &self.reoff {
            if xp[self.cencol] as i32 == *ev {
                self.y[spml + 1] = 1.0;
            }
        }
        if xp[self.cencol] as i32 == self.rgtcen {
            self.y[spml + 1] = -1.0;
        }
        if self.startcol >= 0 && xp[self.t0col] > 0.0 {
            self.y[spml + 2] = xp[self.startcol as usize];
            if self.y[spml + 2] == 1.0 {
                self.y[spml + 3] = xp[self.t0col].ln();
            }
        } else {
            self.y[spml + 2] = 0.0;
        }

        self.xindex += self.m;
        Some(&self.y[..spml + 4])
    }

    /// Compute partials for AFT model (port of ac_tpar8)
    /// This is the main likelihood/gradient/hessian engine.
    pub fn ac_tpar8(&mut self, m: usize, beta: &[f64], db: &mut [f64], ddb: &mut [f64], ll1: &mut f64) -> i32 {
        let scale = self.scale;
        let mlo = self.mlo;
        let ntheta = self.ntheta;
        let spmlo = scale + mlo;
        let spmlo1 = if spmlo == m { spmlo } else { scale + self.mlo1 };

        // Initialize
        *ll1 = 0.0;
        let mut j = 0usize;
        for i in 0..m {
            db[i] = 0.0;
            for _ii in 0..=i {
                ddb[j] = 0.0;
                j += 1;
            }
        }
        self.zrewind();

        // Process each observation
        loop {
            // We need to borrow self mutably for getrec, then use the data
            // Copy observation into local buffer to avoid borrow issues
            let y = {
                if self.xindex >= self.n * self.m {
                    break;
                }
                // Inline getrec logic to get a copy
                let xp_start = self.xindex;
                let xp = &self.x[xp_start..xp_start + self.m];
                let mut y_local = vec![0.0f64; spmlo + 4];

                if self.madj {
                    if scale > 0 {
                        y_local[0] = xp[self.pscale[0]];
                    }
                    for i in 1..scale {
                        y_local[i] = xp[self.pscale[i]] - self.xbar[self.pscale[i]];
                    }
                    if mlo > 0 {
                        y_local[scale] = xp[self.ploc[0]];
                    }
                    for i in 1..mlo {
                        y_local[scale + i] = xp[self.ploc[i]] - self.xbar[self.ploc[i]];
                    }
                } else {
                    for i in 0..scale {
                        y_local[i] = xp[self.pscale[i]];
                    }
                    for i in 0..mlo {
                        y_local[i + scale] = xp[self.ploc[i]];
                    }
                }

                let spml = spmlo;
                y_local[spml] = xp[self.tcol].ln();
                y_local[spml + 1] = 0.0;
                for ev in &self.reoff {
                    if xp[self.cencol] as i32 == *ev {
                        y_local[spml + 1] = 1.0;
                    }
                }
                if xp[self.cencol] as i32 == self.rgtcen {
                    y_local[spml + 1] = -1.0;
                }
                if self.startcol >= 0 && xp[self.t0col] > 0.0 {
                    y_local[spml + 2] = xp[self.startcol as usize];
                    if y_local[spml + 2] == 1.0 {
                        y_local[spml + 3] = xp[self.t0col].ln();
                    }
                } else {
                    y_local[spml + 2] = 0.0;
                }

                self.xindex += self.m;
                y_local
            };

            let x = &y;

            // Compute mu
            let mut mu = 0.0f64;
            for i in spmlo1..spmlo {
                mu += beta[i] * x[i];
            }
            let mustar = mu;
            for i in scale..spmlo1 {
                mu += beta[i] * x[i];
            }

            // Event indicator
            let event: i32 = if x[spmlo + 1] == 1.0 { 1 } else { 0 };
            let start: i32 = if x[spmlo + 2] == 1.0 { 1 } else { 0 };

            // Compute sigma
            let mut u = 0.0f64;
            for i in 0..scale {
                u += beta[i] * x[i];
            }
            let mut dpmu = 0.0f64;
            let mut ddpmu = 0.0f64;
            let mut tem_power = 1.0f64;
            if ntheta > 0 {
                dpmu = beta[spmlo];
                u += dpmu * mustar;
            }
            for i in 1..ntheta {
                let tem2 = tem_power * ((i + 1) * i) as f64;
                tem_power *= mustar;
                dpmu += beta[spmlo + i] * tem_power * (i + 1) as f64;
                ddpmu += beta[spmlo + i] * tem2;
                u += tem_power * mustar * beta[spmlo + i];
            }
            if u > 39.0 || u < -39.0 {
                return -1;
            }
            let sigma = u.exp();
            if event == 1 {
                *ll1 -= u;
            }

            // Compute standardized residual
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
                if self.dist.call(0, u0, &mut f0).is_err() {
                    return -1;
                }
                for i in 0..3 {
                    fdif[i] -= f0[i];
                }
                u0f2pf1 = u0 * f0[2] + f0[1];
            }
            *ll1 += fdif[0];

            // 1st and 2nd partials for gamma's (scale params)
            j = 0;
            let mut tem = f[1] * u_val - f0[1] * u0;
            if event == 1 { tem += 1.0; }
            let tem2_val = u_val * uf2pf1 - u0 * u0f2pf1;
            for i in 0..scale {
                db[i] -= tem * x[i];
                let tem3 = x[i] * tem2_val;
                for ii in 0..=i {
                    ddb[j] -= x[ii] * tem3;
                    j += 1;
                }
            }

            // Partials for etas not in sigma
            let tem_a = fdif[1] / sigma;
            let tem2_a = (uf2pf1 - u0f2pf1) / sigma;
            let tem3_a = fdif[2] / sigma / sigma;
            for i in scale..spmlo1 {
                db[i] -= tem_a * x[i];
                let tem4 = tem2_a * x[i];
                for ii in 0..scale {
                    ddb[j] -= x[ii] * tem4;
                    j += 1;
                }
                let tem4 = tem3_a * x[i];
                for ii in scale..=i {
                    ddb[j] -= tem4 * x[ii];
                    j += 1;
                }
            }

            // Partials for etas coupled to sigma via theta
            if spmlo > spmlo1 {
                let odsput = 1.0 / sigma + dpmu * u_val;
                let mut tem2 = f[2] / sigma / sigma + uf2pf1 * dpmu * (2.0 / sigma + u_val * dpmu)
                    - (f[1] * u_val + event as f64) * ddpmu;
                let mut tem3 = odsput * uf2pf1;
                let mut tem = odsput * f[1];
                let mut odspu0t = 0.0f64;
                if start == 1 {
                    odspu0t = 1.0 / sigma + dpmu * u0;
                    tem2 -= f0[2] / sigma / sigma + u0f2pf1 * dpmu * (2.0 / sigma + u0 * dpmu)
                        - f0[1] * u0 * ddpmu;
                    tem3 -= odspu0t * u0f2pf1;
                    tem -= odspu0t * f0[1];
                }
                tem += dpmu * event as f64;
                for i in spmlo1..spmlo {
                    db[i] -= tem * x[i];
                    let tem4 = tem3 * x[i];
                    for ii in 0..scale {
                        ddb[j] -= x[ii] * tem4;
                        j += 1;
                    }
                    let tem4 = x[i] * (f[2] * odsput - f0[2] * odspu0t + fdif[1] * dpmu) / sigma;
                    for ii in scale..spmlo1 {
                        ddb[j] -= x[ii] * tem4;
                        j += 1;
                    }
                    let tem4 = x[i] * tem2;
                    for ii in spmlo1..=i {
                        ddb[j] -= x[ii] * tem4;
                        j += 1;
                    }
                }

                // Partials for thetas
                let mut tem_ms = mustar;
                let mut tem2_t = f[1] * u_val - f0[1] * u0;
                let mut tem5 = -tem2_t - event as f64;
                for i in 0..ntheta {
                    db[spmlo + i] -= tem_ms * (tem2_t + event as f64);
                    let tem3_t = tem_ms * uf2pf1;
                    let tem03 = tem_ms * u0f2pf1;
                    let mut tem4 = u_val * tem3_t - u0 * tem03;
                    for ii in 0..scale {
                        ddb[j] -= x[ii] * tem4;
                        j += 1;
                    }
                    tem4 = (tem3_t - tem03) / sigma;
                    for ii in scale..spmlo1 {
                        ddb[j] -= x[ii] * tem4;
                        j += 1;
                    }
                    tem4 = odsput * tem3_t - odspu0t * tem03 + tem5;
                    for ii in spmlo1..spmlo {
                        ddb[j] -= x[ii] * tem4;
                        j += 1;
                    }
                    let mut tem3_theta = tem_ms * (uf2pf1 * u_val - u0f2pf1 * u0);
                    for ii in 0..=i {
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
