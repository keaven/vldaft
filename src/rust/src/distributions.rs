// distributions.rs — Port of loscreg8.c distribution functions
// Each function takes (event, x) and writes [f0, f1, f2] to the output array.
// Returns Ok(()) on success, Err(-1) on numerical failure.

use libm::{exp, log, sqrt, atan, lgamma};

/// Weibull distribution: extreme value density/survival/left-truncation
pub fn weibull(event: i32, x: f64, f: &mut [f64; 3]) -> Result<(), i32> {
    if x > 39.0 || x < -39.0 {
        return Err(-1);
    }
    let eu = exp(x);
    if event >= 0 {
        f[0] = -eu;
        f[1] = -eu;
        f[2] = -eu;
        if event == 1 {
            f[0] += x;
            f[1] += 1.0;
        }
    } else {
        let emeu = exp(-eu);
        let omemeu = 1.0 - emeu;
        f[0] = log(omemeu);
        f[1] = eu * emeu / omemeu;
        f[2] = f[1] * (1.0 - eu - f[1]);
    }
    Ok(())
}

/// Gamma distribution
pub fn gamma_dist(event: i32, x: f64, f: &mut [f64; 3], nu: f64, enu: f64) -> Result<(), i32> {
    if x > 39.0 {
        return Err(-1);
    }
    let eu = exp(x);
    if event == 1 {
        f[0] = x * enu - eu - lgamma(enu);
        f[1] = enu - eu;
        f[2] = -eu;
    } else if event == 0 {
        let lcg = lgamma(enu);
        let cg = exp(lcg);
        let icg = if (enu <= eu && eu <= 1.0) || eu < enu {
            let val = cg - pgamma_reg(eu, enu, true) * cg;
            if val <= 0.0 { return Err(-1); }
            val
        } else {
            let val = pgamma_reg(eu, enu, false) * cg;
            if val <= 0.0 { return Err(-1); }
            val
        };
        f[0] = log(icg) - lcg;
        f[1] = -exp(x * enu - eu) / icg;
        f[2] = f[1] * (enu - eu - f[1]);
    } else {
        if x < -39.0 { return Err(-1); }
        let lcg = lgamma(enu);
        let cg = exp(lcg);
        let icg = if (enu <= eu && eu <= 1.0) || eu < enu {
            let val = pgamma_reg(eu, enu, true) * cg;
            if val <= 0.0 { return Err(-1); }
            val
        } else {
            let val = cg - pgamma_reg(eu, enu, false) * cg;
            if val <= 0.0 { return Err(-1); }
            val
        };
        f[0] = log(icg) - lcg;
        f[1] = exp(x * enu - eu) / icg;
        f[2] = f[1] * (enu - eu - f[1]);
    }
    Ok(())
}

/// Logistic distribution
pub fn logistic(event: i32, x: f64, f: &mut [f64; 3]) -> Result<(), i32> {
    if x > 39.0 || x < -39.0 {
        return Err(-1);
    }
    let x_adj = if event < 0 { -x } else { x };
    let eu = exp(x_adj);
    let opeu = 1.0 + eu;
    let opc = if event == 1 { 2.0 } else { 1.0 };
    f[0] = if event == 1 { x_adj } else { 0.0 } - opc * log(opeu);
    f[1] = if event == 1 { 1.0 } else { 0.0 } - opc * eu / opeu;
    f[2] = -(opc * eu / opeu) / opeu;
    // C code: f[1] = (event==1) - f[1_temp]; but we computed f[1] directly above
    // Actually let me re-port exactly:
    let tmp1 = opc * eu / opeu;
    f[0] = if event == 1 { x_adj } else { 0.0 } - opc * log(opeu);
    f[2] = -tmp1 / opeu;
    f[1] = if event == 1 { 1.0 } else { 0.0 } - tmp1;
    if event < 0 {
        f[1] = -f[1];
    }
    Ok(())
}

/// Normal distribution
pub fn normal_dist(event: i32, u: f64, f: &mut [f64; 3]) -> Result<(), i32> {
    if event == 1 {
        f[0] = log_dnorm(u);
        f[1] = -u;
        f[2] = -1.0;
        return Ok(());
    }
    let u_adj = if event < 0 { -u } else { u };

    if u_adj < -8.0 {
        f[0] = 0.0;
        f[1] = 0.0;
        f[2] = 0.0;
        return Ok(());
    } else if u_adj < 0.0 {
        let phi_val = pnorm_cdf(u_adj);
        let phi = dnorm_pdf(u_adj);
        f[0] = log(phi_val);
        f[1] = -phi / phi_val;
        f[2] = -f[1] * (u_adj + f[1]);
        if event < 0 { f[1] = -f[1]; }
    } else {
        let log_surv = log_pnorm_upper(u_adj);
        let phi = dnorm_pdf(u_adj);
        f[0] = log_surv;
        f[1] = -phi / exp(log_surv);
        f[2] = -f[1] * (u_adj + f[1]);
        if event < 0 { f[1] = -f[1]; }
    }
    Ok(())
}

/// Cauchy distribution
pub fn cauchy_dist(event: i32, x: f64, f: &mut [f64; 3]) -> Result<(), i32> {
    let x2 = x * x;
    let odopx2 = 1.0 / (1.0 + x2);
    if event == 1 {
        f[0] = -1.144729886 + log(odopx2);
        f[1] = -2.0 * x * odopx2;
        f[2] = (4.0 * x2 * odopx2 - 2.0) * odopx2;
        return Ok(());
    }
    let x_adj = if event < 0 { -x } else { x };
    let t = 0.5 + atan(-x_adj) * 0.318309886;
    if t <= 0.0 {
        return Err(-1);
    }
    let odopx2_adj = 1.0 / (1.0 + x_adj * x_adj);
    f[0] = log(t);
    f[1] = -odopx2_adj * 0.318309886 / t;
    f[2] = f[1] * (-x_adj * 2.0 * odopx2_adj - f[1]);
    if event < 0 {
        f[1] = -f[1];
    }
    Ok(())
}

// ---- Helper functions for normal distribution ----

/// Standard normal PDF
fn dnorm_pdf(x: f64) -> f64 {
    const INV_SQRT_2PI: f64 = 0.3989422804014327;
    INV_SQRT_2PI * exp(-0.5 * x * x)
}

/// log of standard normal PDF
fn log_dnorm(x: f64) -> f64 {
    const LOG_INV_SQRT_2PI: f64 = -0.9189385332046727;
    LOG_INV_SQRT_2PI - 0.5 * x * x
}

/// Standard normal CDF using rational approximation (Abramowitz & Stegun)
pub fn pnorm_cdf(x: f64) -> f64 {
    // Use the error function approach
    0.5 * erfc(-x / std::f64::consts::SQRT_2)
}

/// log of upper tail P(Z > x)
fn log_pnorm_upper(x: f64) -> f64 {
    if x < 6.0 {
        log(1.0 - pnorm_cdf(x))
    } else {
        // Mills ratio approximation for large x
        log_dnorm(x) - log(x) + log(1.0 - 1.0 / (x * x))
    }
}

/// Complementary error function
fn erfc(x: f64) -> f64 {
    // Use a good rational approximation
    if x >= 0.0 {
        erfc_positive(x)
    } else {
        2.0 - erfc_positive(-x)
    }
}

/// erfc for x >= 0 using Horner form of rational approximation
fn erfc_positive(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t * (0.254829592
        + t * (-0.284496736
            + t * (1.421413741
                + t * (-1.453152027
                    + t * 1.061405429))));
    poly * exp(-x * x)
}

/// Regularized incomplete gamma P(a, x) or Q(a, x)
/// lower_tail=true: P(a,x) = integral from 0 to x
/// lower_tail=false: Q(a,x) = 1 - P(a,x)
fn pgamma_reg(x: f64, a: f64, lower_tail: bool) -> f64 {
    if x <= 0.0 {
        return if lower_tail { 0.0 } else { 1.0 };
    }
    // Use series or continued fraction
    if x < a + 1.0 {
        let p = gamma_series(a, x);
        if lower_tail { p } else { 1.0 - p }
    } else {
        let q = gamma_cf(a, x);
        if lower_tail { 1.0 - q } else { q }
    }
}

/// Regularized lower incomplete gamma by series expansion
fn gamma_series(a: f64, x: f64) -> f64 {
    let mut sum = 1.0 / a;
    let mut term = 1.0 / a;
    for n in 1..200 {
        term *= x / (a + n as f64);
        sum += term;
        if term.abs() < sum.abs() * 1e-15 { break; }
    }
    sum * exp(-x + a * log(x) - lgamma(a))
}

/// Regularized upper incomplete gamma by continued fraction (Lentz)
fn gamma_cf(a: f64, x: f64) -> f64 {
    let mut f = 1e-30_f64;
    let mut c = 1e-30_f64;
    let mut d: f64;
    let b0 = x + 1.0 - a;
    d = 1.0 / b0;
    f = d;
    for n in 1..200 {
        let nf = n as f64;
        let an = nf * (a - nf);
        let bn = x + 2.0 * nf + 1.0 - a;
        d = bn + an * d;
        if d.abs() < 1e-30 { d = 1e-30; }
        d = 1.0 / d;
        c = bn + an / c;
        if c.abs() < 1e-30 { c = 1e-30; }
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-15 { break; }
    }
    f * exp(-x + a * log(x) - lgamma(a))
}
