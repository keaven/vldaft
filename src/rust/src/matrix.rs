//! Matrix algebra for the Newton-Raphson solver.
//!
//! All symmetric matrices are stored packed, lower-triangle by rows:
//! index(i, j) = i * (i + 1) / 2 + j for j <= i.

/// Invert a lower-triangular matrix stored packed.
pub fn triinv(x: &[f64], y: &mut [f64], work: &mut [f64], n: usize) {
    let mut buf = vec![0.0f64; n];
    for i in (0..n).rev() {
        for j in 0..i { work[j] = 0.0; }
        work[i] = 1.0;
        // ubak with work aliased as both z and y: take a copy first.
        let row = i + 1;
        buf[..row].copy_from_slice(&work[..row]);
        ubak(x, work, &buf[..row], 1, row);
        let k = i * (i + 1) / 2;
        for j in 0..=i { y[k + j] = work[j]; }
    }
}

/// Multiply a packed lower-triangular matrix by its transpose.
pub fn trimult(y: &[f64], z: &mut [f64], n: usize) {
    for i in 0..n {
        for j in i..n {
            let k = j * (j + 1) / 2 + i;
            let mut m2 = j * (j + 1) / 2 + j;
            z[k] = y[k] * y[m2];
            let mut k2 = k;
            for l in (j + 1)..n {
                k2 += l;
                m2 += l;
                z[k] += y[k2] * y[m2];
            }
        }
    }
}

/// Plain Cholesky decomposition. Returns 1 if not positive definite.
pub fn chol(s: &[f64], t: &mut [f64], n: usize) -> i32 {
    let mut ii: i64 = -1;
    for i in 0..n {
        let inc = i + 1;
        let i0 = (ii + 1) as usize;
        ii += inc as i64;
        let ii_u = ii as usize;
        let mut t1 = s[ii_u];
        for ik in i0..ii_u { t1 -= t[ik] * t[ik]; }
        if t1 <= 0.0 { return 1; }
        t1 = t1.sqrt();
        t[ii_u] = t1;
        let mut ji = ii_u;
        let mut inc2 = inc;
        for _j in (i + 1)..n {
            ji += inc2;
            inc2 += 1;
            let jk_start = ji - i;
            let mut ik = i0;
            let mut t2 = s[ji];
            for k in 0..i {
                t2 -= t[ik] * t[jk_start + k];
                ik += 1;
            }
            t[ji] = t2 / t1;
        }
    }
    0
}

/// Modified (Gill-Murray-Wright) Cholesky decomposition.
/// Returns 1 if the diagonal had to be inflated, 0 otherwise.
pub fn cholmod(s: &[f64], t: &mut [f64], d: &mut [f64], n: usize) -> i32 {
    let macheps = f64::EPSILON;
    let mut retval = 0i32;

    let mut j = 0usize;
    let mut psi = 0.0f64;
    let mut bup = 0.0f64;
    let mut t2_max = 0.0f64;
    for i in 0..n {
        for _ii in 0..i {
            psi += 2.0 * s[j] * s[j];
            let t1 = s[j].abs();
            j += 1;
            bup = bup.max(t1);
        }
        let t1 = s[j].abs();
        t2_max = t2_max.max(t1);
        psi += s[j] * s[j];
        j += 1;
    }
    psi = psi.sqrt();
    psi = if psi > 1.0 { psi * macheps } else { macheps };
    bup /= n as f64;
    bup = bup.max(t2_max).max(macheps);

    let mut ii: i64 = -1;
    for i in 0..n {
        let i0 = (ii + 1) as usize;
        ii += (i + 1) as i64;
        let ii_u = ii as usize;
        let mut t1 = s[ii_u];
        for ik in i0..ii_u { t1 -= t[ik] * t[ik] * d[ik - i0]; }
        if t1 <= 0.0 { t1 = -t1; retval = 1; }
        if t1 < psi { t1 = psi; retval = 1; }
        let mut ji = ii_u;
        let mut tj = 0.0f64;
        for j2 in (i + 1)..n {
            ji += j2;
            let jk_start = ji - i;
            let mut ik = i0;
            let mut t2 = s[ji];
            for k in 0..i {
                t2 -= t[ik] * t[jk_start + k] * d[k];
                ik += 1;
            }
            t[ji] = t2;
            tj = tj.max(t2.abs());
        }
        let tem = tj * tj / bup;
        if t1 < tem { t1 = tem; retval = 1; }
        d[i] = t1;
        t[ii_u] = 1.0;
        let mut ji2 = ii_u;
        for j2 in (i + 1)..n {
            ji2 += j2;
            t[ji2] /= t1;
        }
    }
    retval
}

/// Solve an upper-triangular system (packed storage).
pub fn ubak(t: &[f64], z: &mut [f64], y: &[f64], m: usize, n: usize) -> i32 {
    let mut iy = n * m - 1;
    let mut it2 = n * (n + 1) / 2;
    for i1 in 0..n {
        it2 -= 1;
        for i2 in 0..m {
            let mut tem = y[iy];
            let mut it = it2;
            let mut iz = n * m - 1 - i2;
            for i3 in 0..i1 {
                tem -= z[iz] * t[it];
                iz -= m;
                it -= n - i3 - 1;
            }
            if tem != 0.0 && t[it] == 0.0 { return 1; }
            z[iy] = tem / t[it];
            if iy == 0 { return 0; }
            iy -= 1;
        }
    }
    0
}

/// Solve a lower-triangular system (packed storage).
pub fn lbak(t: &[f64], z: &mut [f64], y: &[f64], m: usize, n: usize) -> i32 {
    let mut iy = 0usize;
    for i1 in 0..n {
        for i2 in 0..m {
            let mut tem = y[iy];
            let mut it = (i1 + 1) * i1 / 2;
            let mut iz = i2;
            for _i3 in 0..i1 {
                tem -= z[iz] * t[it];
                iz += m;
                it += 1;
            }
            if tem != 0.0 && t[it] == 0.0 { return 1; }
            z[iy] = tem / t[it];
            iy += 1;
        }
    }
    0
}

/// Multiply a packed symmetric matrix (m x m) by a rectangular m x n matrix.
pub fn symmult(x: &[f64], y: &[f64], z: &mut [f64], m: usize, n: usize) {
    let mut l = 0usize;
    for i in 0..m {
        let item = i * (i + 1) / 2;
        for j in 0..n {
            z[l] = 0.0;
            let mut ii2 = item;
            let mut jj = j;
            for _k in 0..i {
                z[l] += x[ii2] * y[jj];
                ii2 += 1;
                jj += n;
            }
            let ij = m - i;
            let mut it = i + 1;
            for _k in 0..ij {
                z[l] += x[ii2] * y[jj];
                jj += n;
                ii2 += it;
                it += 1;
            }
            l += 1;
        }
    }
}

/// General rectangular matrix product: z (k x n) = x (k x m) * y (m x n).
pub fn mmult1(x: &[f64], y: &[f64], z: &mut [f64], k: usize, m: usize, n: usize) {
    let mut i = 0usize;
    for kk in 0..k {
        for nn in 0..n {
            z[i] = 0.0;
            let mut l = nn;
            let mut j2 = kk * m;
            for _mm in 0..m {
                z[i] += x[j2] * y[l];
                j2 += 1;
                l += n;
            }
            i += 1;
        }
    }
}
