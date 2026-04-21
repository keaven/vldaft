/* Define compiler, matrix dimensions, machine epsilon */
#include "rj8def.h"

/* triinv: compute inverse of a triangular matrix */
void triinv(double *x, double *y, double *work, int n)
{
    int i, j, k;
    for (i = n - 1; i >= 0; i--)
    {   for (j = 0; j < i; j++) work[j] = 0;
        work[i] = 1.;
        ubak(x, work, work, 1, i + 1);
        k = i * (i + 1) / 2;
        for (j = 0; j <= i; j++) y[k + j] = work[j];
    }
}


/* trimult: triangular matrix times its transpose */
void trimult(double *y, double *z, int n)
{
    int i, j, k, l, k2, m2;
    for (i = 0; i < n; i++)
    {   for (j = i; j < n; j++)
        {   k = j * (j + 1) / 2 + i;
            m2 = j * (j + 1) / 2 + j;
            z[k] = y[k] * y[m2];
            k2 = k;
            for (l = j + 1; l < n; l++)
            {   k2 += l;
                m2 += l;
                z[k] += y[k2] * y[m2];
            }
        }
    }
}


/* chol: compute cholesky decomposition */
int chol(double *s, double *t, int n)
{
    double t1, t2;
    int i, j, k, ii, ji, ik, jk, inc, i0;
    ii = -1;
    for (i = 0; i < n; i++)
    {   inc = i + 1;
        i0 = ii + 1;
        ii += inc;
        t1 = s[ii];
        for (ik = i0; ik < ii; ik++)
            t1 -= t[ik] * t[ik];
        if (t1 <= 0.) return(1);
        t1 = sqrt(t1);
        t[ii] = t1;
        ji = ii;
        for (j = i + 1; j < n; j++)
        {   ji += inc;
            inc++;
            jk = ji - i;
            ik = i0;
            t2 = s[ji];
            for (k = 0; k < i; k++)
            {   t2 -= t[ik] * t[jk];
                ik++;
                jk++;
            }
            t[ji] = t2 / t1;
        }
    }
    return(0);
}


/* cholmod: modified cholesky decomposition */
int cholmod(double *s, double *t, double *d, int n)
{
    double t1, t2, tj, tem, psi, bup;
    int i, j, k, ii, ji, ik, jk, i0, retval;

    j = 0; psi = 0.; bup = 0.; t2 = 0.; retval = 0;
    for (i = 0; i < n; i++)
    {   for (ii = 0; ii < i; ii++)
        {   psi += 2. * s[j] * s[j];
            t1 = (s[j] > 0. ? s[j] : -s[j]);
            j++;
            bup = (t1 > bup ? t1 : bup);
        }
        t1 = (s[j] > 0 ? s[j] : -s[j]);
        t2 = (t1 > t2 ? t1 : t2);
        psi += s[j] * s[j++];
    }
    psi = sqrt(psi);
    psi = (psi > 1. ? psi * MACHEPS : MACHEPS);
    bup /= n;
    bup = (bup > t2 ? bup : t2);
    bup = (bup > MACHEPS ? bup : MACHEPS);

    ii = -1;
    for (i = 0; i < n; i++)
    {   i0 = ii + 1; ii += i + 1;
        t1 = s[ii]; for (ik = i0; ik < ii; ik++) t1 -= t[ik] * t[ik] * d[ik - i0];
        if (t1 <= 0.) { t1 = -t1; retval = 1; }
        if (t1 < psi) { t1 = psi; retval = 1; }
        ji = ii; tj = 0.;
        for (j = i + 1; j < n; j++)
        {   ji += j; jk = ji - i; ik = i0;
            t2 = s[ji]; for (k = 0; k < i; k++) { t2 -= t[ik] * t[jk] * d[k]; ik++; jk++; }
            t[ji] = t2;
            tem = (t2 > 0. ? t2 : -t2);
            tj = (tem > tj ? tem : tj);
        }
        tem = tj * tj / bup;
        if (t1 < tem) { t1 = tem; retval = 1; }
        d[i] = t1; t[ii] = 1.;
        ji = ii; for (j = i + 1; j < n; j++) { ji += j; t[ji] /= t1; }
    }
    return(retval);
}


/* ubak: find solution to upper triangular set of equations */
int ubak(double *t, double *z, double *y, int m, int n)
{
    double tem;
    int i1, i2, i3, it, it2, iy, iz;
    iy = n * m - 1;
    it2 = n * (n + 1) / 2;
    for (i1 = 0; i1 < n; i1++)
    {   it2 -= 1;
        for (i2 = 0; i2 < m; i2++)
        {   tem = y[iy];
            it = it2;
            iz = n * m - 1 - i2;
            for (i3 = 0; i3 < i1; i3++)
            {   tem -= z[iz] * t[it];
                iz -= m;
                it -= (n - i3 - 1);
            }
            if (tem != 0 && t[it] == 0) return(1);
            z[iy] = tem / t[it];
            iy -= 1;
        }
    }
    return(0);
}


/* lbak: find solution to lower triangular set of equations */
int lbak(double *t, double *z, double *y, int m, int n)
{
    double tem;
    int i1, i2, i3, it, iz, iy;
    iy = 0;
    for (i1 = 0; i1 < n; i1++)
    {   for (i2 = 0; i2 < m; i2++)
        {   tem = y[iy];
            it = (i1 + 1) * i1 / 2;
            iz = i2;
            for (i3 = 0; i3 < i1; i3++)
            {   tem -= z[iz] * t[it];
                iz += m;
                it++;
            }
            if (tem != 0 && t[it] == 0) return(1);
            z[iy] = tem / t[it];
            iy++;
        }
    }
    return(0);
}


/* spdinv: inverse of symmetric, positive definite matrix */
int spdinv(double *t, double *y, int n)
{
    double z1[MAXCOV];
    if (chol(t, y, n) == 1)
        return(1);
    triinv(y, y, z1, n);
    trimult(y, y, n);
    return(0);
}


/* symmult: multiplies symmetric matrix times rectangular */
void symmult(double *x, double *y, double *z, int m, int n)
{
    int i, j, k, l, ij, item, ii2, jj, it;
    l = 0;
    for (i = 0; i < m; i++)
    {   item = i * (i + 1) / 2;
        for (j = 0; j < n; j++)
        {   z[l] = 0;
            ii2 = item;
            jj = j;
            for (k = 0; k < i; k++)
            {   z[l] += x[ii2++] * y[jj];
                jj += n;
            }
            ij = m - i;
            it = i + 1;
            for (k = 0; k < ij; k++)
            {   z[l] += x[ii2] * y[jj];
                jj += n;
                ii2 += it++;
            }
            l++;
        }
    }
}


/* mmult1: multiplies two rectangular matrices */
void mmult1(double *x, double *y, double *z, int k, int m, int n)
{
    int kk, mm, nn, i, j2, l;
    i = -1;
    for (kk = 0; kk < k; kk++)
    {   for (nn = 0; nn < n; nn++)
        {   z[++i] = 0;
            l = nn;
            j2 = kk * m;
            for (mm = 0; mm < m; mm++)
            {   z[i] += x[j2++] * y[l];
                l += n;
            }
        }
    }
}


/* sextract: extracts a portion of a symmetric matrix into a symmetric matrix */
void sextract(int m, double *symin, int n, double *symout, int *rows)
{
    int i, ii, j, j2, j3, k;
    k = 0;
    for (i = 0; i < n; i++)
    {   j = rows[i]; j = j * (j + 1) / 2;
        j2 = j + rows[i];
        for (ii = 0; ii < i; ii++)
        {   if (rows[ii] < rows[i]) symout[k++] = symin[j + rows[ii]];
            else { j3 = rows[ii]; j3 = j3 * (j3 + 1) / 2; symout[k++] = symin[j3 + rows[i]]; }
        }
        symout[k++] = symin[j2];
    }
}


/* rextract: extracts a portion of a symmetric matrix into a rectangular matrix */
void rextract(int m, double *symin, int nr, int nc, double *recout, int *rows, int *cols)
{
    int i, ii, j, j2, k;
    k = 0;
    for (i = 0; i < nr; i++)
    {   j = rows[i]; j = j * (j + 1) / 2;
        for (ii = 0; ii < nc; ii++)
        {   if (cols[ii] <= rows[i]) recout[k++] = symin[j + cols[ii]];
            else { j2 = cols[ii]; j2 = j2 * (j2 + 1) / 2; recout[k++] = symin[j2 + rows[i]]; }
        }
    }
}


/* rectri: takes the lower tri. part of sq. matrix into triangular */
void rectri(int m, double *rec, double *tri)
{
    int i, ii, j, k;
    k = 0;
    for (i = 0; i < m; i++)
    {   j = i * m;
        for (ii = 0; ii <= i; ii++) tri[k++] = rec[j++];
    }
}


/* transpose */
void transpose(double *x, double *y, int m, int n)
{
    int i, j, ix, iy;
    ix = 0;
    for (i = 0; i < m; i++)
    {   iy = i;
        for (j = 0; j < n; j++)
        {   y[iy] = x[ix++];
            iy += m;
        }
    }
}


/* wald: Wald's hypothesis test */
double wald(int m, double *ddb, double *cov, int n, int *in, int *out, double *beta)
{
    int i;
    double t1[MAXREC], t2[MAXREC], t3[MAXREC], s1[MAXSS], s2[MAXSS];
    if (n > 0 && n < m)
    {   sextract(m, cov, m - n, s1, out);
        spdinv(s1, s2, m - n);
        for (i = 0; i < m - n; i++) t1[i] = beta[out[i]];
        symmult(s2, t1, t2, m - n, 1);
        mmult1(t1, t2, t3, 1, m - n, 1);
        return(*t3);
    }
    else
    {   symmult(ddb, beta, t1, m, 1);
        mmult1(beta, t1, t2, 1, m, 1);
        return(*t2);
    }
}
