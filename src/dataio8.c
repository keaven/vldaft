/* dataio8.c -- data loading and record iteration for the AFT engine.
   Backing storage is double-precision; the R caller passes a column-major
   numeric matrix and we keep a row-major double copy plus column means. */

#include "rj8def.h"

/* Module-wide model state.
   The C backend is intentionally single-threaded; protect any concurrent
   access at the R level. */
int  m = 0, scale = 0, mlo = 0, mlo1 = 0, ntheta = 0, tcol = 0, cencol = 0,
     *pscale = NULL, *ploc = NULL, npar = 0, nreoff = 1, *reoff = NULL,
     madj = 1, rgtcen = -1, spmlo1 = 0, spmlo = 0,
     startcol = -1, t0col = 0;
long n = 0, xindex = 0;
double *x = NULL, *y = NULL, *xbar = NULL;

/* Distribution function pointer and gamma parameter */
int (*dist)(int, double, double *) = NULL;
double nu = 0.0, enu = 0.0;


/* getrec: emit the next observation record, with covariate centering if
   madj is set. Layout of the returned row (length scale + mlo + 4):

     [0 .. scale-1]              : scale covariates (gamma)
     [scale .. scale+mlo-1]      : location covariates (eta)
     [scale+mlo]                 : log(time)
     [scale+mlo + 1]             : event indicator (1 event, 0 right-cens, -1 left-cens)
     [scale+mlo + 2]             : left-truncation flag (1 if truncated)
     [scale+mlo + 3]             : log(start time), only valid when flag == 1
*/
double *getrec(void)
{
    int i, spml;
    double *xp;
    if (xindex >= n * m) return NULL;
    xp = x + xindex;

    if (madj)
    {   if (scale > 0) y[0] = xp[pscale[0]];
        for (i = 1; i < scale; i++) y[i] = xp[pscale[i]] - xbar[pscale[i]];
        if (mlo > 0) y[scale] = xp[ploc[0]];
        for (i = 1; i < mlo; i++) y[scale + i] = xp[ploc[i]] - xbar[ploc[i]];
    }
    else
    {   for (i = 0; i < scale; i++) y[i]         = xp[pscale[i]];
        for (i = 0; i < mlo;   i++) y[i + scale] = xp[ploc[i]];
    }
    spml = scale + mlo;
    y[spml]     = log(xp[tcol]);
    y[spml + 1] = 0.;
    for (i = 0; i < nreoff; i++) if ((int)xp[cencol] == reoff[i]) y[spml + 1] = 1.;
    if ((int)xp[cencol] == rgtcen) y[spml + 1] = -1.;
    if (startcol >= 0 && xp[t0col] > 0.)
    {   y[spml + 2] = (double)(int)xp[startcol];
        if (y[spml + 2] == 1.) y[spml + 3] = log(xp[t0col]);
    }
    else y[spml + 2] = 0.;
    xindex += m;
    return y;
}

void zrewind(void) { xindex = 0; }


/* setupls: fill in module-wide model configuration. */
int setupls(int ncol, long nrec, int nlo, int nsc, int nth, int time,
            int censor, int *p1, int *p2, int nevent, int *event,
            int stcol, int time0)
{
    m       = ncol;
    n       = nrec;
    scale   = nsc;
    mlo     = nlo;
    ntheta  = nth;
    mlo1    = (ntheta == 0) ? mlo : 1;
    spmlo   = scale + mlo;
    spmlo1  = scale + mlo1;
    npar    = scale + mlo + ntheta;
    pscale  = p2;
    ploc    = p1;
    tcol    = time;
    t0col   = time0;
    startcol = stcol;
    cencol  = censor;
    nreoff  = nevent;
    reoff   = event;
    return 0;
}


/* data_from_r: copy R's column-major double matrix into our row-major store
   and compute per-column means. Uses R_alloc so R manages the lifetime. */
int data_from_r(double *rdata, int nobs, int ncol)
{
    int i, j;

    n = (long)nobs;
    m = ncol;

    /* R_alloc memory is reclaimed at the end of each .Call(), so always
       re-allocate; do not rely on the previous call's pointers. */
    y    = (double *)R_alloc(MAXCOV + 4, sizeof(double));
    x    = (double *)R_alloc((size_t)nobs * ncol, sizeof(double));
    xbar = (double *)R_alloc(ncol, sizeof(double));
    for (j = 0; j < ncol; j++) xbar[j] = 0.;

    for (i = 0; i < nobs; i++)
    {   for (j = 0; j < ncol; j++)
        {   double v = rdata[(size_t)j * nobs + i]; /* R matrices are column-major */
            x[(size_t)i * ncol + j] = v;
            xbar[j] += v;
        }
    }
    for (j = 0; j < ncol; j++) xbar[j] /= (double)n;

    return 0;
}
