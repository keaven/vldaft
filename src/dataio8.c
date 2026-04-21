/* DATAIO8.C - Data I/O and model setup for R interface
   Refactored to receive data from R via .Call() instead of flat files.
*/

#include "rj8def.h"

/* Global variables */
int m = 0, scale = 0, mlo = 0, mlo1 = 0, ntheta = 0, tcol = 0, cencol = 0,
    *pscale = NULL, *ploc = NULL, npar = 0, nreoff = 1, *reoff = NULL,
    madj = 1, rgtcen = -1, spmlo1 = 0, spmlo = 0,
    startcol = -1, t0col = 0, boot = 0, *bootsamp = NULL, obs = 0;
long n = 0, xindex = 0;
float *x = NULL, *xp = NULL, *y = NULL;
double *xbar = NULL;

/* Distribution function pointer and gamma parameter */
int (*dist)(int, double, double *) = NULL;
double nu = 0.0, enu = 0.0;


/* getrec: return next observation record */
float *getrec(void)
{
    int i, spml;
    if (xindex >= n * m) return(NULL);
    xp = x + xindex;

    if (madj)
    {   if (scale > 0) y[0] = xp[pscale[0]];
        for (i = 1; i < scale; i++) y[i] = xp[pscale[i]] - xbar[pscale[i]];
        if (mlo > 0) y[scale] = xp[ploc[0]];
        for (i = 1; i < mlo; i++) y[scale + i] = xp[ploc[i]] - xbar[ploc[i]];
    }
    else
    {   for (i = 0; i < scale; i++) y[i] = xp[pscale[i]];
        for (i = 0; i < mlo; i++) y[i + scale] = xp[ploc[i]];
    }
    spml = scale + mlo;
    y[spml] = (float)log((double)xp[tcol]);
    y[spml + 1] = 0.;
    for (i = 0; i < nreoff; i++) if ((int)xp[cencol] == reoff[i]) y[spml + 1] = 1.;
    if ((int)xp[cencol] == rgtcen) y[spml + 1] = -1.;
    if (startcol >= 0 && xp[t0col] > 0.)
    {   y[spml + 2] = (float)(int)xp[startcol];
        if (y[spml + 2] == 1) y[spml + 3] = (float)log((double)xp[t0col]);
    }
    else y[spml + 2] = 0;
    xindex += m;
    return(y);
}

void zrewind(void) { xindex = 0; obs = 0; }


/* setupls: setup for model */
int setupls(int ncol, long nrec, int nlo, int nsc, int nth, int time,
            int censor, int *p1, int *p2, int nevent, int *event,
            int stcol, int time0)
{
    m = ncol;
    n = nrec;
    scale = nsc;
    mlo = nlo;
    ntheta = nth;
    if (ntheta == 0) mlo1 = mlo;
    else mlo1 = 1;
    spmlo = scale + mlo;
    spmlo1 = scale + mlo1;
    npar = scale + mlo + ntheta;
    pscale = p2;
    ploc = p1;
    tcol = time;
    t0col = time0;
    startcol = stcol;
    cencol = censor;
    nreoff = nevent;
    reoff = event;
    return(0);
}


/* inclvar: select subset of covariates (for backward selection) */
void inclvar(int nslct, int *in)
{
    int i, oldscale, oldmlo;
    oldscale = scale; scale = 0;
    for (i = 0; i < oldscale; i++) if (in[i] == 1) pscale[scale++] = pscale[i];
    oldmlo = mlo;
    mlo = 0;
    for (i = oldscale; i < oldscale + oldmlo; i++)
        if (in[i] == 1) ploc[mlo++] = ploc[i - oldscale];
    if (ntheta > 0)
    {   ntheta = 0;
        if (in[oldscale + oldmlo] == 1) ntheta = 1;
    }
    npar = scale + mlo + ntheta;
    Rprintf("Scale parameters remaining:");
    for (i = 0; i < scale; i++) Rprintf(" %d", pscale[i]); Rprintf("\n");
    Rprintf("Eta parameters remaining:");
    for (i = 0; i < mlo; i++) Rprintf(" %d", ploc[i]); Rprintf("\n");
    if (ntheta == 0) Rprintf("No theta\n");
    else Rprintf("Theta included\n");
}


/* data_from_r: load data from R double matrix into internal float array */
int data_from_r(double *rdata, int nobs, int ncol)
{
    int i, j;

    n = (long)nobs;
    m = ncol;

    /* allocate y (observation record) */
    if (y == NULL)
    {   y = (float *)R_alloc(MAXCOV + 4, sizeof(float));
    }

    /* allocate x (data array) */
    if ((long)nobs * ncol > MAXX)
    {   error("Data too large: n*p = %d exceeds MAXX = %d", nobs * ncol, MAXX);
        return(-1);
    }
    x = (float *)R_alloc((size_t)nobs * ncol, sizeof(float));

    /* allocate xbar */
    xbar = (double *)R_alloc(ncol, sizeof(double));
    for (j = 0; j < ncol; j++) xbar[j] = 0.;

    /* copy R column-major matrix to row-major float array and compute means */
    for (i = 0; i < nobs; i++)
    {   for (j = 0; j < ncol; j++)
        {   x[i * ncol + j] = (float)rdata[j * nobs + i];  /* R is column-major */
            xbar[j] += rdata[j * nobs + i];
        }
    }
    for (j = 0; j < ncol; j++) xbar[j] /= n;

    return(0);
}
