/* a() function: model parameterization for loscpar() */
#include "rj8def.h"

/* work space and underlying dist declared externally */
extern int scale, mlo, mlo1;
extern double nu, enu;
extern int spmlo, spmlo1, ntheta;

int a(int m, float *x, double *beta, int *event, double *u, double *mu,
      double *sigma, double *delta, double *dmudb, double *ddeldb,
      double *d2mdb2, double *d2ddb2)
{
    int i, ii, iii, j, jj, jjj;
    double mustar, tem, tem2, tem3, tem4, tem5, p[MAXCOV];

    /* compute mu, mustar, delta, sigma, event, u */
    *mu = 0.;
    for (i = spmlo1; i < spmlo; i++)
        *mu += beta[i] * x[i];
    mustar = *mu;
    for (i = scale; i < spmlo1; i++) *mu += beta[i] * x[i];
    *delta = 0.;
    for (i = 0; i < scale; i++)
        *delta += beta[i] * x[i];
    tem = mustar;
    for (i = 0; i < ntheta; i++)
    {   *delta += beta[spmlo + i] * tem;
        tem *= mustar;
    }
    if (*delta > 39 || *delta < -39) return(-1);
    *sigma = exp(*delta);
    *event = (int)x[spmlo + 1];
    *u = (x[spmlo] - *mu) / *sigma;

    /* ddeldb for gamma's */
    for (i = 0; i < scale; i++) ddeldb[i] = x[i];

    /* dmudb for etas */
    for (i = scale; i < spmlo; i++) dmudb[i] = x[i];

    /* compute powers of mustar */
    if (ntheta > 0)
    {   tem = mustar;
        tem2 = 1.;
        tem3 = 0.;
        iii = (spmlo1 * (spmlo1 + 1)) / 2;  /* eta cross partials for delta */
        for (i = spmlo1; i < spmlo; i++)
        {   ddeldb[i] = 0.;
            iii += spmlo1;
            for (ii = spmlo1; ii <= i; ii++) d2ddb2[iii++] = 0.;
        }
        for (i = 0; i < ntheta; i++)
        {   j = i + spmlo;
            ddeldb[j] = tem;
            tem4 = tem2 * beta[j];
            iii = (spmlo1 * (spmlo1 + 1)) / 2;  /* eta cross partials for delta */
            jjj = spmlo + i;
            jjj = ((jjj + 1) * jjj) / 2 + spmlo1; /* eta theta cross partials delta */
            for (ii = spmlo1; ii < spmlo; ii++)
            {   ddeldb[ii] += x[ii] * tem4;
                d2ddb2[jjj++] = tem2 * x[ii];
                iii += spmlo1;
                if (i > 0)
                {   tem5 = tem3 * x[ii] * beta[j];
                    for (jj = spmlo1; jj <= ii; jj++) d2ddb2[iii++] += tem5 * x[jj];
                }
            }
            if (i < ntheta - 1)
            {   tem3 = tem2 * (i + 2);
                tem2 = tem * (i + 2);
                tem *= mustar;
            }
        }
    }
    return(0);
}
