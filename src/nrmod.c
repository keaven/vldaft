/* nrmod.c -- Newton-Raphson solver with modified Cholesky safeguards. */

#include "rj8def.h"


/* nrsolmod: iterate Newton-Raphson starting from beta, returning the number of
   iterations on success or a negative error code:
     -2 : partials could not be computed
     -3 : final partials failed
     -4 : failed to converge within maxiter
     -5 : Hessian not positive-definite at the solution
     -6 : unable to increase the likelihood (all halvings/ascent failed)
     -7 : unable to compute the likelihood
*/

int nrsolmod(int (*par)(int, double*, double*, double*, double*), int m,
             double *finit, double *fmax, double *score,
             double *beta, double *cov, double acc, int maxiter,
             double *wk, double *delb, double *db, double *ddb, int maxhalv)
{
    int i, ii, j, k, conv;
    double fold, steep[MAXCOV], steepdel[MAXCOV], ftem;

    if ((i = (*par)(m, beta, db, ddb, finit)) != 0) return(i);
    *fmax = *finit;

    for (i = 0; i < m; i++) { steep[i] = beta[i] + db[i]; steepdel[i] = db[i]; }

    conv = nr_itmod(m, beta, db, ddb, delb, cov, acc, wk);

    if (conv == 2) *score = -1.; /* Hessian not PD */
    else
    {   symmult(cov, db, wk, m, 1);
        mmult1(wk, db, score, 1, m, 1);
    }

    for (i = 1; i < maxiter && (conv == 0 || conv == 2); i++)
    {   fold = *fmax;

        if ((conv = (*par)(m, beta, db, ddb, fmax)) == -2) return(-2);
        ftem = fold - *fmax;

        /* Step halving: if the likelihood dropped or partials failed, try
           halving the proposed step; if that also fails, try steepest ascent. */
        for (j = 0; j < maxhalv && (fold > *fmax || conv < 0); j++)
        {   for (k = 0; k < m; k++) { delb[k] /= 2.; beta[k] -= delb[k]; }
            if ((conv = (*par)(m, beta, db, ddb, fmax)) == -2) return(-2);
            if (fold > *fmax || conv < 0)
            {   if ((conv = (*par)(m, steep, db, ddb, fmax)) == -2) return(-2);
                if (fold > *fmax || conv < 0)
                {   for (k = 0; k < m; k++)
                    {   steepdel[k] /= 2.;
                        steep[k]    -= steepdel[k];
                    }
                }
                else
                {   for (k = 0; k < m; k++) beta[k] = steep[k];
                }
            }
        }
        if (fold > *fmax) { wk[0] = ftem; return(-6); }
        if (conv < 0)
        {   *fmax = 0.;
            j = 0;
            for (k = 0; k < m; k++)
            {   db[k] = 0.;
                for (ii = 0; ii <= k; ii++) ddb[j++] = 0.;
            }
            return(-7);
        }

        for (k = 0; k < m; k++) { steep[k] = beta[k] + db[k]; steepdel[k] = db[k]; }
        conv = nr_itmod(m, beta, db, ddb, delb, cov, acc, wk);
    }
    if (conv == 0 || conv == 2) return(-4);

    if ((*par)(m, beta, db, ddb, fmax) != 0) return(-3);
    if (chol(ddb, cov, m) != 0) return(-5);
    triinv(cov, cov, wk, m);
    trimult(cov, cov, m);
    return(i);
}


/* nr_itmod: one Newton-Raphson step using modified Cholesky to guarantee
   a descent direction even when the Hessian is indefinite. */

int nr_itmod(int n, double *b, double *db, double *ddb, double *delb,
             double *ddbinv, double acc, double *wk)
{
    double a;
    int k, j, conv;

    conv = 1;
    if (1 == cholmod(ddb, ddbinv, delb, n)) conv = 3;

    lbak(ddbinv, wk, db, 1, n);
    for (j = 0; j < n; j++) wk[j] /= delb[j];
    ubak(ddbinv, delb, wk, 1, n);

    triinv(ddbinv, ddbinv, wk, n);
    trimult(ddbinv, ddbinv, n);

    j = -1;
    for (k = 0; k < n; k++)
    {   j += k + 1;
        b[k] += delb[k];
        if (conv != 0 && conv != 2)
        {   a = sqrt(ddbinv[j]);
            a = delb[k] / a;
            if (a > acc || a < -acc) conv--;
        }
    }
    return(conv);
}
