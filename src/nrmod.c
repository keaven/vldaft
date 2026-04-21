#include "rj8def.h"

int priter = 0;
void setprt(void) { priter = 1 - priter; }

int prpar = 0;
void setprtpr(void) { prpar = 1 - prpar; }


/* nrsolmod: Newton-Raphson solution with modified Cholesky decomposition */

int nrsolmod(int (*par)(int, double*, double*, double*, double*), int m,
             double *finit, double *fmax, double *score,
             double *beta, double *cov, double acc, int maxiter,
             double *wk, double *delb, double *db, double *ddb, int maxhalv)
{
    int i, ii, j, k, conv;
    double fold, steep[MAXCOV], steepdel[MAXCOV], ftem;

    /* compute 1st and 2nd partials at initial beta value */
    if ((i = (*par)(m, beta, db, ddb, finit)) != 0) return(i);

    *fmax = *finit;       /* record initial obj. fn. as maximum so far */
    if (priter == 1)
    {   Rprintf("iter=0 ");
        Rprintf("ll=%9.2f ", *finit);
        Rprintf("Beta:"); for (i = 0; i < m; i++) Rprintf("%10.3e", beta[i]);
        Rprintf("\n");
    }
    if (prpar == 1)
    {   Rprintf("db:"); for (i = 0; i < m; i++) Rprintf("%10.3e", db[i]);
        Rprintf("\nddb\n");
        j = 0;
        for (i = 0; i < m; i++)
        {   for (ii = 0; ii <= i; ii++) Rprintf("%10.3e", ddb[j++]);
            Rprintf("\n");
        }
    }

    /* compute steepest ascent variables */
    for (i = 0; i < m; i++) { steep[i] = beta[i] + db[i]; steepdel[i] = db[i]; }

    /* do first step of Newton-Raphson iteration */
    conv = nr_itmod(m, beta, db, ddb, delb, cov, acc, wk);

    /* compute score statistic */
    if (conv == 2) *score = -1.; /* 2nd partials not pos. def. */
    else
    {   symmult(cov, db, wk, m, 1);
        mmult1(wk, db, score, 1, m, 1);
    }

    /* complete iterations */
    for (i = 1; i < maxiter && (conv == 0 || conv == 2); i++)
    {   fold = *fmax;     /* record obj. fn. */

        /* compute partials */
        if ((conv = (*par)(m, beta, db, ddb, fmax)) == -2) return(-2);
        ftem = fold - *fmax;
        if (priter == 1)
        {   Rprintf("iter=%d ll=%9.2f ", i, *fmax);
            Rprintf("Beta:"); for (ii = 0; ii < m; ii++) Rprintf("%10.3e", beta[ii]);
            Rprintf("\n");
        }

        /* cut delb in half if obj. fn. decreases or partials error */
        for (j = 0; j < maxhalv && (fold > *fmax || conv < 0); j++)
        {   for (k = 0; k < m; k++) { delb[k] /= 2.; beta[k] -= delb[k]; }
            if ((conv = (*par)(m, beta, db, ddb, fmax)) == -2) return(-2);
            /* if that isn't better, try steepest ascent */
            if (fold > *fmax || conv < 0)
            {   if ((conv = (*par)(m, steep, db, ddb, fmax)) == -2) return(-2);
                if (priter == 1) Rprintf("Tried steepest ascent; ");
                if (fold > *fmax || conv < 0)
                {   for (k = 0; k < m; k++)
                    {   steepdel[k] /= 2;
                        steep[k] -= steepdel[k];
                    }
                    if (priter == 1) Rprintf("it didn't work.\n");
                }
                else
                {   if (priter == 1) Rprintf("it worked.\n");
                    for (k = 0; k < m; k++) beta[k] = steep[k];
                }
            }
            if (priter == 1)
            {   Rprintf("iter=%d ll=%9.2f ", i, *fmax);
                Rprintf("Beta:"); for (ii = 0; ii < m; ii++) Rprintf("%10.3e", beta[ii]);
                Rprintf("\n");
            }
        }
        if (fold > *fmax) { wk[0] = ftem; return(-6); } /* couldn't increase likelihood */
        if (conv < 0) /* couldn't compute partials and/or likelihood */
        {   *fmax = 0.;
            j = 0;
            for (k = 0; k < m; k++)
            {   db[k] = 0.;
                for (ii = 0; ii <= k; ii++) ddb[j++] = 0.;
            }
            return(-7);
        }

        /* compute steepest ascent variables */
        for (k = 0; k < m; k++) { steep[k] = beta[k] + db[k]; steepdel[k] = db[k]; }

        /* do next step of Newton-Raphson iteration */
        conv = nr_itmod(m, beta, db, ddb, delb, cov, acc, wk);
    }
    if (conv == 0 || conv == 2) return(-4);  /* check for convergence */

    /* compute final partials, obj. fn., inverse minus 2nd partials */
    if ((*par)(m, beta, db, ddb, fmax) != 0) return(-3);
    if (chol(ddb, cov, m) != 0) return(-5);
    triinv(cov, cov, wk, m);
    trimult(cov, cov, m);
    return(i);  /* return number of iterations */
}


/* nritmod: Single Newton-Raphson iteration with modified Cholesky decomposition */

int nr_itmod(int n, double *b, double *db, double *ddb, double *delb,
             double *ddbinv, double acc, double *wk)
{
    double a;
    int k, j, conv;

    /* perform Cholesky decomposition on second partial matrix */
    conv = 1;
    if (1 == cholmod(ddb, ddbinv, delb, n)) conv = 3;

    /* use Cholesky decomposition to solve for change in b */
    lbak(ddbinv, wk, db, 1, n);
    for (j = 0; j < n; j++) wk[j] /= delb[j];
    ubak(ddbinv, delb, wk, 1, n);

    /* find inverse of ddb */
    triinv(ddbinv, ddbinv, wk, n);
    trimult(ddbinv, ddbinv, n);

    /* update b and check for convergence */
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
