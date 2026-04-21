#include "rj8def.h"

/* Forward declaration */
static int backdrop(int m, double *beta, double *ddb, double *cov,
                    int *opt, double pdrop, int *incl, int maxdf);

/* bw: backward elimination variable selection */

int bw(int (*par)(int, double*, double*, double*, double*), int m, double acc,
       int maxiter, int maxhalv, double pdrop, int maxdf, int nforce,
       int *force, int *nin, int *incl, double *beta, double *cov,
       double *fmax, double *fslct, double *waldmax, double *waldslct)
{
    int i, ii, j, iter, ndrop, out[MAXCOV], opt[MAXCOV], temincl[MAXCOV];
    double finit, score, wk[MAXCOV], delb[MAXCOV], db[MAXCOV], ddb[MAXSS];

    /* get vector of nonforced covariates */
    *nin = m;
    for (i = 0; i < m; i++) { opt[i] = 1; incl[i] = i; }
    for (i = 0; i < nforce; i++) opt[force[i]] = 0;
    j = 0; for (i = 0; i < m; i++) if (opt[i] == 1) out[j++] = i;

    /* Estimate full model */
    iter = nrsolmod(ac_tpar8, m, &finit, fmax, &score, beta, cov, acc, maxiter,
                    wk, delb, db, ddb, maxhalv);
    if (iter < 0)
    {   *waldmax = 0.;
        *waldslct = 0.;
        *fslct = 0.;
        return(iter);
    }

    /* Wald statistic for full model vs. forced covariates only */
    *waldmax = wald(m, ddb, cov, nforce, force, out, beta);

    /* Call backdrop to select covariates */
    ndrop = backdrop(m, beta, ddb, cov, opt, pdrop, temincl, maxdf);

    /* Loop to estimate subsequent models */
    *nin = m;
    while (ndrop > 0 && *nin > nforce)
    {
        /* Reduce covariate vector */
        j = 0;
        for (i = 0; i < *nin; i++)
        {   if (temincl[i] == 1)
            {   beta[j] = beta[i];
                incl[j] = incl[i];
                opt[j++] = opt[i];
            }
        }
        inclvar(*nin, temincl);
        *nin = *nin - ndrop;

        /* Estimate new model without dropped covariates */
        iter = nrsolmod(ac_tpar8, *nin, &finit, fslct, &score, beta, cov, acc, maxiter,
                        wk, delb, db, ddb, maxhalv);
        /* try again with beta=0 if halving didn't work */
        if (iter == -6 || iter == -7)
        {   for (i = 0; i < *nin; i++) beta[i] = 0.;
            iter = nrsolmod(ac_tpar8, *nin, &finit, fslct, &score, beta, cov, acc, maxiter,
                            wk, delb, db, ddb, maxhalv);
        }
        if (iter < 0) { *waldslct = 0.; return(iter); }

        /* Call backdrop to select covariates */
        ndrop = backdrop(*nin, beta, ddb, cov, opt, pdrop, temincl, maxdf);
    }

    /* compute Wald for selected model vs. forced covariates only & return */
    j = 0; ii = 0;
    for (i = 0; i < *nin; i++)
    {   if (opt[i] == 1) out[j++] = i;
        else temincl[ii++] = i;
    }
    *waldslct = wald(*nin, ddb, cov, nforce, temincl, out, beta);
    return(iter);
}


/* backdrop: determine set of variables to drop */

static int backdrop(int m, double *beta, double *ddb, double *cov,
                    int *opt, double pdrop, int *incl, int maxdf)
{
    int i, ii, df, j, nincl, indmin, ntest, ndrop, drop[MAXCOV], in[MAXCOV], temin[MAXCOV];
    double zmin, pmax, tem;

    /* Set all variables to be included in model */
    for (i = 0; i < m; i++) { in[i] = i; incl[i] = 1; }
    nincl = m;

    /* Compute smallest beta/(s.e.) for non-forced covariates */
    j = 0;
    ntest = 0;
    zmin = 10000.;
    indmin = 0;
    for (i = 0; i < m; i++)
    {   if (opt[i] == 1)
        {   tem = beta[i] / sqrt(cov[j]);
            tem = (tem < 0 ? -tem : tem);
            if (tem < zmin) { indmin = i; zmin = tem; }
            ntest++;
        }
        j += i + 2;
    }

    /* compute p-value for this covariate */
    if (ntest > 0) pmax = 2 * pnorm(zmin, 0.0, 1.0, 0, 0);
    else pmax = 0.;

    /* Loop to test whether further covariates should be dropped */
    ndrop = 0; df = 1;
    while (pmax > pdrop && ntest > 0)
    {
        /* Drop most recently selected covariate */
        drop[ndrop] = indmin;
        nincl--; ndrop++;
        for (i = 0; i < m && in[i] != indmin; i++);
        for (j = i; j < nincl; j++) in[j] = in[j + 1];
        incl[indmin] = 0;

        /* Select next variable based on smallest Wald statistic */
        j = 0; zmin = 10000.; ntest = 0;
        for (i = 0; i < m; i++)
        {   if (opt[i] == 1 && incl[i] == 1)
            {   drop[ndrop] = i;
                j = 0; for (ii = 0; ii < nincl; ii++) if (in[ii] != i) temin[j++] = in[ii];
                tem = wald(m, ddb, cov, nincl - 1, temin, drop, beta);
                ntest++;
                if (tem < zmin) { zmin = tem; indmin = i; }
            }
        }

        /* Compute Wald p-val */
        if (df < maxdf) df++;
        pmax = pchisq(zmin, (double)df, 0, 0);
    }

    return(ndrop);
}
