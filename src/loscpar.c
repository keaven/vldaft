/* Computes partials for parametric failure time regression model with
   location and scale a function of covariates using the a() interface. */

#include "rj8def.h"

extern int (*dist)(int, double, double *);

int loscpar(int m, double *beta, double *db, double *ddb, double *ll1)
{
    int i, ii, j, event;
    double u, f[3], mu, sigma, delta, f2ds2, uf2pf1, uf2pf1ds, uuf2pf1, mf1ds;
    double mf1upe, tem1, tem2;
    double dmudb[MAXCOV], ddeldb[MAXCOV], d2mdb2[MAXSS], d2ddb2[MAXSS];
    float *x;

    /* initialize likelihood and partials; go to top of dataset */
    *ll1 = 0.;
    j = 0;
    for (i = 0; i < m; i++)
    {   db[i] = 0.; dmudb[i] = 0.; ddeldb[i] = 0.;
        for (ii = 0; ii <= i; ii++) { ddb[j] = 0.; d2mdb2[j] = 0.; d2ddb2[j++] = 0.; }
    }
    zrewind();

    /* loop to process data */
    while ((x = getrec()) != NULL) {

        /* call a() and (*dist)() to compute needed values */
        i = a(m, x, beta, &event, &u, &mu, &sigma, &delta, dmudb, ddeldb, d2mdb2, d2ddb2);
        if (i < 0)
            return(i);
        if (0 > (*dist)(event, u, f))
            return(-1);

        /* compute needed constants */
        f2ds2 = f[2] / sigma / sigma;
        uf2pf1 = u * f[2] + f[1];
        uf2pf1ds = uf2pf1 / sigma;
        uuf2pf1 = u * uf2pf1;
        mf1ds = -f[1] / sigma;
        mf1upe = -f[1] * u - event;

        /* accumulate log-likelihood and partials */
        *ll1 += f[0] - event * delta;
        j = 0;
        for (i = 0; i < m; i++)
        {   db[i] += mf1ds * dmudb[i] + mf1upe * ddeldb[i];
            tem1 = dmudb[i] * f2ds2 + ddeldb[i] * uf2pf1ds;
            tem2 = dmudb[i] * uf2pf1ds + ddeldb[i] * uuf2pf1;
            for (ii = 0; ii <= i; ii++)
            {   ddb[j] -= dmudb[ii] * tem1 + ddeldb[ii] * tem2
                          + mf1ds * d2mdb2[j] + mf1upe * d2ddb2[j];
                j++;
            }
        }
    }
    return(0);
}


int loscinfl(int m, int spmlo, double *beta, char *fname, double *cov)
{
    int i, ii, j, event;
    double u, f[3], mu, sigma, delta, mf1ds, mf1upe, ldot[MAXCOV], ival[MAXCOV], sum;
    double dmudb[MAXCOV], ddeldb[MAXCOV], d2mdb2[MAXSS], d2ddb2[MAXSS];
    float *x;
    FILE *fp;

    /* open file for printing */
    if (NULL == (fp = fopen(fname, "w"))) return(0);

    /* initialize partials; go to top of dataset */
    for (i = 0; i < m; i++) ldot[i] = 0.;
    zrewind();

    /* loop to process data */
    while ((x = getrec()) != NULL) {

        /* call a() and (*dist)() to compute needed values */
        i = a(m, x, beta, &event, &u, &mu, &sigma, &delta, dmudb, ddeldb, d2mdb2, d2ddb2);
        if (i < 0)
            return(i);
        if (0 > (*dist)(event, u, f))
            return(-1);

        /* compute needed constants */
        mf1ds = -f[1] / sigma;
        mf1upe = -f[1] * u - event;

        /* accumulate log-likelihood and partials */
        for (i = 0; i < m; i++) ldot[i] = mf1ds * dmudb[i] + mf1upe * ddeldb[i];

        /* print x values, influence */
        for (i = 0; i < spmlo + 2; i++)
        {   j = (int)x[i];
            if ((double)j == x[i]) fprintf(fp, "%d ", j);
            else fprintf(fp, "%7.4lf ", (double)x[i]);
        }
        fprintf(fp, "%7.4lf %7.4lf", mu, sigma);
        symmult(cov, ldot, ival, m, 1);
        sum = 0.;
        for (i = 0; i < m; i++)
        {   sum += ldot[i] * ival[i];
            fprintf(fp, "%12.8lf ", ival[i]);
        }
        fprintf(fp, "%12.8lf\n", sum);
    }

    /* close output file */
    fclose(fp);
    return(1);
}
