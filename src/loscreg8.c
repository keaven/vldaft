/* Computes partials for accelerated failure time regression model with
   scale a function of covariates.
   m is the number of parameters to be estimated.
   scale is number of scale parameters to be estimated.
   mlo is the number of location parameters eta's to be estimated.
   scale+mlo must be either m-1 or m; if m-1, a theta parameter is
   estimated, relating the log of the scale parameter to the location
   parameter.
   The subroutines getrec and zrewind() must be defined by
   the user elsewhere.
 */

#include "rj8def.h"

/* work space and underlying dist declared externally */
extern int scale, mlo, mlo1, ntheta;
extern int (*dist)(int, double, double *);
extern double nu, enu; /* gamma distribution parameter (enu) and its log */


/* define distribution functions */

int weibull(int event, double x, double *f)
{
    double eu;
    if (x > 39. || x < -39.) return(-1);
    eu = exp(x);
    if (event >= 0)
    {   f[0] = -eu; f[1] = -eu; f[2] = -eu;
        if (event == 1) { f[0] += x; f[1] += 1.; }
    }
    else
    {   double emeu, omemeu;
        emeu = exp(-eu); omemeu = 1. - emeu;
        f[0] = log(omemeu);
        f[1] = eu * emeu / omemeu;
        f[2] = f[1] * (1. - eu - f[1]);
    }
    return(0);
}


int gamma_dist(int event, double x, double *f)
{
    double lcg, icg, cg;
    double eu;
    if (x > 39.)
        return(-1);
    eu = exp(x);
    if (event == 1)
    {   f[0] = x * enu - eu - lgammafn(enu);
        f[1] = enu - eu;
        f[2] = -eu;
    }
    else if (event == 0)
    {   lcg = lgammafn(enu);
        cg = exp(lcg);
        if ((enu <= eu && eu <= 1.) || eu < enu)
        {   /* lower incomplete gamma: pgamma gives regularized P(a,x) */
            icg = cg - pgamma(eu, enu, 1.0, 1, 0) * cg;
            if (icg <= 0.)
                return(-1);
        }
        else
        {   /* upper incomplete gamma: pgamma with lower.tail=FALSE */
            icg = pgamma(eu, enu, 1.0, 0, 0) * cg;
            if (icg <= 0.)
                return(-1);
        }
        f[0] = log(icg) - lcg;
        f[1] = -exp(x * enu - eu) / icg;
        f[2] = f[1] * (enu - eu - f[1]);
    }
    else if (x < -39.)
        return(-1);
    else
    {   lcg = lgammafn(enu);
        cg = exp(lcg);
        if ((enu <= eu && eu <= 1.) || eu < enu)
        {   icg = pgamma(eu, enu, 1.0, 1, 0) * cg;
            if (icg <= 0.) return(-1);
        }
        else
        {   icg = cg - pgamma(eu, enu, 1.0, 0, 0) * cg;
            if (icg <= 0.) return(-1);
        }
        f[0] = log(icg) - lcg;
        f[1] = exp(x * enu - eu) / icg;
        f[2] = f[1] * (enu - eu - f[1]);
    }
    return(0);
}


int logistic(int event, double x, double *f)
{
    double eu, opeu, opc;
    if (x > 39. || x < -39.) return(-1);
    if (event < 0) x = -x;
    eu = exp(x); opeu = 1 + eu; opc = (event == 1) + 1;
    f[0] = (event == 1) * x - opc * log(opeu);
    f[1] = opc * eu / opeu;
    f[2] = -f[1] / opeu;
    f[1] = (event == 1) - f[1];
    if (event < 0) f[1] = -f[1];
    return(0);
}


int normal_dist(int event, double u, double *f)
{
    double phi, logS;

    if (event == 1)
    {   f[0] = dnorm(u, 0.0, 1.0, 1);  /* log(phi(u)) */
        f[1] = -u;
        f[2] = -1;
        return(0);
    }
    else if (event < 0) u = -u;

    if (u < -8.)
    {   f[0] = 0.; f[1] = 0.; f[2] = 0.; return(0);
    }
    else if (u < 0.)
    {   /* Phi(-u) for u<0, i.e., Phi(|u|) which is the upper tail at u */
        /* S(u) = P(Z > u) = Phi(-u); for u<0, S(u) = Phi(|u|) > 0.5 */
        double Phi_val = pnorm(u, 0.0, 1.0, 1, 0);  /* Phi(u), lower tail */
        phi = dnorm(u, 0.0, 1.0, 0);
        f[0] = log(Phi_val);
        f[1] = -phi / Phi_val;
        f[2] = -f[1] * (u + f[1]);
        if (event < 0) f[1] = -f[1];
        return(0);
    }
    else
    {   /* u >= 0: use log-scale for numerical stability */
        double log_surv = pnorm(u, 0.0, 1.0, 0, 1);  /* log(1-Phi(u)) */
        phi = dnorm(u, 0.0, 1.0, 0);
        f[0] = log_surv;
        f[1] = -phi / exp(log_surv);
        f[2] = -f[1] * (u + f[1]);
        if (event < 0) f[1] = -f[1];
        return(0);
    }
}


int cauchy_dist(int event, double x, double *f)
{
    double odopx2;
    double t, x2;
    x2 = x * x;
    odopx2 = 1. / (1. + x2);
    if (event == 1)
    {   f[0] = -1.144729886 + log(odopx2);
        f[1] = -2. * x * odopx2;
        f[2] = (4. * x2 * odopx2 - 2.) * odopx2;
        return(0);
    }
    else
    {   if (event < 0) x = -x;
        t = .5 + atan(-x) * .318309886;
        if (t <= 0.) return(-1);
        f[0] = log(t);
        f[1] = -odopx2 * .318309886 / t;
        f[2] = f[1] * (-x * 2. * odopx2 - f[1]);
        if (event < 0) f[1] = -f[1];
        return(0);
    }
}


int ac_tpar8(int m, double *beta, double *db, double *ddb, double *ll1)
{
    int i, ii, j, spmlo, spmlo1, event, start;
    double u, f[3], uf2pf1, tem, tem2, tem3, tem4, tem5, mu, mustar, sigma, odsput, dpmu,
           ddpmu, u0, f0[3], fdif[3], u0f2pf1, odspu0t, tem03;
    float *x;
    spmlo = scale + mlo;
    if (spmlo == m) spmlo1 = spmlo;
    else spmlo1 = scale + mlo1;

    /* initialize likelihood and partials; go to top of dataset */
    *ll1 = 0.;
    j = 0;
    for (i = 0; i < m; i++)
    {   db[i] = 0.;
        for (ii = 0; ii <= i; ii++) ddb[j++] = 0.;
    }
    zrewind();

    /* loop to process data */
    while ((x = getrec()) != NULL) {

        /* compute mu */
        mu = 0.;
        for (i = spmlo1; i < spmlo; i++) mu += beta[i] * x[i];
        mustar = mu;
        for (i = scale; i < spmlo1; i++) mu += beta[i] * x[i];

        /* set event indicator */
        if (x[spmlo + 1] == 1.) event = 1; else event = 0;
        if (x[spmlo + 2] == 1.) start = 1; else start = 0;

        /* compute sigma, dpmu, ddpmu */
        u = 0.; for (i = 0; i < scale; i++) u += beta[i] * x[i];
        dpmu = 0.; ddpmu = 0.;
        if (ntheta > 0) { dpmu = beta[spmlo]; u += dpmu * mustar; ddpmu = 0.; tem = 1.; }
        for (i = 1; i < ntheta; i++)
        {   tem2 = tem * (i + 1) * i;
            tem *= mustar;
            dpmu += beta[spmlo + i] * tem * (i + 1);
            ddpmu += beta[spmlo + i] * tem2;
            u += tem * mustar * beta[spmlo + i];
        }
        if (u > 39 || u < -39)
            return(-1);
        sigma = exp(u);
        if (event == 1) *ll1 -= u; /* if event, subt. log sigma from ll */

        /* compute adjusted failure time, values of dist fn & derivatives */
        u = (x[spmlo] - mu) / sigma;
        if (0 > (*dist)((int)x[spmlo + 1], u, f))
            return(-1);
        uf2pf1 = (u * f[2] + f[1]);
        for (i = 0; i < 3; i++) fdif[i] = f[i];
        u0 = 0.; u0f2pf1 = 0.;
        if (start == 1) /* adjust for non-zero starting time */
        {   u0 = (x[spmlo + 3] - mu) / sigma;
            if (0 > (*dist)(0, u0, f0))
                return(-1);
            for (i = 0; i < 3; i++) fdif[i] -= f0[i];
            u0f2pf1 = (u0 * f0[2] + f0[1]);
        }
        *ll1 += fdif[0];

        /* compute 1st and 2nd partials for gamma's */
        j = 0;
        tem = f[1] * u - f0[1] * u0;
        if (event == 1) tem += 1.;
        tem2 = u * uf2pf1 - u0 * u0f2pf1;
        for (i = 0; i < scale; i++)
        {   db[i] -= tem * x[i];
            tem3 = x[i] * tem2;
            for (ii = 0; ii <= i; ii++) ddb[j++] -= x[ii] * tem3;
        }

        /* 1st and 2nd partials for etas not in sigma */
        tem = fdif[1] / sigma;
        tem2 = (uf2pf1 - u0f2pf1) / sigma;
        tem3 = fdif[2] / sigma / sigma;
        for (i = scale; i < spmlo1; i++)
        {   db[i] -= tem * x[i];
            tem4 = tem2 * x[i];
            for (ii = 0; ii < scale; ii++) ddb[j++] -= x[ii] * tem4;
            tem4 = tem3 * x[i];
            for (ii = scale; ii <= i; ii++) ddb[j++] -= tem4 * x[ii];
        }

        /* 1st and 2nd partials for other etas */
        if (spmlo > spmlo1)
        {   odsput = 1. / sigma + dpmu * u;
            tem2 = f[2] / sigma / sigma + uf2pf1 * dpmu * (2. / sigma + u * dpmu)
                   - (f[1] * u + event) * ddpmu;
            tem3 = odsput * uf2pf1;
            tem = odsput * f[1];
            if (start == 1)
            {   odspu0t = 1. / sigma + dpmu * u0;
                tem2 -= f0[2] / sigma / sigma + u0f2pf1 * dpmu * (2. / sigma + u0 * dpmu)
                        - f0[1] * u0 * ddpmu;
                tem3 -= odspu0t * u0f2pf1;
                tem -= odspu0t * f0[1];
            }
            else odspu0t = 0.;
            tem += dpmu * event;
            for (i = spmlo1; i < spmlo; i++)
            {   db[i] -= tem * x[i];
                tem4 = tem3 * x[i];
                for (ii = 0; ii < scale; ii++) ddb[j++] -= x[ii] * tem4;
                tem4 = x[i] * (f[2] * odsput - f0[2] * odspu0t + fdif[1] * dpmu) / sigma;
                for (ii = scale; ii < spmlo1; ii++) ddb[j++] -= x[ii] * tem4;
                tem4 = x[i] * tem2;
                for (ii = spmlo1; ii <= i; ii++) ddb[j++] -= x[ii] * tem4;
            }

            /* 1st and 2nd partials for thetas */
            tem = mustar;
            tem2 = f[1] * u - f0[1] * u0;
            tem5 = -tem2 - event;
            for (i = 0; i < ntheta; i++)
            {   db[spmlo + i] -= tem * (tem2 + event);
                tem3 = tem * uf2pf1; tem03 = tem * u0f2pf1;
                tem4 = u * tem3 - u0 * tem03;
                for (ii = 0; ii < scale; ii++) ddb[j++] -= x[ii] * tem4;
                tem4 = (tem3 - tem03) / sigma;
                for (ii = scale; ii < spmlo1; ii++) ddb[j++] -= x[ii] * tem4;
                tem4 = odsput * tem3 - odspu0t * tem03 + tem5;
                for (ii = spmlo1; ii < spmlo; ii++) ddb[j++] -= x[ii] * tem4;
                tem3 = tem * (uf2pf1 * u - u0f2pf1 * u0);
                for (ii = 0; ii <= i; ii++)
                {   tem3 *= mustar;
                    ddb[j++] -= tem3;
                }
                if (i + 1 < ntheta) { tem *= mustar; tem5 *= mustar * (i + 2); }
            }
        }
    }
    return(0);
}


int ll8(int m, double *beta, double *ll1)
{
    int i, spmlo, spmlo1, event, start;
    double u, u0, f[3], f0[3], tem, mu, mustar, sigma, dpmu;
    float *x;
    spmlo = scale + mlo;
    if (spmlo == m) spmlo1 = spmlo;
    else spmlo1 = scale + mlo1;

    /* initialize likelihood; go to top of dataset */
    *ll1 = 0.;
    zrewind();

    /* loop to process data */
    while ((x = getrec()) != NULL) {

        /* compute mu */
        mu = 0.;
        for (i = spmlo1; i < spmlo; i++) mu += beta[i] * x[i];
        mustar = mu;
        for (i = scale; i < spmlo1; i++) mu += beta[i] * x[i];

        /* set event indicator */
        if (x[spmlo + 1] == 1.) event = 1; else event = 0;
        if (x[spmlo + 2] == 1.) start = 1; else start = 0;

        /* compute sigma */
        u = 0.; for (i = 0; i < scale; i++) u += beta[i] * x[i];
        dpmu = 0.;
        if (ntheta > 0) { dpmu = beta[spmlo]; u += dpmu * mustar; tem = 1.; }
        for (i = 1; i < ntheta; i++)
        {   tem *= mustar;
            dpmu += beta[spmlo + i] * tem * (i + 1);
            u += tem * mustar * beta[spmlo + i];
        }
        if (u > 39 || u < -39) return(-1);
        sigma = exp(u);
        if (event == 1) *ll1 -= u;

        /* compute adjusted failure time, values of dist fn & derivatives */
        u = (x[spmlo] - mu) / sigma;
        if (0 > (*dist)((int)x[spmlo + 1], u, f)) return(-1);
        if (start == 1)
        {   u0 = (x[spmlo + 3] - mu) / sigma;
            if (0 > (*dist)(0, u0, f0)) return(-1);
            for (i = 0; i < 3; i++) f[i] -= f0[i];
        }
        *ll1 += f[0];
    }
    return(0);
}


int ac_tinfl(int m, double *beta, char *fname, double *cov)
{
    int i, ii, j, spmlo, spmlo1, event, start;
    double ldot[MAXCOV], u, f[3], tem, tem2, mu, mustar, sigma, odsput, dpmu,
           u0, f0[3], fdif[3], odspu0t, ival[MAXCOV], sum;
    float *x;
    FILE *fp;
    spmlo = scale + mlo;
    if (spmlo == m) spmlo1 = spmlo;
    else spmlo1 = scale + mlo1;

    /* open file for printing */
    if (NULL == (fp = fopen(fname, "w"))) return(0);

    /* go to top of dataset */
    zrewind();

    /* loop to process data */
    while ((x = getrec()) != NULL) {

        /* compute mu */
        mu = 0.;
        for (i = spmlo1; i < spmlo; i++) mu += beta[i] * x[i];
        mustar = mu;
        for (i = scale; i < spmlo1; i++) mu += beta[i] * x[i];

        /* set event indicator */
        if (x[spmlo + 1] == 1.) event = 1; else event = 0;
        if (x[spmlo + 2] == 1.) start = 1; else start = 0;

        /* compute sigma, dpmu */
        u = 0.; for (i = 0; i < scale; i++) u += beta[i] * x[i];
        dpmu = 0.;
        if (ntheta > 0) { dpmu = beta[spmlo]; u += dpmu * mustar; tem = 1.; }
        for (i = 1; i < ntheta; i++)
        {   tem2 = tem * (i + 1) * i;
            tem *= mustar;
            dpmu += beta[spmlo + i] * tem * (i + 1);
            u += tem * mustar * beta[spmlo + i];
        }
        if (u > 39 || u < -39) return(-1);
        sigma = exp(u);

        /* compute adjusted failure time, derivatives of dist'n fn */
        u = (x[spmlo] - mu) / sigma;
        if (0 > (*dist)((int)x[spmlo + 1], u, f)) return(-1);
        for (i = 0; i < 3; i++) fdif[i] = f[i];
        u0 = 0.;
        if (start == 1)
        {   u0 = (x[spmlo + 3] - mu) / sigma;
            if (0 > (*dist)(0, u0, f0)) return(-1);
            for (i = 0; i < 3; i++) fdif[i] -= f0[i];
        }

        /* compute 1st partials for scale parameters */
        tem = f[1] * u - f0[1] * u0;
        if (event == 1) tem += 1.;
        for (i = 0; i < scale; i++) ldot[i] = -tem * x[i];

        /* 1st partials for etas not in sigma */
        tem = fdif[1] / sigma;
        for (i = scale; i < spmlo1; i++) ldot[i] = -tem * x[i];

        /* 1st partials for other etas */
        if (spmlo > spmlo1)
        {   odsput = 1. / sigma + dpmu * u;
            tem = odsput * f[1];
            if (start == 1)
            {   odspu0t = 1. / sigma + dpmu * u0;
                tem -= odspu0t * f0[1];
            }
            else odspu0t = 0.;
            tem += dpmu * event;
            for (i = spmlo1; i < spmlo; i++) ldot[i] = -tem * x[i];

            /* 1st partial for thetas */
            tem = mustar;
            tem2 = f[1] * u - f0[1] * u0;
            for (i = 0; i < ntheta; i++)
            {   ldot[spmlo + i] = -tem * (tem2 + event);
                if (i + 1 < ntheta) tem *= mustar;
            }
        }

        /* print x values, influence */
        for (i = 0; i < spmlo + 2; i++)
        {   j = (int)x[i];
            if ((double)j == x[i]) fprintf(fp, "%d ", j);
            else fprintf(fp, "%7.4lf ", (double)x[i]);
        }
        fprintf(fp, "%7.4lf ", u);
        fprintf(fp, "%7.4lf ", mu);
        fprintf(fp, "%7.4lf ", sigma);
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
