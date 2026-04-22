/* loscreg8.c -- partials for the Anderson (1991) AFT regression model.

   ac_tpar8() computes the log-likelihood, gradient, and packed Hessian at
   parameter vector beta for the model

       log T = mu(x) + sigma(z, mu) * epsilon

   where log(sigma) = gamma' z + sum_k theta_k (mu - mean(mu))^k and the
   residual epsilon follows one of the five distributions below. */

#include "rj8def.h"

/* Set up by setupls() in dataio8.c */
extern int scale, mlo, mlo1, ntheta;
extern int (*dist)(int, double, double *);
extern double nu, enu; /* gamma shape parameter (enu = exp(nu)) */


/* ---- Distributions --------------------------------------------------
   Each returns 0 on success, -1 if the standardized residual is out of
   range. f[0] = log-density or log-survival contribution, f[1] = first
   derivative wrt the standardized residual, f[2] = second derivative. */

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
        cg  = exp(lcg);
        if ((enu <= eu && eu <= 1.) || eu < enu)
        {   /* lower incomplete gamma: pgamma gives regularized P(a,x) */
            icg = cg - pgamma(eu, enu, 1.0, 1, 0) * cg;
            if (icg <= 0.) return(-1);
        }
        else
        {   /* upper incomplete gamma: pgamma with lower.tail = FALSE */
            icg = pgamma(eu, enu, 1.0, 0, 0) * cg;
            if (icg <= 0.) return(-1);
        }
        f[0] = log(icg) - lcg;
        f[1] = -exp(x * enu - eu) / icg;
        f[2] = f[1] * (enu - eu - f[1]);
    }
    else if (x < -39.)
        return(-1);
    else
    {   lcg = lgammafn(enu);
        cg  = exp(lcg);
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
    double phi;

    if (event == 1)
    {   f[0] = dnorm(u, 0.0, 1.0, 1);  /* log phi(u) */
        f[1] = -u;
        f[2] = -1.;
        return(0);
    }
    else if (event < 0) u = -u;

    if (u < -8.)
    {   f[0] = 0.; f[1] = 0.; f[2] = 0.; return(0);
    }
    else if (u < 0.)
    {   double Phi_val = pnorm(u, 0.0, 1.0, 1, 0); /* Phi(u) */
        phi = dnorm(u, 0.0, 1.0, 0);
        f[0] = log(Phi_val);
        f[1] = -phi / Phi_val;
        f[2] = -f[1] * (u + f[1]);
        if (event < 0) f[1] = -f[1];
        return(0);
    }
    else
    {   double log_surv = pnorm(u, 0.0, 1.0, 0, 1); /* log(1 - Phi(u)) */
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
    static const double INV_PI     = 0.31830988618379067154;
    static const double LOG_PI_NEG = -1.14472988584940017414; /* -log(pi) */
    double odopx2, t, x2;
    x2 = x * x;
    odopx2 = 1. / (1. + x2);
    if (event == 1)
    {   f[0] = LOG_PI_NEG + log(odopx2);
        f[1] = -2. * x * odopx2;
        f[2] = (4. * x2 * odopx2 - 2.) * odopx2;
        return(0);
    }
    if (event < 0) x = -x;
    t = 0.5 + atan(-x) * INV_PI;
    if (t <= 0.) return(-1);
    f[0] = log(t);
    f[1] = -odopx2 * INV_PI / t;
    f[2] = f[1] * (-x * 2. * odopx2 - f[1]);
    if (event < 0) f[1] = -f[1];
    return(0);
}


/* ---- Likelihood and partials ---------------------------------------- */

int ac_tpar8(int m, double *beta, double *db, double *ddb, double *ll1)
{
    int i, ii, j, sp, sp1, event, start;
    double u, f[3] = {0., 0., 0.}, uf2pf1, tem, tem2, tem3, tem4, tem5,
           mu, mustar, sigma, odsput, dpmu, ddpmu,
           u0, f0[3] = {0., 0., 0.}, fdif[3], u0f2pf1, odspu0t, tem03;
    double *x;
    sp  = scale + mlo;
    sp1 = (sp == m) ? sp : scale + mlo1;

    /* Reset likelihood and partials. */
    *ll1 = 0.;
    j = 0;
    for (i = 0; i < m; i++)
    {   db[i] = 0.;
        for (ii = 0; ii <= i; ii++) ddb[j++] = 0.;
    }
    zrewind();

    while ((x = getrec()) != NULL)
    {
        /* mu = sum of location effects; mustar is the "coupling" sub-sum. */
        mu = 0.;
        for (i = sp1; i < sp; i++) mu += beta[i] * x[i];
        mustar = mu;
        for (i = scale; i < sp1; i++) mu += beta[i] * x[i];

        event = (x[sp + 1] == 1.) ? 1 : 0;
        start = (x[sp + 2] == 1.) ? 1 : 0;

        /* u here is temporarily log(sigma). */
        u = 0.; for (i = 0; i < scale; i++) u += beta[i] * x[i];
        dpmu = 0.; ddpmu = 0.; tem = 1.;
        if (ntheta > 0) { dpmu = beta[sp]; u += dpmu * mustar; }
        for (i = 1; i < ntheta; i++)
        {   tem2  = tem * (i + 1) * i;
            tem  *= mustar;
            dpmu  += beta[sp + i] * tem * (i + 1);
            ddpmu += beta[sp + i] * tem2;
            u     += tem * mustar * beta[sp + i];
        }
        if (u > 39 || u < -39) return(-1);
        sigma = exp(u);
        if (event == 1) *ll1 -= u;

        u = (x[sp] - mu) / sigma;
        if (0 > (*dist)((int)x[sp + 1], u, f)) return(-1);
        uf2pf1 = u * f[2] + f[1];
        for (i = 0; i < 3; i++) fdif[i] = f[i];
        u0 = 0.; u0f2pf1 = 0.;
        if (start == 1)
        {   u0 = (x[sp + 3] - mu) / sigma;
            if (0 > (*dist)(0, u0, f0)) return(-1);
            for (i = 0; i < 3; i++) fdif[i] -= f0[i];
            u0f2pf1 = (u0 * f0[2] + f0[1]);
        }
        *ll1 += fdif[0];

        /* Partials for gamma (scale) parameters. */
        j = 0;
        tem = f[1] * u - f0[1] * u0;
        if (event == 1) tem += 1.;
        tem2 = u * uf2pf1 - u0 * u0f2pf1;
        for (i = 0; i < scale; i++)
        {   db[i] -= tem * x[i];
            tem3 = x[i] * tem2;
            for (ii = 0; ii <= i; ii++) ddb[j++] -= x[ii] * tem3;
        }

        /* Partials for etas that are not coupled into sigma. */
        tem  = fdif[1] / sigma;
        tem2 = (uf2pf1 - u0f2pf1) / sigma;
        tem3 = fdif[2] / (sigma * sigma);
        for (i = scale; i < sp1; i++)
        {   db[i] -= tem * x[i];
            tem4 = tem2 * x[i];
            for (ii = 0; ii < scale; ii++) ddb[j++] -= x[ii] * tem4;
            tem4 = tem3 * x[i];
            for (ii = scale; ii <= i; ii++) ddb[j++] -= tem4 * x[ii];
        }

        /* Partials for the remaining etas (coupled into sigma via theta). */
        if (sp > sp1)
        {   odsput = 1. / sigma + dpmu * u;
            tem2   = f[2] / (sigma * sigma) + uf2pf1 * dpmu * (2. / sigma + u * dpmu)
                     - (f[1] * u + event) * ddpmu;
            tem3 = odsput * uf2pf1;
            tem  = odsput * f[1];
            odspu0t = 0.;
            if (start == 1)
            {   odspu0t = 1. / sigma + dpmu * u0;
                tem2 -= f0[2] / (sigma * sigma) + u0f2pf1 * dpmu * (2. / sigma + u0 * dpmu)
                        - f0[1] * u0 * ddpmu;
                tem3 -= odspu0t * u0f2pf1;
                tem  -= odspu0t * f0[1];
            }
            tem += dpmu * event;
            for (i = sp1; i < sp; i++)
            {   db[i] -= tem * x[i];
                tem4 = tem3 * x[i];
                for (ii = 0; ii < scale; ii++) ddb[j++] -= x[ii] * tem4;
                tem4 = x[i] * (f[2] * odsput - f0[2] * odspu0t + fdif[1] * dpmu) / sigma;
                for (ii = scale; ii < sp1; ii++) ddb[j++] -= x[ii] * tem4;
                tem4 = x[i] * tem2;
                for (ii = sp1; ii <= i; ii++) ddb[j++] -= x[ii] * tem4;
            }

            /* Partials for theta_k. */
            tem  = mustar;
            tem2 = f[1] * u - f0[1] * u0;
            tem5 = -tem2 - event;
            for (i = 0; i < ntheta; i++)
            {   db[sp + i] -= tem * (tem2 + event);
                tem3  = tem * uf2pf1;
                tem03 = tem * u0f2pf1;
                tem4  = u * tem3 - u0 * tem03;
                for (ii = 0; ii < scale; ii++) ddb[j++] -= x[ii] * tem4;
                tem4 = (tem3 - tem03) / sigma;
                for (ii = scale; ii < sp1; ii++) ddb[j++] -= x[ii] * tem4;
                tem4 = odsput * tem3 - odspu0t * tem03 + tem5;
                for (ii = sp1; ii < sp; ii++) ddb[j++] -= x[ii] * tem4;
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
