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

static void residual_partials(int m, double sigma, double u,
                              double *mu1, double *v1, double v2[MAXCOV][MAXCOV],
                              double *u1, double u2[MAXCOV][MAXCOV])
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        u1[i] = -mu1[i] / sigma - u * v1[i];
        for (j = 0; j < m; j++)
        {
            u2[i][j] = (mu1[i] * v1[j] + mu1[j] * v1[i]) / sigma
                       + u * v1[i] * v1[j] - u * v2[i][j];
        }
    }
}

int ac_tpar8(int m, double *beta, double *db, double *ddb, double *ll1)
{
    int i, ii, j, sp, sp1, event, start;
    double *x;
    sp  = scale + mlo;
    sp1 = (sp == m) ? sp : scale + mlo1;

    *ll1 = 0.;
    j = 0;
    for (i = 0; i < m; i++)
    {   db[i] = 0.;
        for (ii = 0; ii <= i; ii++) ddb[j++] = 0.;
    }
    zrewind();

    while ((x = getrec()) != NULL)
    {
        double mu = 0., mustar = 0., log_sigma = 0., sigma, dpmu = 0., ddpmu = 0.;
        double mu1[MAXCOV], v1[MAXCOV], u1[MAXCOV], u01[MAXCOV], u21[MAXCOV];
        double v2[MAXCOV][MAXCOV], u2[MAXCOV][MAXCOV], u02[MAXCOV][MAXCOV], u22[MAXCOV][MAXCOV];
        double f[3] = {0., 0., 0.}, f0[3] = {0., 0., 0.}, f2u[3] = {0., 0., 0.};
        double u, u0 = 0., uhi = 0.;
        double q1 = 0., q2 = 0., q11 = 0., q22 = 0., q12 = 0., explicit_v = 0.;

        for (i = 0; i < MAXCOV; i++)
        {
            mu1[i] = v1[i] = u1[i] = u01[i] = u21[i] = 0.;
            for (j = 0; j < MAXCOV; j++) v2[i][j] = u2[i][j] = u02[i][j] = u22[i][j] = 0.;
        }

        for (i = sp1; i < sp; i++) mustar += beta[i] * x[i];
        mu = mustar;
        for (i = scale; i < sp1; i++) mu += beta[i] * x[i];

        for (i = 0; i < scale; i++) log_sigma += beta[i] * x[i];
        if (ntheta > 0)
        {
            double pow_m = 1.;
            dpmu = beta[sp];
            log_sigma += beta[sp] * mustar;
            for (i = 1; i < ntheta; i++)
            {
                double prev_pow = pow_m;
                pow_m *= mustar;
                dpmu += beta[sp + i] * (i + 1) * pow_m;
                ddpmu += beta[sp + i] * (i + 1) * i * prev_pow;
                log_sigma += beta[sp + i] * pow_m * mustar;
            }
        }
        if (log_sigma > 39. || log_sigma < -39.) return(-1);
        sigma = exp(log_sigma);

        for (i = 0; i < scale; i++) v1[i] = x[i];
        for (i = scale; i < sp; i++)
        {
            mu1[i] = x[i];
            if (ntheta > 0 && i >= sp1) v1[i] = dpmu * x[i];
        }
        for (i = 0; i < ntheta; i++)
        {
            double pow_i = 1.;
            for (j = 0; j <= i; j++) pow_i *= mustar;
            v1[sp + i] = pow_i;
        }
        if (ntheta > 0 && sp > sp1)
        {
            for (i = sp1; i < sp; i++)
                for (j = sp1; j < sp; j++)
                    v2[i][j] = ddpmu * x[i] * x[j];
            for (i = sp1; i < sp; i++)
            {
                double pow_i = 1.;
                for (j = 0; j < ntheta; j++)
                {
                    double val = (j + 1) * pow_i * x[i];
                    v2[i][sp + j] = v2[sp + j][i] = val;
                    pow_i *= mustar;
                }
            }
        }

        event = (int)x[sp + 1];
        start = (x[sp + 2] == 1.) ? 1 : 0;
        u = (x[sp] - mu) / sigma;
        residual_partials(m, sigma, u, mu1, v1, v2, u1, u2);

        if (event == 2)
        {
            double a, b, den, ratio;
            if (x[sp + 4] <= x[sp]) return(-1);
            uhi = (x[sp + 4] - mu) / sigma;
            residual_partials(m, sigma, uhi, mu1, v1, v2, u21, u22);
            if (0 > (*dist)(0, u, f) || 0 > (*dist)(0, uhi, f2u)) return(-1);
            ratio = exp(f2u[0] - f[0]);
            if (ratio >= 1.) return(-1);
            *ll1 += f[0] + log1p(-ratio);
            a = exp(f[0]);
            b = exp(f2u[0]);
            den = a - b;
            if (den <= 0.) return(-1);
            q1 = a * f[1] / den;
            q2 = -b * f2u[1] / den;
            q11 = a * (f[2] + f[1] * f[1]) / den - q1 * q1;
            q22 = -b * (f2u[2] + f2u[1] * f2u[1]) / den - q2 * q2;
            q12 = a * f[1] * b * f2u[1] / (den * den);
        }
        else
        {
            if (0 > (*dist)(event, u, f)) return(-1);
            explicit_v = (event == 1) ? -1. : 0.;
            *ll1 += f[0] + explicit_v * log_sigma;
            q1 = f[1];
            q11 = f[2];
        }

        if (start == 1)
        {
            u0 = (x[sp + 3] - mu) / sigma;
            residual_partials(m, sigma, u0, mu1, v1, v2, u01, u02);
            if (0 > (*dist)(0, u0, f0)) return(-1);
            *ll1 -= f0[0];
        }

        j = 0;
        for (i = 0; i < m; i++)
        {
            double grad = q1 * u1[i] + explicit_v * v1[i];
            if (event == 2) grad += q2 * u21[i];
            if (start == 1) grad -= f0[1] * u01[i];
            db[i] += grad;

            for (ii = 0; ii <= i; ii++)
            {
                double second = q11 * u1[i] * u1[ii] + q1 * u2[i][ii]
                                + explicit_v * v2[i][ii];
                if (event == 2)
                {
                    second += q22 * u21[i] * u21[ii] + q2 * u22[i][ii]
                              + q12 * (u1[i] * u21[ii] + u21[i] * u1[ii]);
                }
                if (start == 1)
                {
                    second -= f0[2] * u01[i] * u01[ii] + f0[1] * u02[i][ii];
                }
                ddb[j++] -= second;
            }
        }
    }
    return(0);
}
