/* vldaft_cure_eval.c -- compiled log-likelihood / score evaluator for
   promotion-time cure models built on the vldaft baseline distributions. */

#include "rj8def.h"

extern double nu, enu; /* baseline gamma shape used by gamma_dist() */

static double mat_get(const double *x, int nrow, int row, int col)
{
    return x[(size_t)col * nrow + row];
}

static void fill_baseline(int (*dist_fn)(int, double, double *),
                          double w, double dist_nu,
                          double *logf, double *dlogf,
                          double *logS, double *dlogS, double *S, double *F)
{
    double f_event[3] = {0., 0., 0.};
    double f_surv[3]  = {0., 0., 0.};

    if (dist_fn == &gamma_dist)
    {
        nu = dist_nu;
        enu = exp(dist_nu);
    }

    if (dist_fn(1, w, f_event) != 0 || dist_fn(0, w, f_surv) != 0)
        error("vldaft_cure_eval: unable to evaluate baseline distribution");

    *logf  = f_event[0];
    *dlogf = f_event[1];
    *logS  = f_surv[0];
    *dlogS = f_surv[1];
    *S     = exp(*logS);
    *F     = 1. - *S;
}

static int choose_dist(int dist_code, int (**dist_fn)(int, double, double *))
{
    switch (dist_code)
    {
        case 1: *dist_fn = &weibull;     return 0;
        case 2: *dist_fn = &logistic;    return 0;
        case 3: *dist_fn = &normal_dist; return 0;
        case 4: *dist_fn = &cauchy_dist; return 0;
        case 5: *dist_fn = &gamma_dist;  return 0;
        default: return -1;
    }
}

SEXP vldaft_cure_eval(SEXP s_data, SEXP s_time_col, SEXP s_event_col,
                      SEXP s_time2_col, SEXP s_start_col, SEXP s_t0_col,
                      SEXP s_loc_cols, SEXP s_scale_cols, SEXP s_cure_cols,
                      SEXP s_theta, SEXP s_mlo1, SEXP s_dist, SEXP s_base_nu,
                      SEXP s_par, SEXP s_fixed_sigma)
{
    int nobs, nloc, nsc, ncu, theta, mlo1, dist_code, fixed_sigma;
    int time_col, event_col, time2_col, start_col, t0_col;
    int i, j;
    int *loc_cols, *scale_cols, *cure_cols;
    int (*dist_fn)(int, double, double *) = NULL;
    double *x, *par, *grad;
    double total_ll = 0.0, base_nu;
    int total_par;
    SEXP result, names, s_grad;

    nobs = nrows(s_data);
    x = REAL(s_data);
    time_col = INTEGER(s_time_col)[0];
    event_col = INTEGER(s_event_col)[0];
    time2_col = INTEGER(s_time2_col)[0];
    start_col = INTEGER(s_start_col)[0];
    t0_col = INTEGER(s_t0_col)[0];
    loc_cols = INTEGER(s_loc_cols);
    scale_cols = INTEGER(s_scale_cols);
    cure_cols = INTEGER(s_cure_cols);
    nloc = length(s_loc_cols);
    nsc = length(s_scale_cols);
    ncu = length(s_cure_cols);
    theta = INTEGER(s_theta)[0];
    mlo1 = INTEGER(s_mlo1)[0];
    dist_code = INTEGER(s_dist)[0];
    base_nu = REAL(s_base_nu)[0];
    fixed_sigma = INTEGER(s_fixed_sigma)[0];
    par = REAL(s_par);

    if (choose_dist(dist_code, &dist_fn) != 0)
        error("vldaft_cure_eval: unknown distribution code %d", dist_code);

    total_par = nloc + ncu + theta + (fixed_sigma ? 0 : nsc);
    if (length(s_par) != total_par)
        error("vldaft_cure_eval: parameter length mismatch (%d != %d)",
              length(s_par), total_par);

    PROTECT(result = allocVector(VECSXP, 2));
    PROTECT(names  = allocVector(STRSXP, 2));
    PROTECT(s_grad = allocVector(REALSXP, total_par));
    grad = REAL(s_grad);
    for (j = 0; j < total_par; j++) grad[j] = 0.0;

    for (i = 0; i < nobs; i++)
    {
        double logt = log(mat_get(x, nobs, i, time_col));
        int event = (int)mat_get(x, nobs, i, event_col);
        int has_start = (start_col >= 0) ? ((int)mat_get(x, nobs, i, start_col)) : 0;
        double t0 = (has_start && t0_col >= 0) ? mat_get(x, nobs, i, t0_col) : 0.0;
        double mu = 0.0, log_sigma_eff = 0.0, sigma_eff = 1.0;
        double cure_eta = 0.0, cure_p, tau;
        double dlogse_dmustar = 0.0;
        double theta_pows[MAXCOV];
        double xloc[MAXCOV], xsc[MAXCOV], xcu[MAXCOV];
        double mustar = 0.0, w;
        double logf, dlogf, logS, dlogS, S, F;
        double logf2, dlogf2, logS2, dlogS2, S2, F2, w2 = 0.0;
        double logf0, dlogf0, logS0, dlogS0, S0, F0;
        double cure_surv = 0.0, cure_surv2 = 0.0, left_denom = 1.0, interval_denom = 1.0;
        int p_scale_offset = 0;
        int p_loc_offset;
        int p_cure_offset;
        int p_theta_offset;
        int mlo_start = (mlo1 < 0) ? 0 : mlo1;

        if (!fixed_sigma)
        {
            for (j = 0; j < nsc; j++)
            {
                xsc[j] = mat_get(x, nobs, i, scale_cols[j]);
                log_sigma_eff += par[j] * xsc[j];
            }
            p_scale_offset = nsc;
        }
        for (j = 0; j < nloc; j++)
        {
            xloc[j] = mat_get(x, nobs, i, loc_cols[j]);
            mu += par[p_scale_offset + j] * xloc[j];
        }

        p_loc_offset = p_scale_offset;
        p_cure_offset = p_scale_offset + nloc;
        p_theta_offset = p_cure_offset + ncu;

        for (j = 0; j < ncu; j++)
        {
            xcu[j] = mat_get(x, nobs, i, cure_cols[j]);
            cure_eta += par[p_cure_offset + j] * xcu[j];
        }

        mustar = 0.0;
        if (theta > 0)
            for (j = mlo_start; j < nloc; j++) mustar += par[p_loc_offset + j] * xloc[j];

        theta_pows[0] = mustar;
        if (theta > 0)
        {
            log_sigma_eff += par[p_theta_offset] * mustar;
            dlogse_dmustar = par[p_theta_offset];
        }
        for (j = 1; j < theta; j++)
        {
            theta_pows[j] = theta_pows[j - 1] * mustar;
            log_sigma_eff += par[p_theta_offset + j] * theta_pows[j];
            dlogse_dmustar += par[p_theta_offset + j] * (j + 1) * theta_pows[j - 1];
        }

        if (!fixed_sigma)
        {
            if (log_sigma_eff > 39.0 || log_sigma_eff < -39.0)
                error("vldaft_cure_eval: log(sigma) out of range");
            sigma_eff = exp(log_sigma_eff);
        }

        w = (logt - mu) / sigma_eff;
        fill_baseline(dist_fn, w, base_nu, &logf, &dlogf, &logS, &dlogS, &S, &F);

        cure_p = 1.0 / (1.0 + exp(-cure_eta));
        if (cure_p <= DBL_EPSILON) cure_p = DBL_EPSILON;
        if (cure_p >= 1.0 - DBL_EPSILON) cure_p = 1.0 - DBL_EPSILON;
        tau = -log(cure_p);
        cure_surv = exp(-tau * F);
        if (event < 0)
        {
            left_denom = 1.0 - cure_surv;
            if (left_denom <= DBL_EPSILON) left_denom = DBL_EPSILON;
        }
        if (event == 2)
        {
            if (time2_col < 0) error("vldaft_cure_eval: interval upper time missing");
            w2 = (log(mat_get(x, nobs, i, time2_col)) - mu) / sigma_eff;
            fill_baseline(dist_fn, w2, base_nu, &logf2, &dlogf2, &logS2, &dlogS2, &S2, &F2);
            cure_surv2 = exp(-tau * F2);
            interval_denom = cure_surv - cure_surv2;
            if (interval_denom <= DBL_EPSILON) interval_denom = DBL_EPSILON;
        }
        else
        {
            logf2 = dlogf2 = logS2 = dlogS2 = 0.0;
            S2 = 1.0;
            F2 = 0.0;
        }

        if (event == 1)
            total_ll += log(tau) + logf - logt - log_sigma_eff - tau * F;
        else if (event == 2)
            total_ll += log(interval_denom);
        else if (event < 0)
            total_ll += log(left_denom);
        else
            total_ll += -tau * F;

        if (has_start && t0 > 0.0)
        {
            double w0 = (log(t0) - mu) / sigma_eff;
            fill_baseline(dist_fn, w0, base_nu, &logf0, &dlogf0, &logS0, &dlogS0, &S0, &F0);
            total_ll += tau * F0;
        }
        else
        {
            logf0 = dlogf0 = logS0 = dlogS0 = 0.0;
            S0 = 1.0;
            F0 = 0.0;
        }

        if (!fixed_sigma)
        {
            for (j = 0; j < nsc; j++)
            {
                double dw = -w * xsc[j];
                double dw2 = -w2 * xsc[j];
                double event_score = dlogf * dw - xsc[j];
                double surv_deriv = S * dlogS * dw;
                double score;
                if (event == 1)
                    score = event_score + tau * surv_deriv;
                else if (event == 2)
                    score = (cure_surv * tau * surv_deriv -
                             cure_surv2 * tau * S2 * dlogS2 * dw2) / interval_denom;
                else if (event < 0)
                    score = -cure_surv * tau * surv_deriv / left_denom;
                else
                    score = tau * surv_deriv;
                if (has_start && t0 > 0.0)
                    score -= tau * S0 * dlogS0 * (-((log(t0) - mu) / sigma_eff) * xsc[j]);
                grad[j] += score;
            }
        }

        for (j = 0; j < nloc; j++)
        {
            double dlogse = 0.0;
            double dw;
            double event_score, score;
            if (theta > 0 && j >= mlo_start)
                dlogse = dlogse_dmustar * xloc[j];
            dw = -xloc[j] / sigma_eff - w * dlogse;
            event_score = dlogf * dw - dlogse;
            {
                double dw2 = -xloc[j] / sigma_eff - w2 * dlogse;
                double surv_deriv = S * dlogS * dw;
                if (event == 1)
                    score = event_score + tau * surv_deriv;
                else if (event == 2)
                    score = (cure_surv * tau * surv_deriv -
                             cure_surv2 * tau * S2 * dlogS2 * dw2) / interval_denom;
                else if (event < 0)
                    score = -cure_surv * tau * surv_deriv / left_denom;
                else
                    score = tau * surv_deriv;
            }
            if (has_start && t0 > 0.0)
            {
                double w0 = (log(t0) - mu) / sigma_eff;
                double dw0 = -xloc[j] / sigma_eff - w0 * dlogse;
                score -= tau * S0 * dlogS0 * dw0;
            }
            grad[p_loc_offset + j] += score;
        }

        for (j = 0; j < ncu; j++)
        {
            double score;
            if (event == 1)
                score = xcu[j] * ((1.0 - cure_p) * (F - 1.0 / tau));
            else if (event == 2)
                score = xcu[j] * (cure_surv * (1.0 - cure_p) * F -
                                  cure_surv2 * (1.0 - cure_p) * F2) / interval_denom;
            else if (event < 0)
                score = -xcu[j] * (cure_surv * (1.0 - cure_p) * F / left_denom);
            else
                score = xcu[j] * ((1.0 - cure_p) * F);
            if (has_start && t0 > 0.0)
                score -= xcu[j] * ((1.0 - cure_p) * F0);
            grad[p_cure_offset + j] += score;
        }

        for (j = 0; j < theta; j++)
        {
            double mustar_pow = theta_pows[j];
            double dlogse = mustar_pow;
            double dw = -w * dlogse;
            double dw2 = -w2 * dlogse;
            double event_score = dlogf * dw - dlogse;
            double surv_deriv = S * dlogS * dw;
            double score;
            if (event == 1)
                score = event_score + tau * surv_deriv;
            else if (event == 2)
                score = (cure_surv * tau * surv_deriv -
                         cure_surv2 * tau * S2 * dlogS2 * dw2) / interval_denom;
            else if (event < 0)
                score = -cure_surv * tau * surv_deriv / left_denom;
            else
                score = tau * surv_deriv;
            if (has_start && t0 > 0.0)
            {
                double w0 = (log(t0) - mu) / sigma_eff;
                double dw0 = -w0 * dlogse;
                score -= tau * S0 * dlogS0 * dw0;
            }
            grad[p_theta_offset + j] += score;
        }
    }

    SET_VECTOR_ELT(result, 0, ScalarReal(total_ll));
    SET_STRING_ELT(names, 0, mkChar("loglik"));
    SET_VECTOR_ELT(result, 1, s_grad);
    SET_STRING_ELT(names, 1, mkChar("gradient"));
    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(3);
    return result;
}
