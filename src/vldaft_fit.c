/* vldaft_fit.c -- R-callable entry point for AFT model fitting. */

#include "rj8def.h"

extern int  m, scale, mlo, mlo1, ntheta, npar, madj, spmlo, spmlo1;
extern int *pscale, *ploc, *reoff, nreoff, rgtcen, startcol, t0col;
extern long n;
extern double *x, *y, *xbar;
extern int (*dist)(int, double, double *);
extern double nu, enu;


/*
 * vldaft_fit: fit the AFT regression model.
 *
 * All arguments come from R as SEXPs:
 *   s_data        numeric matrix (n x p), column-major
 *   s_time_col    integer, 0-based index of the time column
 *   s_event_col   integer, 0-based index of the event-indicator column
 *   s_event_vals  integer vector of values meaning "event"
 *   s_rgtcen_val  integer value meaning "right-censored" (-1 if none)
 *   s_start_col   integer, column with the delayed-entry flag (-1 if none)
 *   s_t0_col      integer, column with left-truncation start time (-1 if none)
 *   s_loc_cols    integer vector of location-covariate columns
 *   s_scale_cols  integer vector of scale-covariate columns
 *   s_ntheta      integer, order of the theta polynomial
 *   s_mlo1        integer, number of location covariates NOT coupled via theta
 *                 (-1 to leave the default chosen by setupls())
 *   s_dist        integer, 1=weibull, 2=logistic, 3=normal, 4=cauchy, 5=gamma
 *   s_nu          numeric, gamma-shape parameter (only when s_dist == 5)
 *   s_init        numeric vector of initial parameter values, or NULL
 *   s_acc         numeric, convergence accuracy
 *   s_maxiter     integer, max Newton-Raphson iterations
 *   s_maxhalv     integer, max step halvings
 *   s_adjust      integer, 1 to center covariates, 0 otherwise
 *
 * Returns a named list with:
 *   coefficients, vcov, loglik, loglik_init, score, iter, npar, nobs
 */
SEXP vldaft_fit(SEXP s_data, SEXP s_time_col, SEXP s_event_col,
                SEXP s_event_vals, SEXP s_rgtcen_val,
                SEXP s_start_col, SEXP s_t0_col,
                SEXP s_loc_cols, SEXP s_scale_cols,
                SEXP s_ntheta, SEXP s_mlo1,
                SEXP s_dist, SEXP s_nu,
                SEXP s_init, SEXP s_acc, SEXP s_maxiter, SEXP s_maxhalv,
                SEXP s_adjust)
{
    int nobs, ncol, nlo, nsc, nth, iter, i;
    int *loc_cols, *scale_cols;
    double acc, finit, fmax, score_stat;
    double beta[MAXCOV], cov[MAXSS], wk[MAXCOV],
           delb[MAXCOV], db[MAXCOV], ddb[MAXSS];
    int maxiter_val, maxhalv_val, dist_choice, ninit, total_par;
    SEXP result, names;

    nobs = nrows(s_data);
    ncol = ncols(s_data);

    /* Validate parameter count against compile-time capacity. */
    nlo = length(s_loc_cols);
    nsc = length(s_scale_cols);
    nth = INTEGER(s_ntheta)[0];
    total_par = nlo + nsc + nth;
    if (total_par > MAXCOV)
        error("vldaft: too many parameters (%d > MAXCOV=%d). "
              "Reduce the model or rebuild with a larger MAXCOV.",
              total_par, MAXCOV);

    if (data_from_r(REAL(s_data), nobs, ncol) != 0)
        error("vldaft: failed to load data matrix");

    dist_choice = INTEGER(s_dist)[0];
    switch (dist_choice)
    {   case 1: dist = &weibull;      break;
        case 2: dist = &logistic;     break;
        case 3: dist = &normal_dist;  break;
        case 4: dist = &cauchy_dist;  break;
        case 5:
            dist = &gamma_dist;
            nu  = REAL(s_nu)[0];
            enu = exp(nu);
            break;
        default:
            error("vldaft: unknown distribution code %d", dist_choice);
    }

    madj   = INTEGER(s_adjust)[0];
    nreoff = length(s_event_vals);
    reoff  = INTEGER(s_event_vals);
    rgtcen = INTEGER(s_rgtcen_val)[0];

    loc_cols   = INTEGER(s_loc_cols);
    scale_cols = INTEGER(s_scale_cols);

    setupls(ncol, (long)nobs, nlo, nsc, nth,
            INTEGER(s_time_col)[0], INTEGER(s_event_col)[0],
            loc_cols, scale_cols,
            nreoff, reoff,
            INTEGER(s_start_col)[0], INTEGER(s_t0_col)[0]);

    /* Override mlo1 only if the caller asked for a specific split. */
    if (INTEGER(s_mlo1)[0] >= 0)
    {   mlo1   = INTEGER(s_mlo1)[0];
        spmlo1 = scale + mlo1;
    }

    for (i = 0; i < npar; i++) beta[i] = 0.;
    if (s_init != R_NilValue && length(s_init) > 0)
    {   ninit = length(s_init);
        if (ninit > npar) ninit = npar;
        for (i = 0; i < ninit; i++) beta[i] = REAL(s_init)[i];
    }

    acc         = REAL(s_acc)[0];
    maxiter_val = INTEGER(s_maxiter)[0];
    maxhalv_val = INTEGER(s_maxhalv)[0];

    iter = nrsolmod(ac_tpar8, npar, &finit, &fmax, &score_stat,
                    beta, cov, acc, maxiter_val,
                    wk, delb, db, ddb, maxhalv_val);

    /* Assemble the return list. */
    PROTECT(result = allocVector(VECSXP, 8));
    PROTECT(names  = allocVector(STRSXP, 8));

    {   SEXP s_beta = PROTECT(allocVector(REALSXP, npar));
        for (i = 0; i < npar; i++) REAL(s_beta)[i] = beta[i];
        SET_VECTOR_ELT(result, 0, s_beta);
        SET_STRING_ELT(names, 0, mkChar("coefficients"));

        /* Expand packed lower-triangular cov into a symmetric dense matrix. */
        SEXP s_cov = PROTECT(allocMatrix(REALSXP, npar, npar));
        {   int k = 0, j;
            for (i = 0; i < npar; i++)
            {   for (j = 0; j <= i; j++)
                {   REAL(s_cov)[i + j * npar] = cov[k];
                    REAL(s_cov)[j + i * npar] = cov[k];
                    k++;
                }
            }
        }
        SET_VECTOR_ELT(result, 1, s_cov);
        SET_STRING_ELT(names, 1, mkChar("vcov"));

        SET_VECTOR_ELT(result, 2, PROTECT(ScalarReal(fmax)));
        SET_STRING_ELT(names,  2, mkChar("loglik"));
        SET_VECTOR_ELT(result, 3, PROTECT(ScalarReal(finit)));
        SET_STRING_ELT(names,  3, mkChar("loglik_init"));
        SET_VECTOR_ELT(result, 4, PROTECT(ScalarReal(score_stat)));
        SET_STRING_ELT(names,  4, mkChar("score"));
        SET_VECTOR_ELT(result, 5, PROTECT(ScalarInteger(iter)));
        SET_STRING_ELT(names,  5, mkChar("iter"));
        SET_VECTOR_ELT(result, 6, PROTECT(ScalarInteger(npar)));
        SET_STRING_ELT(names,  6, mkChar("npar"));
        SET_VECTOR_ELT(result, 7, PROTECT(ScalarInteger(nobs)));
        SET_STRING_ELT(names,  7, mkChar("nobs"));
    }

    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(10);
    return result;
}
