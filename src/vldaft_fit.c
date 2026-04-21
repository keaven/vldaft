/* vldaft_fit.c - R-callable entry point for AFT model fitting */

#include "rj8def.h"

/* External variables from DATAIO8.C */
extern int m, scale, mlo, mlo1, ntheta, npar, madj, spmlo, spmlo1;
extern int *pscale, *ploc, *reoff, nreoff, rgtcen, startcol, t0col;
extern long n;
extern float *x, *y;
extern double *xbar;
extern int (*dist)(int, double, double *);
extern double nu, enu;

/* External from DATAIO8.C */
extern int data_from_r(double *rdata, int nobs, int ncol);

/*
 * vldaft_fit: Fit AFT regression model
 *
 * Arguments (all SEXP):
 *   s_data      - numeric matrix (n x p), data in column-major order
 *   s_time_col  - integer, 0-based column index for time variable
 *   s_event_col - integer, 0-based column index for event indicator
 *   s_event_vals- integer vector, values indicating an event
 *   s_rgtcen_val- integer, value indicating right censoring (-1 if none)
 *   s_start_col - integer, column for delayed entry indicator (-1 if none)
 *   s_t0_col    - integer, column for starting time (-1 if none)
 *   s_loc_cols  - integer vector, 0-based column indices for location covariates
 *   s_scale_cols- integer vector, 0-based column indices for scale covariates
 *   s_ntheta    - integer, number of theta parameters
 *   s_mlo1      - integer, number of location covariates that also couple to sigma via theta
 *   s_dist      - integer, distribution: 1=weibull, 2=logistic, 3=normal, 4=cauchy, 5=gamma
 *   s_nu        - numeric, gamma distribution shape parameter (only for dist=5)
 *   s_init      - numeric vector, initial parameter values (or NULL for zeros)
 *   s_acc       - numeric, convergence accuracy
 *   s_maxiter   - integer, max Newton-Raphson iterations
 *   s_maxhalv   - integer, max step halvings
 *   s_adjust    - integer, 1=center covariates, 0=don't
 *
 * Returns: list with coefficients, covariance, log-likelihood, score, convergence info
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
    double beta[MAXCOV], cov[MAXSS], wk[MAXCOV], delb[MAXCOV], db[MAXCOV], ddb[MAXSS];
    int maxiter_val, maxhalv_val;
    int dist_choice;
    SEXP result, names;

    /* Get dimensions */
    nobs = nrows(s_data);
    ncol = ncols(s_data);

    /* Load data into internal float array */
    if (data_from_r(REAL(s_data), nobs, ncol) != 0)
        error("Failed to load data");

    /* Set distribution */
    dist_choice = INTEGER(s_dist)[0];
    switch (dist_choice)
    {   case 1: dist = &weibull; break;
        case 2: dist = &logistic; break;
        case 3: dist = &normal_dist; break;
        case 4: dist = &cauchy_dist; break;
        case 5:
            dist = &gamma_dist;
            nu = REAL(s_nu)[0];
            enu = exp(nu);
            break;
        default: error("Unknown distribution: %d", dist_choice);
    }

    /* Set adjustment flag */
    madj = INTEGER(s_adjust)[0];

    /* Set event values */
    nreoff = length(s_event_vals);
    reoff = INTEGER(s_event_vals);
    rgtcen = INTEGER(s_rgtcen_val)[0];

    /* Get covariate column indices */
    nlo = length(s_loc_cols);
    nsc = length(s_scale_cols);
    nth = INTEGER(s_ntheta)[0];
    loc_cols = INTEGER(s_loc_cols);
    scale_cols = INTEGER(s_scale_cols);

    /* Set up model */
    setupls(ncol, (long)nobs, nlo, nsc, nth,
            INTEGER(s_time_col)[0], INTEGER(s_event_col)[0],
            loc_cols, scale_cols,
            nreoff, reoff,
            INTEGER(s_start_col)[0], INTEGER(s_t0_col)[0]);

    /* Override mlo1 if specified */
    if (INTEGER(s_mlo1)[0] >= 0)
    {   mlo1 = INTEGER(s_mlo1)[0];
        spmlo1 = scale + mlo1;
    }

    /* Initialize beta */
    for (i = 0; i < npar; i++) beta[i] = 0.;
    if (s_init != R_NilValue && length(s_init) > 0)
    {   int ninit = length(s_init);
        if (ninit > npar) ninit = npar;
        for (i = 0; i < ninit; i++) beta[i] = REAL(s_init)[i];
    }

    /* Get control parameters */
    acc = REAL(s_acc)[0];
    maxiter_val = INTEGER(s_maxiter)[0];
    maxhalv_val = INTEGER(s_maxhalv)[0];

    /* Fit model */
    iter = nrsolmod(ac_tpar8, npar, &finit, &fmax, &score_stat,
                    beta, cov, acc, maxiter_val, wk, delb, db, ddb, maxhalv_val);

    /* Build result list */
    PROTECT(result = allocVector(VECSXP, 8));
    PROTECT(names = allocVector(STRSXP, 8));

    /* coefficients */
    SEXP s_beta = PROTECT(allocVector(REALSXP, npar));
    for (i = 0; i < npar; i++) REAL(s_beta)[i] = beta[i];
    SET_VECTOR_ELT(result, 0, s_beta);
    SET_STRING_ELT(names, 0, mkChar("coefficients"));

    /* covariance matrix (symmetric storage -> full matrix) */
    SEXP s_cov = PROTECT(allocMatrix(REALSXP, npar, npar));
    {   int k = 0;
        for (i = 0; i < npar; i++)
        {   for (int j2 = 0; j2 <= i; j2++)
            {   REAL(s_cov)[i + j2 * npar] = cov[k];
                REAL(s_cov)[j2 + i * npar] = cov[k];
                k++;
            }
        }
    }
    SET_VECTOR_ELT(result, 1, s_cov);
    SET_STRING_ELT(names, 1, mkChar("vcov"));

    /* log-likelihood */
    SEXP s_ll = PROTECT(ScalarReal(fmax));
    SET_VECTOR_ELT(result, 2, s_ll);
    SET_STRING_ELT(names, 2, mkChar("loglik"));

    /* initial log-likelihood */
    SEXP s_llinit = PROTECT(ScalarReal(finit));
    SET_VECTOR_ELT(result, 3, s_llinit);
    SET_STRING_ELT(names, 3, mkChar("loglik_init"));

    /* score statistic */
    SEXP s_score = PROTECT(ScalarReal(score_stat));
    SET_VECTOR_ELT(result, 4, s_score);
    SET_STRING_ELT(names, 4, mkChar("score"));

    /* convergence code */
    SEXP s_conv = PROTECT(ScalarInteger(iter));
    SET_VECTOR_ELT(result, 5, s_conv);
    SET_STRING_ELT(names, 5, mkChar("iter"));

    /* number of parameters */
    SEXP s_npar = PROTECT(ScalarInteger(npar));
    SET_VECTOR_ELT(result, 6, s_npar);
    SET_STRING_ELT(names, 6, mkChar("npar"));

    /* number of observations */
    SEXP s_nobs = PROTECT(ScalarInteger(nobs));
    SET_VECTOR_ELT(result, 7, s_nobs);
    SET_STRING_ELT(names, 7, mkChar("nobs"));

    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(10);
    return(result);
}
