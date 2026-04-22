#ifndef RJ8DEF_H
#define RJ8DEF_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

/* Upper bounds for the Newton-Raphson work arrays; also enforced in R/vldaft.R.
   - MAXCOV: maximum number of parameters (scale + location + theta)
   - MAXSS : packed lower-triangular storage size = MAXCOV * (MAXCOV + 1) / 2
   If you bump MAXCOV, keep MAXSS in sync. */
#define MAXCOV 30
#define MAXSS  ((MAXCOV * (MAXCOV + 1)) / 2)

/* Positive-definite floor used by cholmod(); matches double-precision eps. */
#define MACHEPS DBL_EPSILON

/* Function prototypes - matrix.c */
void triinv(double *x, double *y, double *work, int n);
void trimult(double *y, double *z, int n);
int  chol(double *s, double *t, int n);
int  cholmod(double *s, double *t, double *d, int n);
int  ubak(double *t, double *z, double *y, int m, int n);
int  lbak(double *t, double *z, double *y, int m, int n);
void symmult(double *x, double *y, double *z, int m, int n);
void mmult1(double *x, double *y, double *z, int k, int m, int n);

/* Function prototypes - nrmod.c */
int nrsolmod(int (*par)(int, double*, double*, double*, double*), int m,
             double *finit, double *fmax, double *score,
             double *beta, double *cov, double acc, int maxiter,
             double *wk, double *delb, double *db, double *ddb, int maxhalv);
int nr_itmod(int n, double *b, double *db, double *ddb, double *delb,
             double *ddbinv, double acc, double *wk);

/* Function prototypes - loscreg8.c */
int weibull    (int event, double x, double *f);
int gamma_dist (int event, double x, double *f);
int logistic   (int event, double x, double *f);
int normal_dist(int event, double u, double *f);
int cauchy_dist(int event, double x, double *f);
int ac_tpar8   (int m, double *beta, double *db, double *ddb, double *ll1);

/* Function prototypes - vldaft_cure_eval.c */
SEXP vldaft_cure_eval(SEXP s_data, SEXP s_time_col, SEXP s_event_col,
                      SEXP s_start_col, SEXP s_t0_col,
                      SEXP s_loc_cols, SEXP s_scale_cols, SEXP s_cure_cols,
                      SEXP s_theta, SEXP s_mlo1, SEXP s_dist, SEXP s_base_nu,
                      SEXP s_par, SEXP s_fixed_sigma);

/* Function prototypes - dataio8.c */
double *getrec(void);
void    zrewind(void);
int     setupls(int ncol, long nrec, int nlo, int nsc, int nth, int time,
                int censor, int *p1, int *p2, int nevent, int *event,
                int stcol, int time0);
int     data_from_r(double *rdata, int nobs, int ncol);

#endif /* RJ8DEF_H */
