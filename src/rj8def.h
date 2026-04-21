#ifndef RJ8DEF_H
#define RJ8DEF_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>

#define MAXX 3500000
#define MAXSS 210
#define MAXCOV 30
#define MAXREC 40000
#define MACHEPS 3.e-39

/* Function prototypes - MATRIX.C */
void triinv(double *x, double *y, double *work, int n);
void trimult(double *y, double *z, int n);
int chol(double *s, double *t, int n);
int cholmod(double *s, double *t, double *d, int n);
int ubak(double *t, double *z, double *y, int m, int n);
int lbak(double *t, double *z, double *y, int m, int n);
int spdinv(double *t, double *y, int n);
void symmult(double *x, double *y, double *z, int m, int n);
void mmult1(double *x, double *y, double *z, int k, int m, int n);
void sextract(int m, double *symin, int n, double *symout, int *rows);
void rextract(int m, double *symin, int nr, int nc, double *recout, int *rows, int *cols);
void rectri(int m, double *rec, double *tri);
void transpose(double *x, double *y, int m, int n);
double wald(int m, double *ddb, double *cov, int n, int *in, int *out, double *beta);

/* Function prototypes - NRMOD.C */
void setprt(void);
void setprtpr(void);
int nrsolmod(int (*par)(int, double*, double*, double*, double*), int m,
             double *finit, double *fmax, double *score,
             double *beta, double *cov, double acc, int maxiter,
             double *wk, double *delb, double *db, double *ddb, int maxhalv);
int nr_itmod(int n, double *b, double *db, double *ddb, double *delb,
             double *ddbinv, double acc, double *wk);

/* Function prototypes - LOSCREG8.C */
int weibull(int event, double x, double *f);
int gamma_dist(int event, double x, double *f);
int logistic(int event, double x, double *f);
int normal_dist(int event, double u, double *f);
int cauchy_dist(int event, double x, double *f);
int ac_tpar8(int m, double *beta, double *db, double *ddb, double *ll1);
int ll8(int m, double *beta, double *ll1);
int ac_tinfl(int m, double *beta, char *fname, double *cov);

/* Function prototypes - LOSCPAR.C */
int loscpar(int m, double *beta, double *db, double *ddb, double *ll1);

/* Function prototypes - APOLYTH.C */
int a(int m, float *x, double *beta, int *event, double *u, double *mu,
      double *sigma, double *delta, double *dmudb, double *ddeldb,
      double *d2mdb2, double *d2ddb2);

/* Function prototypes - DATAIO8.C */
float *getrec(void);
void zrewind(void);
int setupls(int ncol, long nrec, int nlo, int nsc, int nth, int time,
            int censor, int *p1, int *p2, int nevent, int *event,
            int stcol, int time0);
void inclvar(int nslct, int *in);

/* Function prototypes - BW.C */
int bw(int (*par)(int, double*, double*, double*, double*), int m, double acc,
       int maxiter, int maxhalv, double pdrop, int maxdf, int nforce,
       int *force, int *nin, int *incl, double *beta, double *cov,
       double *fmax, double *fslct, double *waldmax, double *waldslct);

#endif /* RJ8DEF_H */
