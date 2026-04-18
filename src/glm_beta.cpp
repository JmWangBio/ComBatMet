/*
 * glm_beta.c
 *
 * Levenberg-Marquardt GLM fitting for beta regression with pre-determined phi.
 */

#include <Rcpp.h>
#include <R_ext/Lapack.h>
using namespace Rcpp;

/* ================================================================== */
/* Internal math helpers                                               */
/* ================================================================== */

/* digamma */
static double digamma_b(double x)
{
    double r = 0.0;
    while (x < 6.0) { r -= 1.0/x; x += 1.0; }
    r += log(x) - 0.5/x;
    double x2 = x*x;
    r -= (1.0/12.0)/x2;
    r += (1.0/120.0)/(x2*x2);
    r -= (1.0/252.0)/(x2*x2*x2);
    return r;
}

/* trigamma */
static double trigamma_b(double x)
{
    double r = 0.0;
    while (x < 6.0) { r += 1.0/(x*x); x += 1.0; }
    r += 1.0/x + 0.5/(x*x);
    double x2 = x*x;
    r += (1.0/6.0)/(x2*x);
    r -= (1.0/30.0)/(x2*x2*x);
    r += (1.0/42.0)/(x2*x2*x2*x);
    return r;
}

/* ================================================================== */
/* X^T W X                                                            */
/* ================================================================== */

static void compute_xtwx_b(int nlib, int ncoef, double *X, double *W, double *out)
{
    const double *xptr1 = X;
    for (int c1 = 0; c1 < ncoef; ++c1, xptr1 += nlib) {
        const double *xptr2 = X;
        for (int c2 = 0; c2 <= c1; ++c2, xptr2 += nlib) {
            out[c2] = 0.0;
            for (int lib = 0; lib < nlib; ++lib)
                out[c2] += xptr1[lib] * xptr2[lib] * W[lib];
        }
        out += ncoef;
    }
}

/* ================================================================== */
/* Beta log-likelihood (negative) for one observation                 */
/* ================================================================== */

static double compute_unit_beta_nll(double y, double mu, double phi)
{
    if (mu <= 0.0 || mu >= 1.0 || phi <= 0.0) return 1e30;
    double a = mu*phi, b = (1.0-mu)*phi;
    return -(lgamma(phi) - lgamma(a) - lgamma(b)
             + (a-1.0)*log(y) + (b-1.0)*log(1.0-y));
}

/* ================================================================== */
/* Forward pass: mu[i] = sigmoid(X[i,] . beta)                        */
/* ================================================================== */

static void autofill(int ncoef, double *beta, int nlib,
                     double *mu, double *dm)
{
    /* Use mu as temporary storage for eta */
    for (int lib = 0; lib < nlib; ++lib) mu[lib] = 0.0;

    const char   trans = 'N';
    const int    incx  = 1, incy = 1;
    const double alpha = 1.0, s_beta = 0.0;
    F77_CALL(dgemv)(&trans, &nlib, &ncoef,
                    &alpha, dm, &nlib,
                    beta,  &incx,
                    &s_beta, mu, &incy FCONE);

    for (int lib = 0; lib < nlib; ++lib)
        mu[lib] = 1.0 / (1.0 + exp(-mu[lib]));
}

/* ================================================================== */
/* fit_leven_vec_beta                                                  */
/* ================================================================== */

static void fit_leven_vec_beta(
    int nlib, double *y, double *phi_vec, double *w,
    int ncoef, double *dm, int maxit, double tol,
    double *zwpt, double *drvt, double *dl, double *db,
    double *xtwx, double *xtwc, double *nbt, double *nmu,
    double *obt, double *omu,
    double *odev, int *oiter, int *ofail)
{
    const double low_value           = 1e-10;
    const double one_millionth       = 1e-6;
    const double supremely_low_value = 1e-13;
    const double ridiculously_low    = 1e-100;
    const char   uplo = 'U';
    const int    nrhs = 1;

    /* Degenerate row check */
    int all_invalid = 1;
    for (int lib = 0; lib < nlib; ++lib)
        if (y[lib] > low_value && y[lib] < 1.0-low_value)
            { all_invalid = 0; break; }

    if (all_invalid) {
        for (int c   = 0; c   < ncoef; ++c)   obt[c]   = NA_REAL;
        for (int lib = 0; lib < nlib;  ++lib)  omu[lib] = 0.5;
        *odev = 0.0; *oiter = 0; *ofail = 0;
        return;
    }

    autofill(ncoef, obt, nlib, omu, dm);

    double dev = 0.0;
    for (int lib = 0; lib < nlib; ++lib)
        dev += w[lib] * compute_unit_beta_nll(y[lib], omu[lib], phi_vec[lib]);

    double max_info = -1.0, lambda = 0.0;
    int info = 0, iter = 0, failed = 0;

    while ((++iter) <= maxit) {

        /* Score and Fisher weights */
        for (int lib = 0; lib < nlib; ++lib) {
            double mui = omu[lib];
            double phi = phi_vec[lib];
            double dmu = mui*(1.0-mui);
            double a   = mui*phi;
            double b   = (1.0-mui)*phi;
            zwpt[lib] = phi*phi*(trigamma_b(a)+trigamma_b(b))*dmu*dmu*w[lib];
            drvt[lib] = phi*(log(y[lib])-log(1.0-y[lib])
                            -digamma_b(a)+digamma_b(b))*dmu*w[lib];
        }

        /* Fisher information */
        compute_xtwx_b(nlib, ncoef, dm, zwpt, xtwx);

        /* Score vector and max_info */
        double *dmc = dm, *xi = xtwx;
        for (int c = 0; c < ncoef; ++c, dmc += nlib, xi += ncoef) {
            dl[c] = 0.0;
            for (int lib = 0; lib < nlib; ++lib)
                dl[c] += drvt[lib]*dmc[lib];
            if (xi[c] > max_info) max_info = xi[c];
        }
        if (iter == 1) {
            lambda = max_info * one_millionth;
            if (lambda < supremely_low_value) lambda = supremely_low_value;
        }

        int lev = 0, low_dev = 0;
        while (++lev) {

            /* Copy lower triangle + augment diagonal, then Cholesky */
            do {
                double *xi2 = xtwx, *xci = xtwc;
                for (int c1 = 0; c1 < ncoef; ++c1, xi2 += ncoef, xci += ncoef) {
                    for (int c2 = 0; c2 <= c1; ++c2) xci[c2] = xi2[c2];
                    xci[c1] += lambda;
                }
                F77_CALL(dpotrf)(&uplo, &ncoef, xtwc, &ncoef, &info FCONE);
                if (info != 0) {
                    lambda *= 10.0;
                    if (lambda <= 0.0) lambda = ridiculously_low;
                } else { break; }
            } while (1);

            for (int c = 0; c < ncoef; ++c) db[c] = dl[c];
            F77_CALL(dpotrs)(&uplo, &ncoef, &nrhs, xtwc, &ncoef,
                             db, &ncoef, &info FCONE);
            if (info != 0) Rcpp::stop("dpotrs failed in fit_leven_vec_beta");

            for (int c = 0; c < ncoef; ++c) nbt[c] = obt[c]+db[c];
            autofill(ncoef, nbt, nlib, nmu, dm);

            double ndev = 0.0;
            for (int lib = 0; lib < nlib; ++lib)
                ndev += w[lib]*compute_unit_beta_nll(y[lib],nmu[lib],phi_vec[lib]);

            if (ndev < supremely_low_value) low_dev = 1;
            if (ndev <= dev || low_dev) {
                for (int c   = 0; c   < ncoef; ++c)   obt[c]   = nbt[c];
                for (int lib = 0; lib < nlib;  ++lib)  omu[lib] = nmu[lib];
                dev = ndev;
                break;
            }
            lambda *= 2.0;
            if (lambda <= 0.0) lambda = ridiculously_low;
            if (lambda/max_info > 1.0/supremely_low_value) { failed = 1; break; }
        }

        double div = 0.0;
        for (int c = 0; c < ncoef; ++c) div += dl[c]*db[c];
        if (failed || low_dev || div < tol) {
            *odev = dev; *oiter = iter; *ofail = failed;
            break;
        }
        if (lev == 1) lambda /= 10.0;
    }
}

// [[Rcpp::export]]
List fit_leven_beta_cpp(NumericMatrix y_r,
                        NumericMatrix phi_r,
                        Nullable<NumericMatrix> weights_r,
                        NumericMatrix design_r,
                        NumericMatrix beta_r,
                        double tol,
                        int maxit)
{
  int ntag  = y_r.nrow();
  int nlib  = y_r.ncol();
  int ncoef = design_r.ncol();

  NumericMatrix coef_out(ntag, ncoef);
  NumericMatrix mu_out(ntag, nlib);
  NumericVector dev_out(ntag);
  IntegerVector iter_out(ntag);
  LogicalVector fail_out(ntag);

  bool wnull = weights_r.isNull();
  NumericMatrix w_mat = wnull ? NumericMatrix(0, 0)
    : NumericMatrix(weights_r);

  // Per-row working buffers
  std::vector<double> yp(nlib), pp(nlib), wp(nlib);
  std::vector<double> bp(ncoef), up(nlib), nm(nlib), nb(ncoef);
  std::vector<double> zw(nlib), dr(nlib), dl(ncoef), db(ncoef);
  std::vector<double> xwx(ncoef*ncoef), xwc(ncoef*ncoef);

  double *dm = design_r.begin();

  for (int tag = 0; tag < ntag; ++tag) {

    // Unpack one tag row
    for (int lib = 0; lib < nlib; ++lib) {
      yp[lib] = y_r[lib*ntag + tag];
      pp[lib] = phi_r[lib*ntag + tag];
      wp[lib] = wnull ? 1.0 : w_mat[lib*ntag + tag];
    }
    for (int c = 0; c < ncoef; ++c)
      bp[c] = beta_r[c*ntag + tag];

    double od; int oi, of;
    fit_leven_vec_beta(nlib, yp.data(), pp.data(), wp.data(),
                       ncoef, dm, maxit, tol,
                       zw.data(), dr.data(), dl.data(), db.data(),
                       xwx.data(), xwc.data(),
                       nb.data(), nm.data(),
                       bp.data(), up.data(), &od, &oi, &of);

    // Write back
    for (int lib = 0; lib < nlib;  ++lib) mu_out[lib*ntag+tag]  = up[lib];
    for (int c   = 0; c   < ncoef; ++c)   coef_out[c*ntag+tag]  = bp[c];
    dev_out[tag]  = od;
    iter_out[tag] = oi;
    fail_out[tag] = of;
  }

  return List::create(
    Named("coefficients")  = coef_out,
    Named("fitted.values") = mu_out,
    Named("deviance")      = dev_out,
    Named("iter")          = iter_out,
    Named("failed")        = fail_out
  );
}
