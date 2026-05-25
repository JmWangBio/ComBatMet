/*
 * glm_betabin.cpp
 *
 * Levenberg-Marquardt GLM fitting for beta-binomial regression with pre-determined phi.
 */

#include <Rcpp.h>
#include <R_ext/Lapack.h>
using namespace Rcpp;

/* ================================================================== */
/* Internal math helpers                                               */
/* ================================================================== */

static double digamma_bb(double x)
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

static double trigamma_bb(double x)
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

static void compute_xtwx_bb(int nlib, int ncoef,
                            double *X, double *W, double *out)
{
  const double *xptr1 = X;
  for (int c1 = 0; c1 < ncoef; ++c1, xptr1 += nlib) {
    const double *xptr2 = X;
    for (int c2 = 0; c2 <= c1; ++c2, xptr2 += nlib) {
      out[c2] = 0.0;
      for (int lib = 0; lib < nlib; ++lib)
        out[c2] += xptr1[lib]*xptr2[lib]*W[lib];
    }
    out += ncoef;
  }
}

static void xtwx_bb(int nlib, int ncoef,
                    const double *X, const std::vector<double> &W,
                    std::vector<double> &out)
{
  const double *xptr1 = X;
  double *op = out.data();
  for (int c1 = 0; c1 < ncoef; ++c1, xptr1 += nlib) {
    const double *xptr2 = X;
    for (int c2 = 0; c2 <= c1; ++c2, xptr2 += nlib) {
      op[c2] = 0.0;
      for (int lib = 0; lib < nlib; ++lib)
        op[c2] += xptr1[lib]*xptr2[lib]*W[lib];
    }
    op += ncoef;
  }
}

/* ================================================================== */
/* Beta-binomial log-likelihood (negative) for one observation        */
/* ================================================================== */

static double bb_nll(double y, double n, double mu, double phi)
{
  if (mu <= 0.0 || mu >= 1.0 || phi <= 0.0 || n <= 0.0) return 1e30;
  double a = mu*phi, b = (1.0-mu)*phi;
  return -(lgamma(y+a)     + lgamma(n-y+b)     - lgamma(n+a+b)
             - lgamma(a)       - lgamma(b)          + lgamma(a+b));
}

/* ================================================================== */
/* Forward pass: mu[i] = sigmoid(X[i,] . beta)                        */
/* ================================================================== */

static void autofill_bb(int ncoef, double *beta, int nlib,
                        double *mu, double *dm)
{
  for (int lib = 0; lib < nlib; ++lib) mu[lib] = 0.0;
  const char   trans = 'N';
  const int    incx  = 1, incy = 1;
  const double alpha = 1.0, sbeta = 0.0;
  F77_CALL(dgemv)(&trans, &nlib, &ncoef,
           &alpha, dm, &nlib,
           beta,  &incx,
           &sbeta, mu, &incy FCONE);
  for (int lib = 0; lib < nlib; ++lib)
    mu[lib] = 1.0/(1.0+exp(-mu[lib]));
}

/* ================================================================== */
/* fit_leven_vec_betabin                                               */
/* ================================================================== */

static void fit_leven_vec_betabin(
    int nlib, double *y, double *n_vec, double *phi_vec, double *w,
    int ncoef, double *dm, int maxit, double tol,
    double *zwpt, double *drvt, double *dl, double *db,
    double *xtwx, double *xtwc, double *nbt, double *nmu,
    double *obt, double *omu,
    double *odev, int *oiter, int *ofail)
{
  const double one_millionth       = 1e-6;
  const double supremely_low_value = 1e-13;
  const double ridiculously_low    = 1e-100;
  const char   uplo = 'U';
  const int    nrhs = 1;
  
  /* Degenerate row check */
  int all_zero = 1;
  for (int lib = 0; lib < nlib; ++lib) {
    if (n_vec[lib] > 0.0) { all_zero = 0; break; }
  }
  
  if (all_zero) {
    for (int c = 0; c < ncoef; ++c) obt[c] = NA_REAL;
    for (int lib = 0; lib < nlib; ++lib) omu[lib] = 0.5;
    *odev = 0.0; *oiter = 0; *ofail = 0;
    return;
  }
  
  autofill_bb(ncoef, obt, nlib, omu, dm);
  
  double dev = 0.0;
  for (int lib = 0; lib < nlib; ++lib)
    dev += w[lib]*bb_nll(y[lib], n_vec[lib], omu[lib], phi_vec[lib]);
  
  double max_info = -1.0, lambda = 0.0;
  int info = 0, iter = 0, failed = 0;
  
  while ((++iter) <= maxit) {
    
    /* Score and observed Fisher weights */
    for (int lib = 0; lib < nlib; ++lib) {
      double mui = omu[lib];
      double phi = phi_vec[lib];
      double ni  = n_vec[lib];
      double dmu = mui*(1.0-mui);
      double a   = mui*phi;
      double b   = (1.0-mui)*phi;
      double yi  = y[lib];
      
      drvt[lib] = phi*(digamma_bb(yi+a) - digamma_bb(ni-yi+b)
                         -digamma_bb(a)    + digamma_bb(b))
        * dmu * w[lib];
        
        double score_raw = digamma_bb(yi+a) - digamma_bb(ni-yi+b)
          - digamma_bb(a)    + digamma_bb(b);
        double curv = phi*phi
        * (trigamma_bb(a) + trigamma_bb(b)
             - trigamma_bb(yi+a) - trigamma_bb(ni-yi+b))
             * dmu*dmu;
             double link_adj = -phi * score_raw * dmu * (1.0 - 2.0*mui);
             zwpt[lib] = (curv + link_adj) * w[lib];
    }
    
    compute_xtwx_bb(nlib, ncoef, dm, zwpt, xtwx);
    
    double *dmc = dm, *xi = xtwx;
    for (int c = 0; c < ncoef; ++c, dmc += nlib, xi += ncoef) {
      dl[c] = 0.0;
      for (int lib = 0; lib < nlib; ++lib)
        dl[c] += drvt[lib]*dmc[lib];
      if (xi[c] > max_info) max_info = xi[c];
    }
    if (iter == 1) {
      lambda = max_info*one_millionth;
      if (lambda < supremely_low_value) lambda = supremely_low_value;
    }
    
    int lev = 0, low_dev = 0;
    while (++lev) {
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
        } else break;
      } while (1);
      
      for (int c = 0; c < ncoef; ++c) db[c] = dl[c];
      F77_CALL(dpotrs)(&uplo, &ncoef, &nrhs, xtwc, &ncoef,
               db, &ncoef, &info FCONE);
      if (info != 0) Rcpp::stop("dpotrs failed in fit_leven_vec_betabin");
      
      for (int c = 0; c < ncoef; ++c) nbt[c] = obt[c]+db[c];
      autofill_bb(ncoef, nbt, nlib, nmu, dm);
      
      double ndev = 0.0;
      for (int lib = 0; lib < nlib; ++lib)
        ndev += w[lib]*bb_nll(y[lib], n_vec[lib], nmu[lib], phi_vec[lib]);
      
      if (ndev < supremely_low_value) low_dev = 1;
      if (ndev <= dev || low_dev) {
        for (int c = 0; c < ncoef; ++c) obt[c] = nbt[c];
        for (int lib = 0; lib < nlib; ++lib) omu[lib] = nmu[lib];
        dev = ndev; break;
      }
      lambda *= 2.0;
      if (lambda <= 0.0) lambda = ridiculously_low;
      if (lambda/max_info > 1.0/supremely_low_value) { failed=1; break; }
    }
    
    double div = 0.0;
    for (int c = 0; c < ncoef; ++c) div += dl[c]*db[c];
    if (failed || low_dev || div < tol) {
      *odev = dev; *oiter = iter; *ofail = failed; break;
    }
    if (lev == 1) lambda /= 10.0;
  }
}

// [[Rcpp::export]]
List fit_leven_betabin_cpp(NumericMatrix y_r,
                           NumericMatrix n_r,
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
  
  bool wnull  = weights_r.isNull();
  NumericMatrix w_mat = wnull ? NumericMatrix(0,0)
    : NumericMatrix(weights_r);
  
  std::vector<double> yp(nlib), np(nlib), pp(nlib), wp(nlib);
  std::vector<double> bp(ncoef), up(nlib), nm(nlib), nb(ncoef);
  std::vector<double> zw(nlib), dr(nlib), dl(ncoef), db(ncoef);
  std::vector<double> xwx(ncoef*ncoef), xwc(ncoef*ncoef);
  
  double *dm = design_r.begin();
  
  for (int tag = 0; tag < ntag; ++tag) {
    for (int lib = 0; lib < nlib; ++lib) {
      yp[lib] = y_r[lib*ntag+tag];
      np[lib] = n_r[lib*ntag+tag];
      pp[lib] = phi_r[lib*ntag+tag];
      wp[lib] = wnull ? 1.0 : w_mat[lib*ntag+tag];
    }
    for (int c = 0; c < ncoef; ++c) bp[c] = beta_r[c*ntag+tag];
    
    double od; int oi, of;
    fit_leven_vec_betabin(nlib, yp.data(), np.data(), pp.data(), wp.data(),
                          ncoef, dm, maxit, tol,
                          zw.data(), dr.data(), dl.data(), db.data(),
                          xwx.data(), xwc.data(),
                          nb.data(), nm.data(),
                          bp.data(), up.data(), &od, &oi, &of);
    
    for (int lib = 0; lib < nlib; ++lib) mu_out[lib*ntag+tag] = up[lib];
    for (int c   = 0; c   < ncoef; ++c)  coef_out[c*ntag+tag] = bp[c];
    dev_out[tag] = od; iter_out[tag] = oi; fail_out[tag] = of;
  }
  
  return List::create(
    Named("coefficients")  = coef_out,
    Named("fitted.values") = mu_out,
    Named("deviance")      = dev_out,
    Named("iter")          = iter_out,
    Named("failed")        = fail_out
  );
}

// [[Rcpp::export]]
NumericVector compute_apl_betabin_cpp(NumericMatrix y_r,
                                      NumericMatrix n_r,
                                      NumericMatrix mu_r,
                                      NumericMatrix phi_r,
                                      Nullable<NumericMatrix> weights_r,
                                      bool do_adjust,
                                      NumericMatrix design_r)
{
  int ntag  = y_r.nrow();
  int nlib  = y_r.ncol();
  int ncoef = design_r.ncol();
  
  bool wnull = weights_r.isNull();
  NumericMatrix w_mat = wnull ? NumericMatrix(0,0)
    : NumericMatrix(weights_r);
  
  const double *dm     = design_r.begin();
  const double log_low = std::log(1e-10);
  const double low_val = 1e-10;
  
  NumericVector output(ntag, 0.0);
  
  std::vector<double> zwpt(nlib);
  std::vector<double> xtwx_buf(ncoef*ncoef);
  
  /* dsytrf workspace -- queried once */
  const char uplo = 'U';
  int info, lwork = -1;
  double tmp;
  std::vector<int>    ipiv(ncoef);
  std::vector<double> work_sytrf;
  if (do_adjust && ncoef > 1) {
    F77_CALL(dsytrf)(&uplo, &ncoef, xtwx_buf.data(), &ncoef,
             ipiv.data(), &tmp, &lwork, &info FCONE);
    lwork = (int)(tmp+0.5); if (lwork < 1) lwork = 1;
    work_sytrf.resize(lwork);
  }
  
  for (int tag = 0; tag < ntag; ++tag) {
    
    std::vector<double> yrow(nlib), nrow(nlib);
    std::vector<double> murow(nlib), phirow(nlib), wrow(nlib);
    for (int lib = 0; lib < nlib; ++lib) {
      yrow[lib]   = y_r[lib*ntag+tag];
      nrow[lib]   = n_r[lib*ntag+tag];
      murow[lib]  = mu_r[lib*ntag+tag];
      phirow[lib] = phi_r[lib*ntag+tag];
      wrow[lib]   = wnull ? 1.0 : w_mat[lib*ntag+tag];
    }
    
    double phi_tag = phirow[0];
    double apl = 0.0;
    
    /* ---------------------------------------------------------- */
    /* Beta-binomial log-likelihood and Fisher weights            */
    /* ---------------------------------------------------------- */
    for (int lib = 0; lib < nlib; ++lib) {
      double mui = murow[lib];
      double ni  = nrow[lib];
      
      if (mui <= 0.0 || mui >= 1.0 || ni <= 0.0) {
        if (do_adjust) zwpt[lib] = 0.0;
        continue;
      }
      
      double a = mui*phi_tag, b = (1.0-mui)*phi_tag;
      double yi = yrow[lib];
      
      double loglik = std::lgamma(yi+a)     + std::lgamma(ni-yi+b)
        - std::lgamma(ni+a+b)
        - std::lgamma(a)         - std::lgamma(b)
        + std::lgamma(a+b);
        apl += wrow[lib]*loglik;
        
        if (do_adjust) {
          double dmu = mui*(1.0-mui);
          double score_raw = digamma_bb(yi+a) - digamma_bb(ni-yi+b)
            - digamma_bb(a)    + digamma_bb(b);
          double curv = phi_tag*phi_tag
          * (trigamma_bb(a) + trigamma_bb(b)
               - trigamma_bb(yi+a) - trigamma_bb(ni-yi+b))
               * dmu*dmu;
               double link_adj = -phi_tag * score_raw * dmu * (1.0 - 2.0*mui);
               zwpt[lib] = (curv + link_adj) * wrow[lib];
        }
    }
    
    /* ---------------------------------------------------------- */
    /* Cox-Reid adjustment                                        */
    /* ---------------------------------------------------------- */
    if (do_adjust) {
      double adj = 0.0;
      if (ncoef == 1) {
        double s = 0.0;
        for (int lib = 0; lib < nlib; ++lib)
          s += zwpt[lib]*dm[lib]*dm[lib];
        adj = (s < low_val) ? log_low/2.0 : std::log(s)/2.0;
      } else {
        xtwx_bb(nlib, ncoef, dm, zwpt, xtwx_buf);
        F77_CALL(dsytrf)(&uplo, &ncoef, xtwx_buf.data(), &ncoef,
                 ipiv.data(), work_sytrf.data(), &lwork, &info
                           FCONE);
        if (info < 0) Rcpp::stop("LDL factorization failed for XtWX");
        for (int i = 0; i < ncoef; ++i) {
          double d = xtwx_buf[i*ncoef+i];
          adj += (std::abs(d) < 1e-10) ? log_low : std::log(std::abs(d));
        }
        adj /= 2.0;
      }
      apl -= adj;
    }
    
    output[tag] = apl;
  }
  return output;
}
