/*
 * compute_apl_beta.cpp
 *
 * Computes the Cox-Reid Adjusted Profile Likelihood (APL) for beta
 * regression.
 */

#include <Rcpp.h>
#include <R_ext/Lapack.h>
using namespace Rcpp;

/* ================================================================== */
/* Special functions                                                   */
/* ================================================================== */

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

static void xtwx_beta(int nlib, int ncoef,
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
        op[c2] += xptr1[lib] * xptr2[lib] * W[lib];
    }
    op += ncoef;
  }
}

// [[Rcpp::export]]
NumericVector compute_apl_beta_cpp(NumericMatrix y_r,
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

  const double *dm = design_r.begin();
  const double log_low = std::log(1e-10);
  const double low_val = 1e-10;

  NumericVector output(ntag, 0.0);

  /* Working buffers reused across rows */
  std::vector<double> zwpt(nlib);
  std::vector<double> xtwx_buf(ncoef * ncoef);

  /* dsytrf workspace -- queried once, reused for every tag */
  const char uplo = 'U';
  int info, lwork = -1;
  double tmp;
  std::vector<int>    ipiv(ncoef);
  std::vector<double> work_sytrf;
  if (do_adjust && ncoef > 1) {
    F77_CALL(dsytrf)(&uplo, &ncoef, xtwx_buf.data(), &ncoef,
                     ipiv.data(), &tmp, &lwork, &info FCONE);
    lwork = (int)(tmp + 0.5);
    if (lwork < 1) lwork = 1;
    work_sytrf.resize(lwork);
  }

  for (int tag = 0; tag < ntag; ++tag) {

    /* Extract row vectors (column-major: [row,col] = col*ntag+row) */
    std::vector<double> yrow(nlib), murow(nlib);
    std::vector<double> phirow(nlib), wrow(nlib);
    for (int lib = 0; lib < nlib; ++lib) {
      yrow[lib]  = y_r[lib*ntag  + tag];
      murow[lib] = mu_r[lib*ntag + tag];
      phirow[lib]= phi_r[lib*ntag+ tag];
      wrow[lib]  = wnull ? 1.0 : w_mat[lib*ntag + tag];
    }

    /* phi is constant across samples for one tag */
    double phi_tag = phirow[0];

    double apl = 0.0;

    /* ---------------------------------------------------------- */
    /* Beta log-likelihood                                         */
    /* ---------------------------------------------------------- */
    for (int lib = 0; lib < nlib; ++lib) {
      double mui = murow[lib];

      /* Skip degenerate fitted values */
      if (mui <= 0.0 || mui >= 1.0) {
        if (do_adjust) zwpt[lib] = 0.0;
        continue;
      }

      double a = mui * phi_tag;
      double b = (1.0 - mui) * phi_tag;

      double loglik = std::lgamma(phi_tag)
        - std::lgamma(a) - std::lgamma(b)
        + (a - 1.0) * std::log(yrow[lib])
        + (b - 1.0) * std::log(1.0 - yrow[lib]);
      apl += wrow[lib] * loglik;

      /* Fisher weight for Cox-Reid term */
      if (do_adjust) {
        double dmu = mui * (1.0 - mui);
        zwpt[lib]  = wrow[lib]
          * phi_tag * phi_tag
          * (trigamma_b(a) + trigamma_b(b))
          * dmu * dmu;
      }
    }

    /* ---------------------------------------------------------- */
    /* Cox-Reid adjustment: -(1/2) * log |det I(beta_hat, phi)|   */
    /* ---------------------------------------------------------- */
    if (do_adjust) {
      double adj = 0.0;

      if (ncoef == 1) {
        /* Scalar case: det = sum_i W_ii * X_i^2 */
        double s = 0.0;
        for (int lib = 0; lib < nlib; ++lib)
          s += zwpt[lib] * dm[lib] * dm[lib];
        adj = (s < low_val) ? log_low / 2.0
          : std::log(s) / 2.0;
      } else {
        xtwx_beta(nlib, ncoef, dm, zwpt, xtwx_buf);
        F77_CALL(dsytrf)(&uplo, &ncoef, xtwx_buf.data(), &ncoef,
                         ipiv.data(), work_sytrf.data(), &lwork, &info FCONE);
        if (info < 0) Rcpp::stop("LDL factorization failed for XtWX");
        for (int i = 0; i < ncoef; ++i) {
          double d = xtwx_buf[i*ncoef + i];
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
