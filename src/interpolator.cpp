/*
 * interpolator.cpp
 *
 * Converted from interpolator.c (edgeR) with minimal changes.
 * fmm_spline() is inlined from fmm_spline.c (also from edgeR).
 *
 */

#include <Rcpp.h>
using namespace Rcpp;

static void fmm_spline(int n, double *x, double *y,
                       double *b, double *c, double *d)
{
    int nm1, i;
    double t;

    /* Adjustment for 1-based arrays */
    x--; y--; b--; c--; d--;

    if (n < 2) {
        return;
    }

    if (n < 3) {
        t    = (y[2] - y[1]);
        b[1] = t / (x[2] - x[1]);
        b[2] = b[1];
        c[1] = c[2] = d[1] = d[2] = 0.0;
        return;
    }

    nm1 = n - 1;

    /* Set up tridiagonal system */
    /* b = diagonal, d = offdiagonal, c = right hand side */

    d[1] = x[2] - x[1];
    c[2] = (y[2] - y[1]) / d[1];
    for (i = 2; i < n; i++) {
        d[i]   = x[i+1] - x[i];
        b[i]   = 2.0 * (d[i-1] + d[i]);
        c[i+1] = (y[i+1] - y[i]) / d[i];
        c[i]   = c[i+1] - c[i];
    }

    /* End conditions */
    b[1] = -d[1];
    b[n] = -d[nm1];
    c[1] = c[n] = 0.0;
    if (n > 3) {
        c[1] = c[3]/(x[4]-x[2]) - c[2]/(x[3]-x[1]);
        c[n] = c[nm1]/(x[n]-x[n-2]) - c[n-2]/(x[nm1]-x[n-3]);
        c[1] = c[1] * d[1] * d[1] / (x[4]-x[1]);
        c[n] = -c[n] * d[nm1] * d[nm1] / (x[n]-x[n-3]);
    }

    /* Gaussian elimination */
    for (i = 2; i <= n; i++) {
        t    = d[i-1] / b[i-1];
        b[i] = b[i] - t * d[i-1];
        c[i] = c[i] - t * c[i-1];
    }

    /* Backward substitution */
    c[n] = c[n] / b[n];
    for (i = nm1; i >= 1; i--)
        c[i] = (c[i] - d[i] * c[i+1]) / b[i];

    /* Compute polynomial coefficients */
    b[n] = (y[n] - y[n-1]) / d[n-1] + d[n-1] * (c[n-1] + 2.0*c[n]);
    for (i = 1; i <= nm1; i++) {
        b[i] = (y[i+1] - y[i]) / d[i] - d[i] * (c[i+1] + 2.0*c[i]);
        d[i] = (c[i+1] - c[i]) / d[i];
        c[i] = 3.0 * c[i];
    }
    c[n] = 3.0 * c[n];
    d[n] = d[nm1];
}

static double find_max(int npts, double *x, double *y,
                       double *b, double *c, double *d)
{
    double maxed    = -1;
    int    maxed_at = -1;

    for (int i = 0; i < npts; ++i) {
        if (maxed_at < 0 || y[i] > maxed) {
            maxed    = y[i];
            maxed_at = i;
        }
    }
    double x_max = x[maxed_at];

    fmm_spline(npts, x, y, b, c, d);

    /* Check the segment to the LEFT of the peak */
    if (maxed_at > 0) {
        double ld = d[maxed_at - 1];
        double lc = c[maxed_at - 1];
        double lb = b[maxed_at - 1];

        /* delta = c^2 - 3*d*b  (discriminant of S'(s) = 0) */
        /* fsquare(lc) from edgeR replaced by (lc)*(lc)      */
        double delta    = (lc)*(lc) - 3*ld*lb;
        int    solvable = (delta < 0) ? 0 : 1;

        if (solvable) {
            double numerator  = -lc - sqrt(delta);
            double chosen_sol = numerator / (3*ld);

            if (chosen_sol > 0 &&
                chosen_sol < x[maxed_at] - x[maxed_at - 1]) {
                double temp = ((ld*chosen_sol + lc)*chosen_sol + lb)
                              * chosen_sol + y[maxed_at - 1];
                if (temp > maxed) {
                    maxed = temp;
                    x_max = chosen_sol + x[maxed_at - 1];
                }
            }
        }
    }

    /* Check the segment to the RIGHT of the peak */
    if (maxed_at < npts - 1) {
        double rd = d[maxed_at];
        double rc = c[maxed_at];
        double rb = b[maxed_at];

        double delta    = (rc)*(rc) - 3*rd*rb;
        int    solvable = (delta < 0) ? 0 : 1;

        if (solvable) {
            double numerator  = -rc - sqrt(delta);
            double chosen_sol = numerator / (3*rd);

            if (chosen_sol > 0 &&
                chosen_sol < x[maxed_at + 1] - x[maxed_at]) {
                double temp = ((rd*chosen_sol + rc)*chosen_sol + rb)
                              * chosen_sol + y[maxed_at];
                if (temp > maxed) {
                    maxed = temp;
                    x_max = chosen_sol + x[maxed_at];
                }
            }
        }
    }

    return x_max;
}

// [[Rcpp::export]]
NumericVector maximize_interpolant_cpp(NumericVector spts,
                                        NumericMatrix ll)
{
    int npts = spts.size();
    int ntag = ll.nrow();

    if (ll.ncol() != npts)
        Rcpp::stop("length(spts) must equal ncol(ll)");

    /* Working arrays for spline coefficients, reused across rows */
    std::vector<double> lptr(npts);
    std::vector<double> b(npts), c(npts), d(npts);
    std::vector<double> x_copy(npts);

    for (int i = 0; i < npts; ++i) x_copy[i] = spts[i];

    NumericVector out(ntag);

    for (int tag = 0; tag < ntag; ++tag) {

        /* Extract row of ll (column-major: [tag, i] = i * ntag + tag) */
        for (int i = 0; i < npts; ++i)
            lptr[i] = ll[i * ntag + tag];

        out[tag] = find_max(npts,
                            x_copy.data(), lptr.data(),
                            b.data(), c.data(), d.data());
    }

    return out;
}
