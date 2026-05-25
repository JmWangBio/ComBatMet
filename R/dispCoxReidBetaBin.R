#  dispCoxReidBetaBin.R
#
#  Estimates a single common precision phi for beta-binomial regression
#  by maximising the sum of APL values across all rows.

dispCoxReidBetaBin <- function(y, n, design = NULL, weights = NULL,
                               interval = c(1e-4, 1e4),
                               tol = 1e-5, subset = 10000)
  # y        count matrix (ntags x nsamples): numCs
  # n        coverage matrix (ntags x nsamples)
  # design   design matrix
  # interval search interval for phi
  # subset   max rows to use for the subset
{
  y <- as.matrix(y)
  n <- as.matrix(n)
  design <- as.matrix(design)
  
  # Systematic subset stratified by average empirical logit proportion
  if (!is.null(subset) && nrow(y) > subset) {
    prop     <- (y + 0.5) / (n + 1)
    prop     <- pmax(pmin(prop, 1 - 1e-6), 1e-6)
    avg_logit <- rowMeans(log(prop / (1 - prop)))
    i <- systematicSubset(subset, avg_logit)
    y       <- y[i, , drop = FALSE]
    n       <- n[i, , drop = FALSE]
    weights <- if (!is.null(weights)) weights[i, , drop = FALSE] else NULL
  }
  
  fun <- function(par, y, n, design, weights) {
    phi_val <- par^4
    sum(adjustedProfileLikBetaBin(phi_val, y = y, n = n,
                                  design = design, weights = weights))
  }
  
  out <- stats::optimize(f        = fun,
                         interval = interval^0.25,
                         y        = y, n = n,
                         design   = design, weights = weights,
                         maximum  = TRUE, tol = tol)
  out$maximum^4
}
