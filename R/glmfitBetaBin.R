#  glmfitBetaBin.R
#
#  Top-level GLM fitting for beta-binomial regression.

glmFitBetaBin <- function(y, n, design, phi,
                          weights = NULL, coef.start = NULL,
                          maxit = 200L, tol = 1e-06)
  # y       count matrix (ntags x nsamples): numCs
  # n       coverage matrix (ntags x nsamples)
  # design  design matrix (nsamples x ncoef)
  # phi     scalar OR matrix of dim(y): precision.
  #
  # Returns list: coefficients, fitted.values, deviance, iter, failed,
  #               df.residual, design, phi, weights, counts, coverage
{
  y <- as.matrix(y)
  n <- as.matrix(n)
  ntags <- nrow(y)
  nsamples <- ncol(y)
  design <- as.matrix(design)
  
  if (length(phi) == 1L) {
    phi_mat <- matrix(as.double(phi), ntags, nsamples)
  } else {
    phi_mat <- as.matrix(phi)
  }
  
  # Fit
  fit <- mglmLevenbergBetaBin(y, n, design = design, phi = phi_mat,
                              weights = weights, coef.start = coef.start,
                              maxit = maxit, tol = tol)
  
  # Attach metadata and return
  fit$counts      <- y
  fit$coverage    <- n
  fit$df.residual <- rep(nsamples - ncol(design), ntags)
  fit$design      <- design
  fit$phi         <- phi_mat
  fit$weights     <- weights
  fit
}
