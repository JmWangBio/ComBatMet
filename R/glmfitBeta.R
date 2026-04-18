#  glmfitBeta.R
#
#  Top-level GLM fitting for beta regression models.

glmFitBeta <- function(y, design, phi,
                       weights=NULL, coef.start=NULL,
                       maxit=200L, tol=1e-06)
  # Fit beta regression GLM for each row of y with pre-determined phi.
  #
  # y          matrix of proportions (ntags x nsamples), each in (0,1)
  # design     design matrix (nsamples x ncoef)
  # phi        scalar OR matrix of dim(y): precision per tag x sample.
  # weights    optional weight matrix of dim(y)
  # coef.start optional starting coefficients (ntags x ncoef)
  #
  # Returns a list with elements:
  #   coefficients   ntags x ncoef  logit-scale fitted coefficients
  #   fitted.values  ntags x nsamples  fitted proportions in (0,1)
  #   deviance       length ntags  negative log-likelihood per tag
  #   iter           length ntags  iterations taken
  #   failed         length ntags  logical: TRUE if fitting failed
  #   df.residual    length ntags  nsamples - ncoef per tag
  #   bv, design, phi, weights
{
  y <- as.matrix(y)
  ntags <- nrow(y)
  nsamples <- ncol(y)
  design <- as.matrix(design)

  # Broadcast scalar to full matrix
  if (length(phi) == 1L) {
    phi <- matrix(as.double(phi), ntags, nsamples)
  } else {
    phi <- as.matrix(phi)
  }

  # Fit
  fit <- mglmLevenbergBeta(y, design=design, phi=phi,
                           weights=weights, coef.start=coef.start,
                           maxit=maxit, tol=tol)

  # Attach metadata and return
  fit$bv <- y
  fit$df.residual <- rep(nsamples - ncol(design), ntags)
  fit$design <- design
  fit$phi <- phi
  fit$weights <- weights
  fit
}
