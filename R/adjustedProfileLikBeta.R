#  adjustedProfileLikBeta.R
#
#  Computes the Cox-Reid Adjusted Profile Likelihood (APL) for beta
#  regression at a pre-determined phi for all rows of y.

adjustedProfileLikBeta <- function(phi, y, design, weights=NULL,
                                   adjust=TRUE, start=NULL,
                                   get.coef=FALSE,
                                   maxit=200L, tol=1e-08)
  # phi       scalar, vector of nrow(y), or matrix of dim(y)
  # y         matrix of proportions (ntags x nsamples)
  # design    design matrix (nsamples x ncoef)
  # weights   optional weight matrix of dim(y)
  # adjust    logical: include Cox-Reid adjustment (default TRUE)
  # start     optional starting coefficients for glmFitBeta
  # get.coef  logical: if TRUE return list(apl=..., beta=...)
  #
  # Returns:
  #   get.coef=FALSE -> numeric vector of nrow(y), APL per tag
  #   get.coef=TRUE  -> list(apl=vector, beta=coefficient matrix)
{
  y <- as.matrix(y)
  ntags <- nrow(y)
  nsamples <- ncol(y)

  # Expand phi to a full ntags x nsamples matrix
  if (length(phi) == 1L) {
    phi_mat <- matrix(as.double(phi), ntags, nsamples)
  } else if (length(phi) == ntags) {
    phi_mat <- matrix(rep(as.double(phi), each=nsamples),
                      nrow=ntags, ncol=nsamples, byrow=TRUE)
  }

  # Fit conditional MLE for each row
  fit <- glmFitBeta(y, design=design, phi=phi_mat,
                    weights=weights, coef.start=start,
                    maxit=as.integer(maxit), tol=as.double(tol))
  mu <- fit$fitted.values

  adjust <- as.logical(adjust)[1]
  design <- as.matrix(design)
  if (!is.double(design)) storage.mode(design) <- "double"

  # Compute APL via C++ function
  apl <- compute_apl_beta_cpp(y, mu, phi_mat, weights, adjust, design)

  if (get.coef) {
    list(apl=apl, beta=fit$coefficients)
  } else {
    apl
  }
}
