#  adjustedProfileLikBetaBin.R
#
#  Computes the Cox-Reid Adjusted Profile Likelihood (APL) for
#  beta-binomial regression at a pre-determined phi for all rows.

adjustedProfileLikBetaBin <- function(phi, y, n, design,
                                      weights = NULL, adjust = TRUE,
                                      start = NULL, get.coef = FALSE,
                                      maxit = 200L, tol = 1e-06)
  # phi       scalar, vector of nrow(y), or matrix of dim(y)
  # y         count matrix (ntags x nsamples): numCs
  # n         coverage matrix (ntags x nsamples)
  # design    design matrix (nsamples x ncoef)
  # weights   optional weight matrix
  # adjust    include Cox-Reid term
  # start     optional starting coefficients
  # get.coef  if TRUE return list(apl=..., beta=...)
{
  y <- as.matrix(y)
  n <- as.matrix(n)
  ntags <- nrow(y)
  nsamples <- ncol(y)
  
  if (length(phi) == 1L) {
    phi_mat <- matrix(as.double(phi), ntags, nsamples)
  } else if (length(phi) == ntags) {
    phi_mat <- matrix(rep(as.double(phi), each = nsamples),
                      nrow = ntags, ncol = nsamples, byrow = TRUE)
  }
  
  fit <- glmFitBetaBin(y, n, design = design, phi = phi_mat,
                       weights = weights, coef.start = start,
                       maxit = as.integer(maxit), tol = as.double(tol))
  mu <- fit$fitted.values
  
  adjust <- as.logical(adjust)[1]
  design <- as.matrix(design)
  if (!is.double(design)) storage.mode(design) <- "double"
  
  apl <- compute_apl_betabin_cpp(y, n, mu, phi_mat, weights,
                                 adjust, design)
  
  if (get.coef) list(apl = apl, beta = fit$coefficients)
  else          apl
}
