#  mglmLevenbergBetaBin.R
#
#  Fit rowwise beta-binomial GLMs with logit link using
#  Levenberg-Marquardt algorithm, with phi (precision) pre-determined.

mglmLevenbergBetaBin <- function(y, n, design, phi,
                                 weights = NULL, coef.start = NULL,
                                 maxit = 200L, tol = 1e-06)
  # y          count matrix (ntags x nsamples): numCs
  # n          coverage matrix (ntags x nsamples)
  # design     design matrix (nsamples x ncoef)
  # phi        matrix of dim(y): precision
  # weights    optional weight matrix of dim(y)
  # coef.start optional starting coefficients (ntags x ncoef)
  #
  # Returns list: coefficients, fitted.values, deviance, iter, failed
{
  y <- as.matrix(y)
  n <- as.matrix(n)
  ntags <- nrow(y)
  nsamples <- ncol(y)
  
  design <- as.matrix(design)
  if (!is.double(design)) storage.mode(design) <- "double"
  
  ncoef <- ncol(design)
  
  phi_mat <- as.matrix(phi)
  if (!is.double(phi_mat)) storage.mode(phi_mat) <- "double"
  
  if (!is.null(weights)) {
    weights <- as.matrix(weights)
    if (!is.double(weights)) storage.mode(weights) <- "double"
  }
  
  # Calculate starting coefficients
  if (is.null(coef.start)) {
    delta    <- 0.5     # add 0.5 to avoid exact 0 or n
    prop0    <- (y + delta) / (n + 2*delta)
    eta0     <- log(prop0 / (1 - prop0))
    coef.start <- t(qr.coef(qr(design), t(eta0)))
    rownames(coef.start) <- rownames(y)
    colnames(coef.start) <- colnames(design)
  }
  
  beta <- as.matrix(coef.start)
  if (!is.double(beta)) storage.mode(beta) <- "double"
  
  # Call C++ fitter
  output <- fit_leven_betabin_cpp(y, n, phi_mat, weights, design,
                                  beta, as.double(tol), as.integer(maxit))
  
  # Name and return
  colnames(output$coefficients)  <- colnames(design)
  rownames(output$coefficients)  <- rownames(y)
  dimnames(output$fitted.values) <- dimnames(y)
  output
}
