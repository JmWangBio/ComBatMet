#  mglmLevenbergBeta.R
#
#  Fit rowwise beta regression GLMs with logit link using
#  Levenberg-Marquardt algorithm, with phi (precision) pre-determined.

mglmLevenbergBeta <- function(y, design, phi, weights=NULL,
                              coef.start=NULL, maxit=200L, tol=1e-06)
  # y          matrix of proportions (ntags x nsamples), each in (0,1)
  # design     design matrix (nsamples x ncoef)
  # phi        matrix of dim(y): precision per tag x sample
  # weights    optional weight matrix of dim(y)
  # coef.start optional starting coefficients (ntags x ncoef)
  #
  # Returns a list: coefficients, fitted.values, deviance, iter, failed
{
  y <- as.matrix(y)
  ntags <- nrow(y)
  nsamples  <- ncol(y)
  
  design <- as.matrix(design)
  if (!is.double(design)) storage.mode(design) <- "double"
  
  ncoef <- ncol(design)
  
  phi <- as.matrix(phi)
  if (!is.double(phi)) storage.mode(phi) <- "double"

  if (!is.null(weights)) {
    weights <- as.matrix(weights)
    if (!is.double(weights)) storage.mode(weights) <- "double"
  }

  # Calculate starting coefficients
  if (is.null(coef.start)) {
    y0 <- y

    # Logit-transform data
    eta0 <- log(y0 / (1 - y0))

    # Project onto design: solve design^T design beta = design^T (eta0^T)
    coef.start <- t(qr.coef(qr(design), t(eta0)))

    # Starting beta: ntags x ncoef
    rownames(coef.start) <- rownames(y)
    colnames(coef.start) <- colnames(design)
  }

  beta <- as.matrix(coef.start)
  if (!is.double(beta)) storage.mode(beta) <- "double"

  # Call C fitter
  output <- fit_leven_beta_cpp(y, phi, weights, design, beta, tol, maxit)

  # Name and return
  colnames(output$coefficients) <- colnames(design)
  rownames(output$coefficients) <- rownames(y)
  dimnames(output$fitted.values) <- dimnames(y)
  output
}
