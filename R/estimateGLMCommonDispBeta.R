#  estimateGLMCommonDispBeta.R
#
#  Estimates a single common precision phi for a matrix of beta
#  regression rows using the Cox-Reid APL method.

estimateGLMCommonDispBeta <- function(y, design=NULL,
                                      weights=NULL,
                                      subset=10000,
                                      interval=c(1e-4, 1e4),
                                      verbose=FALSE)
  # y        matrix of proportions (ntags x nsamples), each in (0,1)
  # design   design matrix
  # weights  optional weight matrix of dim(y)
  # subset   max rows to use in estimation
  # interval search interval for phi
  # verbose  if TRUE print the estimated phi
  #
  # Returns: scalar common precision estimate
{
  y <- as.matrix(y)
  design <- as.matrix(design)

  phi <- dispCoxReidBeta(y, design=design, weights=weights,
                         interval=interval, subset=subset)

  if (verbose)
    cat(sprintf("Common phi = %.4f\n", phi))
  phi
}
