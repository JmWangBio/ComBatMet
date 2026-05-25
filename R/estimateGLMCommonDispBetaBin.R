#  estimateGLMCommonDispBetaBin.R
#
#  Estimates a single common precision phi for beta-binomial regression.

estimateGLMCommonDispBetaBin <- function(y, n, design = NULL,
                                         weights = NULL,
                                         subset = 10000,
                                         interval = c(1e-4, 1e4),
                                         verbose = FALSE)
{
  y <- as.matrix(y)
  n <- as.matrix(n)
  design <- as.matrix(design)
  
  phi <- dispCoxReidBetaBin(y, n, design = design, weights = weights,
                            interval = interval, subset = subset)
  
  if (verbose)
    cat(sprintf("Common phi = %.4f\n", phi))
  phi
}
