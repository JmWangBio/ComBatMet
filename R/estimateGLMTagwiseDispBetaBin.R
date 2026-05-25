#  estimateGLMTagwiseDispBetaBin.R
#
#  Estimates tag-specific precision phi for beta-binomial regression.

estimateGLMTagwiseDispBetaBin <- function(y, n, design = NULL,
                                          phi_common,
                                          weights = NULL,
                                          grid.npts = 11,
                                          grid.range = c(-6, 6))
{
  y <- as.matrix(y)
  n <- as.matrix(n)
  ntags <- nrow(y)
  nsamples <- ncol(y)
  design <- as.matrix(design)
  
  tagwise.phi <- dispCoxReidInterpolateTagwiseBetaBin(
    y, n, design = design, phi_common = phi_common,
    weights = weights, grid.npts = grid.npts,
    grid.range = grid.range)
  
  names(tagwise.phi) <- rownames(y)
  tagwise.phi
}
