#  estimateGLMTagwiseDispBeta.R
#
#  Estimates a tag-specific precision parameter phi for each row of
#  a beta-value matrix using the Cox-Reid APL with grid search and
#  cubic spline interpolation.

estimateGLMTagwiseDispBeta <- function(y, design=NULL,
                                       phi_common,
                                       weights=NULL,
                                       grid.npts=11,
                                       grid.range=c(-6, 6))
  # y           matrix of proportions (ntags x nsamples), each in (0,1)
  # design      design matrix (nsamples x ncoef)
  # phi_common  scalar: common precision from estimateGLMCommonDispBeta
  # weights     optional weight matrix of dim(y)
  # grid.npts   number of APL grid points per tag
  # grid.range  log2-scale range of grid around phi_common
  #
  # Returns: numeric vector of nrow(y), one phi estimate per tag
{
  y <- as.matrix(y)
  ntags <- nrow(y)
  nsamples <- ncol(y)
  design <- as.matrix(design)

  # Delegate to the grid-search/interpolation function
  tagwise.phi <- dispCoxReidInterpolateTagwiseBeta(
    y, design=design, phi_common=phi_common,
    weights=weights, grid.npts=grid.npts,
    grid.range=grid.range)

  names(tagwise.phi) <- rownames(y)
  tagwise.phi
}
