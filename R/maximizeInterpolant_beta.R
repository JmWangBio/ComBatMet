#  maximizeInterpolant_beta.R
#
#  Finds the x-location of the maximum of a cubic spline interpolant
#  for each row of a likelihood (APL) matrix.

maximizeInterpolantBeta <- function(x, y)
  # x   numeric vector of grid points
  # y   matrix: each row is a vector of APL values at the grid points
  #
  # Returns: numeric vector of nrow(y), x-location of the
  #          maximum of the cubic spline interpolant per row
{
  y <- as.matrix(y)

  if (!is.double(x)) storage.mode(x) <- "double"
  if (!is.double(y)) storage.mode(y) <- "double"

  maximize_interpolant_cpp(x, y)
}
