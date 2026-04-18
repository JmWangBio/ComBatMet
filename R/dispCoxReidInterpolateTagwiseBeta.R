#  dispCoxReidInterpolateTagwiseBeta.R
#
#  Estimates tagwise precision phi for beta regression using the
#  Cox-Reid APL evaluated on a per-tag grid and maximized by
#  cubic spline interpolation.

dispCoxReidInterpolateTagwiseBeta <- function(y, design, phi_common,
                                              weights=NULL,
                                              grid.npts=11,
                                              grid.range=c(-6, 6),
                                              maxit.grid=50L,
                                              tol.grid=1e-4)
  # y           matrix of proportions (ntags x nsamples), each in (0,1)
  # design      design matrix (nsamples x ncoef), full column rank
  # phi_common  scalar: common precision used to center the per-tag grid.
  #             The grid spans phi_common * 2^grid.range on a log2 scale.
  # weights     optional weight matrix of dim(y)
  # grid.npts   number of APL grid points per tag
  # grid.range  log2-scale range around phi_common
  #
  # Returns: numeric vector of nrow(y), phi per tag
{
  y <- as.matrix(y)
  ntags <- nrow(y)
  nsamples <- ncol(y)

  design <- as.matrix(design)
  ncoefs <- ncol(design)

  # Grid of phi values on a log2 scale centered at phi_common.
  spline.pts <- seq(from=grid.range[1], to=grid.range[2],
                    length.out=grid.npts)

  # Center-outward traversal order
  # Find the grid point closest to 0 (i.e. closest to phi_common).
  # Then build an order that alternates outward from that center:
  #   center, center+1, center-1, center+2, center-2, ...
  center_i <- which.min(abs(spline.pts))
  left <- seq(center_i - 1, 1, by = -1)
  right <- seq(center_i + 1, grid.npts, by =  1)

  # Interleave: center, then pairs (right[1], left[1]), (right[2], left[2])...
  pairs <- as.vector(rbind(
    c(right, rep(NA, max(0, length(left) - length(right)))),
    c(left,  rep(NA, max(0, length(right) - length(left))))
  ))
  pairs <- pairs[!is.na(pairs)]
  visit_order <- c(center_i, pairs)

  # Evaluate APL at every grid point for all tags simultaneously.
  apl <- matrix(0.0, nrow=ntags, ncol=grid.npts)
  beta_cache <- vector("list", grid.npts)

  for (idx in seq_along(visit_order)) {
    i <- visit_order[idx]
    phi_i <- phi_common * 2^spline.pts[i]

    # Use the closest already-evaluated neighbor.
    warm <- NULL
    if (i > 1 && !is.null(beta_cache[[i - 1]])) warm <- beta_cache[[i - 1]]
    if (i < grid.npts && !is.null(beta_cache[[i + 1]])) {
      # If both neighbors are cached prefer the one whose phi is closer
      if (is.null(warm) ||
          abs(spline.pts[i+1] - spline.pts[i]) <
          abs(spline.pts[i-1] - spline.pts[i]))
        warm <- beta_cache[[i + 1]]
    }

    out <- adjustedProfileLikBeta(phi_i, y=y, design=design,
                                  weights=weights,
                                  start=warm,
                                  get.coef=TRUE,
                                  maxit=maxit.grid,
                                  tol=tol.grid)
    apl[, i] <- out$apl
    last.beta <- out$beta
  }

  # Tag-by-tag spline maximization.
  d <- maximizeInterpolantBeta(spline.pts, apl)

  # Convert from log2-offset back to phi scale
  phi_common * 2^d
}
