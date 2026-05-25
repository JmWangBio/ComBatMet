#  dispCoxReidInterpolateTagwiseBetaBin.R
#
#  Estimates tagwise precision phi for beta-binomial regression using
#  Cox-Reid APL on a per-tag grid, maximised by cubic spline interpolation.

dispCoxReidInterpolateTagwiseBetaBin <- function(y, n, design, phi_common,
                                                 weights = NULL,
                                                 grid.npts = 11,
                                                 grid.range = c(-6, 6),
                                                 maxit.grid = 50L,
                                                 tol.grid = 1e-4)
  # y           count matrix (ntags x nsamples): numCs
  # n           coverage matrix (ntags x nsamples)
  # design      design matrix (nsamples x ncoef)
  # phi_common  scalar: common precision used to center the grid
  # weights     optional weight matrix
  # grid.npts   number of APL grid points per tag
  # grid.range  log2-scale range around phi_common
{
  y <- as.matrix(y)
  n <- as.matrix(n)
  ntags <- nrow(y)
  nsamples <- ncol(y)
  
  design <- as.matrix(design)
  ncoefs <- ncol(design)
  
  spline.pts <- seq(from = grid.range[1], to = grid.range[2],
                    length.out = grid.npts)
  
  #   Center-outward traversal
  center_i <- which.min(abs(spline.pts))
  left  <- seq(center_i - 1, 1,         by = -1)
  right <- seq(center_i + 1, grid.npts, by =  1)
  pairs <- as.vector(rbind(
    c(right, rep(NA, max(0, length(left)  - length(right)))),
    c(left,  rep(NA, max(0, length(right) - length(left))))
  ))
  pairs <- pairs[!is.na(pairs)]
  visit_order <- c(center_i, pairs)
  
  apl        <- matrix(0.0, nrow = ntags, ncol = grid.npts)
  beta_cache <- vector("list", grid.npts)
  
  for (idx in seq_along(visit_order)) {
    i     <- visit_order[idx]
    phi_i <- phi_common * 2^spline.pts[i]
    
    warm <- NULL
    if (i > 1        && !is.null(beta_cache[[i-1]])) warm <- beta_cache[[i-1]]
    if (i < grid.npts && !is.null(beta_cache[[i+1]])) {
      if (is.null(warm) ||
          abs(spline.pts[i+1] - spline.pts[i]) <
          abs(spline.pts[i-1] - spline.pts[i]))
        warm <- beta_cache[[i+1]]
    }
    
    out <- adjustedProfileLikBetaBin(phi_i, y = y, n = n,
                                     design = design,
                                     weights = weights,
                                     start = warm,
                                     get.coef = TRUE,
                                     maxit = maxit.grid,
                                     tol = tol.grid)
    apl[, i]        <- out$apl
    beta_cache[[i]] <- out$beta
  }
  
  d <- maximizeInterpolantBeta(spline.pts, apl)
  phi_common * 2^d
}
