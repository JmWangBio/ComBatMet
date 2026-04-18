#  dispCoxReidBeta.R
#
#  Estimates a single common precision phi for beta regression by
#  maximising the sum of APL values across all rows of y.

dispCoxReidBeta <- function(y, design=NULL, weights=NULL,
                            interval=c(1e-4, 1e4),
                            tol=1e-5, subset=10000)
  # y        matrix of proportions (ntags x nsamples), each in (0,1)
  # design   design matrix
  # weights  optional weight matrix of dim(y)
  # interval search interval for phi
  # tol      optimization tolerance passed to optimize()
  # subset   maximum number of rows to use for the subset
  #
  # Returns: scalar common precision estimate
{
  y <- as.matrix(y)
  design <- as.matrix(design)

  # Systematic subsetting when y is large
  if (!is.null(subset) && nrow(y) > subset) {
    avg_logit <- rowMeans(log(y / (1 - y)))
    i <- systematicSubset(subset, avg_logit)
    y <- y[i, , drop=FALSE]
    weights <- if (!is.null(weights)) weights[i, , drop=FALSE] else NULL
  }

  # Objective function: sum of APL across all rows at a given phi
  fun <- function(par, y, design, weights) {
    phi_val <- par^4
    sum(adjustedProfileLikBeta(phi_val, y=y, design=design,
                               weights=weights))
  }

  out <- stats::optimize(f = fun,
                         interval = interval^0.25,
                         y = y,
                         design = design,
                         weights = weights,
                         maximum = TRUE,
                         tol = tol)

  out$maximum^4
}
