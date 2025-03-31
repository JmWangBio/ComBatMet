
####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times) {
  return(matrix(rep(vec, n_times), ncol = n_times, byrow = FALSE))
}

#### Monte Carlo integration functions for beta distribution
monte_carlo_int_beta <- function(dat, mu, gamma, phi, delta, feature.subset.n, ncores) {
  cl <- parallel::makeCluster(ncores)
  
  MC_int <- function(i) {
    pos_res <- c(gamma.star = NA, 
                 delta.star = NA)
    
    m <- mu[-i, !is.na(dat[i, ])]
    p <- phi[-i, !is.na(dat[i, ])]
    x <- dat[i, !is.na(dat[i, ])]
    gamma_sub <- gamma[-i]
    delta_sub <- delta[-i]
    
    # take a subset of features to do integration - save time
    if (!is.null(feature.subset.n) & is.numeric(feature.subset.n) & length(feature.subset.n) == 1) {
      mcint_ind <- sample(1:(nrow(dat) - 1), feature.subset.n, replace = FALSE)
      m <- m[mcint_ind, ]
      p <- p[mcint_ind, ]
      gamma_sub <- gamma_sub[mcint_ind]
      delta_sub <- delta_sub[mcint_ind]
      G_sub <- feature.subset.n
    } else {
      G_sub <- nrow(dat) - 1
    }
    
    LH <- sapply(1:G_sub, function(j) {
      prod(stats::dbeta(x, shape1 = m[j, ] * (p[j, ] * exp(delta_sub[j])), 
                        shape2 = (1 - m[j, ]) * (p[j, ] * exp(delta_sub[j]))))
    })
    LH[is.na(LH)] <- 0
    if (sum(LH) == 0 | is.na(sum(LH))){
      pos_res <- c(gamma.star = as.numeric(gamma[i]),
                   delta.star = as.numeric(delta[i]))
    } else {
      pos_res["gamma.star"] <- sum(gamma_sub * LH, na.rm = TRUE) / 
        sum(LH[!is.na(gamma_sub)])
      pos_res["delta.star"] <- sum(delta_sub * LH, na.rm = TRUE) / 
        sum(LH[!is.na(delta_sub)])
    }
    return(pos_res)
  }
  
  # run in parallel
  pos_res_lst <- parallel::parLapply(cl, 1:nrow(dat), MC_int)
  parallel::stopCluster(cl)
  
  pos_res_mat <- do.call(rbind, pos_res_lst)
  res <- list(gamma_star = pos_res_mat[, "gamma.star"], 
              delta_star = pos_res_mat[, "delta.star"]) 
  return(res)
} 

#### Monte Carlo integration functions for beta-binomial distributions
monte_carlo_int_betabin <- function(numCs.dat, coverage.dat, mu, gamma, 
                                    phi, delta, feature.subset.n, ncores) {
  cl <- parallel::makeCluster(ncores)
  
  MC_int <- function(i) {
    pos_res <- c(gamma.star = NA, 
                 delta.star = NA)
    
    m <- mu[-i, !is.na(numCs.dat[i, ]) & !is.na(coverage.dat[i, ])]
    p <- phi[-i, !is.na(numCs.dat[i, ]) & !is.na(coverage.dat[i, ])]
    x <- numCs.dat[i, !is.na(numCs.dat[i, ]) & !is.na(coverage.dat[i, ])]
    y <- coverage.dat[i, !is.na(numCs.dat[i, ]) & !is.na(coverage.dat[i, ])]
    gamma_sub <- gamma[-i]
    delta_sub <- delta[-i]
    
    # take a subset of features to do integration - save time
    if (!is.null(feature.subset.n) & is.numeric(feature.subset.n) & length(feature.subset.n) == 1) {
      mcint_ind <- sample(1:(nrow(numCs.dat) - 1), feature.subset.n, replace = FALSE)
      m <- m[mcint_ind, ]
      p <- p[mcint_ind, ]
      gamma_sub <- gamma_sub[mcint_ind]
      delta_sub <- delta_sub[mcint_ind]
      G_sub <- feature.subset.n
    } else {
      G_sub <- nrow(numCs.dat)-1
    }
    
    LH <- sapply(1:G_sub, function(j) {
      prod(updog::dbetabinom(x, 
                             size = y, 
                             mu = m[j, ], 
                             rho = 1 / (p[j, ] * exp(delta_sub[j]) + 1),
                             log = FALSE))
    })
    LH[is.na(LH)] <- 0
    if (sum(LH) == 0 | is.na(sum(LH))){
      pos_res <- c(gamma.star = as.numeric(gamma[i]),
                   delta.star = as.numeric(delta[i]))
    } else {
      pos_res["gamma.star"] <- sum(gamma_sub * LH, na.rm = TRUE) / 
        sum(LH[!is.na(gamma_sub)])
      pos_res["delta.star"] <- sum(delta_sub * LH, na.rm = TRUE) / 
        sum(LH[!is.na(delta_sub)])
    }
    return(pos_res)
  }
  
  # run in parallel
  pos_res_lst <- parallel::parLapply(cl, 1:nrow(numCs.dat), MC_int)
  parallel::stopCluster(cl)
  
  pos_res_mat <- do.call(rbind, pos_res_lst)
  res <- list(gamma_star = pos_res_mat[, "gamma.star"], 
              delta_star = pos_res_mat[, "delta.star"])	
  return(res)
} 

#### Match quantiles for beta distribution
match_quantiles_beta <- function(bv_sub, old_mu, old_phi, new_mu, new_phi) {
  new_bv_sub <- matrix(NA, nrow = nrow(bv_sub), ncol = ncol(bv_sub))
  for (a in 1:nrow(bv_sub)) {
    for (b in 1:ncol(bv_sub)) {
      tmp_p <- stats::pbeta(bv_sub[a, b], shape1 = old_mu[a, b] * old_phi[a, b], 
                            shape2 = (1 - old_mu[a, b]) * old_phi[a, b])
      if (is.na(tmp_p)) {
        new_bv_sub[a, b] <- bv_sub[a, b]  
      } else {
        new_bv_sub[a, b] <- stats::qbeta(tmp_p, shape1 = new_mu[a, b] * new_phi[a, b], 
                                         shape2 = (1 - new_mu[a, b]) * new_phi[a, b])
      }
    }
  }
  return(new_bv_sub)
}

#### Match quantiles for beta-binomial distribution
match_quantiles_betabin <- function(numCs_sub, coverage_sub, old_mu, 
                                    old_phi, new_mu, new_phi) {
  new_numCs_sub <- matrix(NA, nrow = nrow(numCs_sub), ncol = ncol(numCs_sub))
  for (a in 1:nrow(numCs_sub)) {
    for (b in 1:ncol(numCs_sub)) {
      tmp_p <- updog::pbetabinom(numCs_sub[a, b],
                                 size = coverage_sub[a, b],
                                 mu = old_mu[a, b],
                                 rho = 1 / (1 + old_phi[a, b]),
                                 log_p = FALSE)
      if (is.na(tmp_p) | abs(tmp_p - 1) < 1e-4) {
        new_numCs_sub[a, b] <- numCs_sub[a, b]  
        # for outlier count, if p==1, will return incorrect values -> use original count instead
      } else {
        new_numCs_sub[a, b] <- updog::qbetabinom(tmp_p,
                                                 size = coverage_sub[a, b],
                                                 mu = new_mu[a, b],
                                                 rho = 1 / (1 + new_phi[a, b]))
      }
    }
  }
  return(new_numCs_sub)
}
