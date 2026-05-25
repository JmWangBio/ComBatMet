# utils_parallel.R
#
# Internal worker functions for parallel execution in ComBat_met
# and ComBat_biseq.

.check_one_tag <- function(inp, batches_ind, n_batch) {
  bv_chunk <- inp$bv_chunk
  design <- inp$design
  mean_only_vec <- inp$mean_only_vec
  n_tags <- nrow(bv_chunk)
  
  zv  <- logical(n_tags)
  zvb <- logical(n_tags)
  
  for (ki in seq_len(n_tags)) {
    full_mat <- cbind(design, bv_chunk[ki, ])
    nona     <- which(stats::complete.cases(full_mat))
    
    if (length(nona) == 0 ||
        qr(full_mat[nona, ])$rank < ncol(full_mat)) {
      zv[ki] <- TRUE
      next
    }
    if (!mean_only_vec[ki]) {
      for (i in seq_along(batches_ind)) {
        idx <- intersect(batches_ind[[i]], nona)
        sub <- full_mat[idx, c(i, (n_batch + 1):ncol(full_mat)), drop = FALSE]
        if (qr(sub)$rank < ncol(full_mat) - n_batch + 1) {
          zvb[ki] <- TRUE
          break
        }
      }
    }
  }
  list(zv = zv, zvb = zvb)
}

.est_phi_task <- function(inp, verbose = FALSE) {
  if (verbose) {
    cat(paste0("Batch ", inp$j, "\n"))
  }
  phi_t <- estimateGLMTagwiseDispBeta(inp$bv_b, 
                                      design = inp$des_b,
                                      phi_common = inp$phi_c)
  list(j = inp$j, chunk_tags = inp$chunk_tags, phi_t = phi_t)
}

.adjust_chunk <- function(inp, batches_ind, n_batch) {
  res <- matrix(NA_real_, nrow = nrow(inp$sub_bv), ncol = ncol(inp$sub_bv))
  for (kk in seq_len(n_batch)) {
    cols    <- batches_ind[[kk]]
    bv_sub  <- inp$sub_bv[, cols, drop = FALSE]
    old_mu  <- inp$sub_mu_hat_mat[, cols, drop = FALSE]
    old_phi <- inp$sub_phi_hat_mat[, cols, drop = FALSE] *
      exp(inp$sub_delta_hat_mat[, kk])
    new_mu  <- inp$sub_mu_star_mat[, cols, drop = FALSE]
    new_phi <- inp$sub_phi_star_mat[, cols, drop = FALSE]
    res[, cols] <- match_quantiles_beta(bv_sub  = bv_sub,
                                        old_mu  = old_mu,
                                        old_phi = old_phi,
                                        new_mu  = new_mu,
                                        new_phi = new_phi)
  }
  res
}

# ------------------------------------------------------------------ #
# Beta-binomial worker functions (ComBat_biseq)                      #
# ------------------------------------------------------------------ #

.check_one_tag_bb <- function(inp, batches_ind, n_batch) {
  numCs_chunk    <- inp$numCs_chunk
  cov_chunk      <- inp$cov_chunk
  design         <- inp$design
  mean_only_vec  <- inp$mean_only_vec
  n_tags         <- nrow(numCs_chunk)
  
  zv  <- logical(n_tags)
  
  for (ki in seq_len(n_tags)) {
    # zero-coverage samples treated as NA
    cov_row   <- cov_chunk[ki, ]
    numCs_row <- numCs_chunk[ki, ]
    cov_row[cov_row == 0]     <- NA
    numCs_row[is.na(cov_row)] <- NA
    
    full_mat <- as.data.frame(cbind(design,
                                    coverage = cov_row,
                                    numCs    = numCs_row))
    nona <- which(stats::complete.cases(full_mat))
    
    if (length(nona) == 0) {
      zv[ki] <- TRUE
      next
    }
    
    # check model variance using empirical proportion
    full_ratio <- full_mat$numCs / full_mat$coverage
    if (qr(cbind(full_mat[, 1:(ncol(full_mat)-2)],
                 full_ratio)[nona, ])$rank < ncol(full_mat) - 1) {
      zv[ki] <- TRUE
    }
  }
  list(zv = zv)
}

.est_phi_task_bb <- function(inp, verbose = FALSE) {
  if (verbose) cat(paste0("Batch ", inp$j, "\n"))
  phi_t <- estimateGLMTagwiseDispBetaBin(inp$numCs_b,
                                         n          = inp$cov_b,
                                         design     = inp$des_b,
                                         phi_common = inp$phi_c)
  list(j = inp$j, chunk_tags = inp$chunk_tags, phi_t = phi_t)
}

.adjust_chunk_bb <- function(inp, batches_ind, n_batch) {
  res <- matrix(NA_real_,
                nrow = nrow(inp$sub_numCs),
                ncol = ncol(inp$sub_numCs))
  for (kk in seq_len(n_batch)) {
    cols         <- batches_ind[[kk]]
    numCs_sub    <- inp$sub_numCs[,    cols, drop = FALSE]
    coverage_sub <- inp$sub_coverage[, cols, drop = FALSE]
    old_mu       <- inp$sub_mu_hat[,   cols, drop = FALSE]
    old_phi      <- inp$sub_phi_hat[,  cols, drop = FALSE] *
      exp(inp$sub_delta_hat[, kk])
    new_mu       <- inp$sub_mu_star[,  cols, drop = FALSE]
    new_phi      <- inp$sub_phi_star[, cols, drop = FALSE]
    res[, cols]  <- match_quantiles_betabin(
      numCs_sub    = numCs_sub,
      coverage_sub = coverage_sub,
      old_mu       = old_mu,
      old_phi      = old_phi,
      new_mu       = new_mu,
      new_phi      = new_phi)
  }
  res
}
