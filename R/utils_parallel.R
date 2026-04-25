# utils_parallel.R
#
# Internal worker functions for parallel execution in ComBat_met.

.check_one_tag <- function(k, bv, design, batches_ind,
                           n_batch, mean.only.vec) {
  full_mat <- cbind(design, bv[k, ])
  nona     <- which(stats::complete.cases(full_mat))
  zv  <- FALSE
  zvb <- FALSE
  if (length(nona) == 0 ||
      qr(full_mat[nona, ])$rank < ncol(full_mat)) {
    zv <- TRUE
    return(c(zv = 1L, zvb = 0L))
  }
  if (!mean.only.vec[k]) {
    for (i in seq_along(batches_ind)) {
      idx <- intersect(batches_ind[[i]], nona)
      sub <- full_mat[idx, c(i, (n_batch + 1):ncol(full_mat)), drop = FALSE]
      if (qr(sub)$rank < ncol(full_mat) - n_batch + 1) {
        zvb <- TRUE
        break
      }
    }
  }
  c(zv = as.integer(zv), zvb = as.integer(zvb))
}

.est_phi_task <- function(task, bv, batches_ind, mod, n_batches,
                          batchmod, design, phi_common_per_batch,
                          verbose = FALSE) {
  if (verbose) {
    cat(paste0("Batch ", task$j, "\n"))
  }
  j          <- task$j
  chunk_tags <- task$chunk_tags
  idx_all    <- batches_ind[[j]]
  bv_b       <- bv[chunk_tags, idx_all, drop = FALSE]
  des_b      <- if (n_batches[j] > (ncol(design) - ncol(batchmod) + 1) &&
                    qr(mod[idx_all, , drop = FALSE])$rank == ncol(mod)) 
    mod[idx_all, , drop = FALSE] 
  else
    matrix(1, n_batches[j], 1)
  phi_t <- estimateGLMTagwiseDispBeta(bv_b, design = des_b,
                                      phi_common = phi_common_per_batch[j])
  list(j = j, chunk_tags = chunk_tags, phi_t = phi_t)
}

.adjust_chunk <- function(rows, bv, mu_hat_mat, phi_hat_mat,
                          delta_hat_mat, mu_star_mat, phi_star_mat,
                          batches_ind, n_batch) {
  res <- matrix(NA_real_, nrow = length(rows), ncol = ncol(bv))
  for (kk in seq_len(n_batch)) {
    cols    <- batches_ind[[kk]]
    bv_sub  <- bv[rows, cols, drop = FALSE]
    old_mu  <- mu_hat_mat[rows, cols, drop = FALSE]
    old_phi <- phi_hat_mat[rows, cols, drop = FALSE] *
      exp(delta_hat_mat[rows, kk])
    new_mu  <- mu_star_mat[rows, cols, drop = FALSE]
    new_phi <- phi_star_mat[rows, cols, drop = FALSE]
    res[, cols] <- match_quantiles_beta(bv_sub  = bv_sub,
                                        old_mu  = old_mu,
                                        old_phi = old_phi,
                                        new_mu  = new_mu,
                                        new_phi = new_phi)
  }
  res
}
