
## Author: Junmin Wang
## Date: August 6th, 2024
## This function calculates and returns a list of explained variations by batch and condition combinations.
## The code is adapted from the batchqc_explained_variation() function in the BatchQC package.
## bv: a matrix of beta values
## condition: condition covariate of interest
## batch: batch covariate

calc_explained_variation <- function(bv, condition, batch) {
  # convert beta values to M values
  mv <- log(bv / (1 - bv))
  
  # create design matrices
  if (nlevels(as.factor(condition)) > 1) {
    cond_mod <- model.matrix(~as.factor(condition))
  } else {
    cond_mod <- matrix(rep(1, ncol(bv)), ncol = 1)
  }
  batch_mod <- model.matrix(~as.factor(batch))
  mod <- cbind(cond_mod, batch_mod[, -1])
  mod00 <- matrix(rep(1, ncol(bv)), ncol = 1)
  
  # calculate residuals
  resid_00 <- mv - mv %*% mod00 %*% solve(t(mod00) %*% mod00) %*% t(mod00)
  rss_00 <- rowSums(resid_00 * resid_00)
  
  resid_full <- mv - mv %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  rss_full <- rowSums(resid_full * resid_full)
  
  resid_batch_only <- mv - mv %*% batch_mod %*% 
    solve(t(batch_mod) %*% batch_mod) %*% t(batch_mod)
  rss_batch_only <- rowSums(resid_batch_only * resid_batch_only)
  
  resid_cond_only <- mv - mv %*% cond_mod %*% 
    solve(t(cond_mod) %*% cond_mod) %*% t(cond_mod)
  rss_cond_only <- rowSums(resid_cond_only * resid_cond_only)
  
  # calculate R^2
  r2_full <- 1 - rss_full / rss_00
  r2_cond_only <- 1 - rss_cond_only / rss_00
  r2_batch_only <- 1 - rss_batch_only / rss_00
 
  # store result
  res <- list(ConditionAndBatch = r2_full,
              Condition = r2_cond_only,
              Batch = r2_batch_only)
  
  return(res)
}
