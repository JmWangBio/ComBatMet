
# Custom wrapper function for RUVm
apply_RUVm_DE <- function(bvmat, group, k = 2, pseudo_beta = 1e-4) {
  # Check for valid inputs
  if (k <= 0) {
    stop("k must be a positive integer.")
  }
  
  if (length(group) != ncol(bvmat)) {
    stop("group must have the same length as the number of columns in bvmat.")
  }
  
  # Convert extreme 0 or 1 values to pseudo-beta
  if (pseudo_beta <= 0 | pseudo_beta >= 0.5) {
    stop("invalid pseudo beta-values.")
  }
  bvmat[bvmat == 0] <- pseudo_beta
  bvmat[bvmat == 1] <- 1 - pseudo_beta
  
  # Convert beta-values to M-values
  mvmat <- log2(bvmat / (1 - bvmat))
  
  # Keep rows with nonzero variances
  non_zero_var_idx <- rowVars(mvmat) > 0
  mvmat_nona <- mvmat[non_zero_var_idx, ]
  
  # Set up the design matrix
  design <- model.matrix(~group)
  
  # Fit the model using limma on M-values
  fit <- lmFit(mvmat_nona, design)
  
  # Apply empirical Bayes moderation
  fit <- eBayes(fit)

  # Print limma results
  limma_results <- topTable(fit, coef = 2, num = Inf)
  
  # Identify negative controls
  neg_control <- rownames(mvmat_nona) %in% rownames(limma_results[limma_results$adj.P.Val > 0.9, ])
  neg_control[is.na(neg_control)] <- FALSE
  
  # Perform RUV adjustment
  ruv_fit_nona <- RUVfit(Y = mvmat_nona, X = group, ctl = neg_control)
  
  # Apply the adjustment
  ruv_adj_fit_nona <- RUVadj(Y = mvmat_nona, fit = ruv_fit_nona)
  
  # Result after RUV correction
  ruv_adj_res <- data.frame(matrix(NA, nrow = 1000, ncol = 9))
  colnames(ruv_adj_res) <- c("F.p", "F.p.BH", "p_X1", "p.BH_X1", 
                             "b_X1", "sigma2", "var.b_X1", "fit.ctl", 
                             "mean")
  ruv_adj_res_noa <- topRUV(ruv_adj_fit_nona, num = Inf)
  ruv_adj_res[rownames(ruv_adj_res_noa), ] <- ruv_adj_res_noa

  # Reformat the result
  final_res <- data.frame(matrix(NA, nrow = 1000, ncol = 4))
  colnames(final_res) <- c("id", "t", "pval", "adj.pval")
  final_res$id <- 1:1000
  final_res$t <- ruv_adj_res$b_X1 / sqrt(ruv_adj_res$var.b_X1)
  final_res$pval <- ruv_adj_res$p_X1
  final_res$adj.pval <- ruv_adj_res$p.BH_X1
  final_res$t <- ifelse(is.na(final_res$pval), 0, final_res$t)
  final_res$pval <- ifelse(is.na(final_res$pval), 1, final_res$pval)
  final_res$adj.pval <- ifelse(is.na(final_res$adj.pval), 1, final_res$adj.pval)
  
  return(final_res)
}
