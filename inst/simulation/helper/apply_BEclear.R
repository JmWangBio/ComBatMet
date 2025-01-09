
# Custom wrapper function for BEclear
apply_BEclear <- function(bvmat, batch, group) {
  # Validate inputs
  if (length(batch) != ncol(bvmat) ||
      length(group) != ncol(bvmat)) {
    stop("The length of batch and group must match the number of samples in bvmat.")
  }
  
  # Get unique groups
  unique_grps <- unique(group)
  
  # Initialize a matrix to store the corrected data
  corrected_bvmat <- matrix(NA, nrow = nrow(bvmat), ncol = ncol(bvmat))
  colnames(corrected_bvmat) <- colnames(bvmat)
  rownames(corrected_bvmat) <- rownames(bvmat)
  
  # Loop through each group and apply BEclear
  for (grp in unique_grps) {
    # Subset data for the current group
    grp_samples <- which(group == grp)
    bvmat_grp <- bvmat[, grp_samples, drop = FALSE]
    batch_grp <- batch[grp_samples]
    
    # Apply BEclear to the subset
    res_grp <- correctBatchEffect(data = bvmat_grp, 
                                  samples = data.frame(sample_id = paste0("beta", grp_samples),
                                                       batch_id = batch_grp))
    corrected_bvmat_grp <- res_grp$correctedPredictedData
    
    # Store corrected back back into the final matrix
    corrected_bvmat[, grp_samples] <- corrected_bvmat_grp
  }
  
  return(corrected_bvmat)
}
