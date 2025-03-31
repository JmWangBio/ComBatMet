
## Author: Junmin Wang
## Date: July 18th, 2024
## This script can be used for differential methylation analysis via minfi.

diff_met_wf_minfi <- function(bv, group, pseudo_beta = 1e-4) {
  ## convert extreme 0 or 1 values to pseudo-beta
  if (pseudo_beta <= 0 | pseudo_beta >= 0.5) {
    stop("Invalid pseudo beta-values.")
  }
  bv[bv <= 0] <- pseudo_beta
  bv[bv >= 1] <- 1 - pseudo_beta
  
  ## convert beta values to M values
  mv <- log(bv / (1 - bv))
  
  ## dmpFinder
  dmp <- minfi::dmpFinder(mv, pheno=group, type="categorical")
  
  t_stat <- dmp$f
  pval <- dmp$pval
  adj.pval <- dmp$qval
  
  ## create output
  res <- data.frame(id = as.numeric(rownames(dmp)),
                    t = t_stat,
                    pval = pval,
                    adj.pval = adj.pval) %>%
    arrange(id)
  
  ## correct the p-values for zero-variance rows
  res[as.numeric(which(rowVars(bv) == 0)), c("pval", "adj.pval")] <- 1
  
  return(res)
}
