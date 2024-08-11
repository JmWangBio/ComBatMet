
## Author: Junmin Wang
## Date: July 18th, 2024
## IMPORTANT NOTE: To run this script successfully, you need to have the emdbook R package installed.
## This function is adapted from the dataSim() function in the methylKit R package.
## It can be used to simulate bisulfite sequencing data with two treatment conditions in two batches.

dataSim_custom <- function(replicates, sites, treatment, batch,
                           percentage = 10, effect = 10, batch.effect = 5,
                           alpha = 0.4, beta = 0.5,
                           theta = 10, disp.batch.effect = 1,
                           sample.ids = NULL, assembly = "hg18",
                           context = "CpG", add.info = FALSE)
{
  if (length(treatment) != replicates) {
    stop("treatment must be of same length as requested number of replicates")
  }
  if (replicates < 1) {
    stop("number of replicates has to be >= 1.")
  }
  if (replicates == 1) {
    warning("single replicate methylBase cannot be used downstream.")
  }
  if (!is.null(sample.ids)) {
    if (length(sample.ids) != replicates) {
      stop("sample.ids must be of same length as requested number of replicates")
    }
  }
  else {
    sample.ids <- ifelse(treatment == 1, paste0("test", cumsum(treatment)), 
                         paste0("ctrl", cumsum(!treatment)))
  }
  
  ## prepare parameters for simulation (x should be between batch.effect and 1 - batch.effect)
  raw <- matrix(ncol = replicates * 3, nrow = sites)
  index <- seq(1, replicates * 3, by = 3)
  x <- batch.effect / 100 + (1 - 2 * batch.effect / 100) * rbeta(sites, alpha, beta)
  
  ## check if baseline methylation percentage + effect + batch.effect < 1 for enough sites
  if (sum(x < 1 - effect / 100 - batch.effect / 100) >= floor(sites * percentage / 100)) {
    treatment_indices <- sample((1:sites)[x < 1 - effect / 100 - batch.effect / 100], 
                                size = floor(sites * percentage / 100))
  } else {
    stop("Error encountered during simulation. Either effect size or alpha and beta needs to be adjusted.")
  }
  
  ## split sites into two groups, based on the sign of batch effects
  batch_sites_split <- sample(c(1, 2), sites, replace = TRUE)
  batch_up_sites <- which(batch_sites_split == 1)
  batch_down_sites <- which(batch_sites_split == 2)
  
  ## loop through each replicate
  for (i in 1:replicates) {
    ## simulate coverage
    coverage <- rnbinom(sites, 1, 0.01)
    coverage <- ifelse(coverage < 10, 10, coverage)
    raw[, index[i]] <- coverage
    
    ## x: baseline methylation percentage
    ## y: baseline methylation percentage + treatment effect +/- batch effect
    ## theta_rep: dispersion parameter for the replicate
    y <- x
    theta_rep <- theta
    if (treatment[i] == 1) { # treatment: either 0 or 1
      effects <- effect / 100
      y[treatment_indices] <- y[treatment_indices] + effects
    }
    if (batch[i] == 1) { # batch: either 0 or 1
      batch.effects.up <- batch.effect / 100
      batch.effects.down <- batch.effect / 100
      y[batch_up_sites] <- y[batch_up_sites] + batch.effects.up
      y[batch_down_sites] <- y[batch_down_sites] - batch.effects.down
      theta_rep <- theta * disp.batch.effect
    }
    
    ## simulate methylated reads
    raw[, index[i] + 1] <- ceiling(coverage * emdbook::rbetabinom(
      n = sites, prob = y, size = 50, theta = theta_rep
    ) / 50)
    
    ## calculate unmethylated reads
    raw[, index[i] + 2] <- coverage - raw[, index[i] + 1]
  }
  
  ## simulate chr, start, end, strand info
  df <- raw
  df = as.data.frame(df)
  info <- data.frame(chr = rep("chr1", times = sites), start = 1:sites, 
                     end = 2:(sites + 1), strand = rep("+", times = sites))
  df <- cbind(info, df)
  
  ## create a methylBase object
  coverage.ind = seq(5, by = 3, length.out = replicates)
  numCs.ind = coverage.ind + 1
  numTs.ind = coverage.ind + 2
  names(df)[coverage.ind] = paste(c("coverage"), 1:replicates,
                                  sep = "")
  names(df)[numCs.ind] = paste(c("numCs"), 1:replicates, sep = "")
  names(df)[numTs.ind] = paste(c("numTs"), 1:replicates, sep = "")
  obj = new("methylBase", df, sample.ids = sample.ids, assembly = assembly, 
            context = context, treatment = treatment, coverage.index = coverage.ind, 
            numCs.index = numCs.ind, numTs.index = numTs.ind, destranded = FALSE, 
            resolution = "base")
  
  ## create output
  if (!add.info) {
    obj
  }
  else {
    addinfo = treatment_indices
    list(obj, addinfo, batch_up_sites, batch_down_sites)
  }
}
