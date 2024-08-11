
## Author: Junmin Wang
## Date: July 18th, 2024
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, methylKit, methylSig, DSS, emdbook, limma, and sva R packages installed.
## Make sure to change the path of the output to where you want to save it (line 197).
## Also, please make sure to provide the correct paths to "dataSim_custom.R" and "diff_met_wf.R" (lines 21 and 22).
## This script can be used to simulate bisulfite sequencing data with varying mean and dispersion batch effects.
## The simulated data are adjusted for batch effects using different methods followed by diff. met. analysis.
## Beware that running 1000 simulations requires lots of computing power, so using a high performance computing environment is strongly recommended.
## To run fewer simulations, change "Nsims" from 1000 to 1 (line 29).

## load libraries
library(dplyr)
library(methylKit)
library(methylSig)
library(DSS)
library(limma)
library(ComBatMet)

## load simulation functions
source("path/to/dataSim_custom.R")
source("path/to/diff_met_wf.R")

## variables
batch_effect_lst <- c(0, 2, 5, 10)
disp_batch_effect_lst <- c(1, 2, 5, 10)

## fixed parameters
Nsims <- 1000
Nsites <- 1000
replicates <- 20
treatment <- c(rep(0, replicates / 2), rep(1, replicates / 2))
batch <- rep(rep(c(0, 1), each = replicates / 4), 2)

## placeholder for results
pval.df.all <- data.frame()

## loop through simulation indices
for (j in 1:Nsims) {
  pval.df <- data.frame()
  
  ## loop through possible values of dispersion batch effects
  for (disp_batch_effect in disp_batch_effect_lst) {
    ## loop through possible values of mean batch effects
    for (batch_effect in batch_effect_lst) {
      cat(sprintf("\ndisp_batch_effect: %s; batch_effect: %s\n\n", 
                  disp_batch_effect, batch_effect))
      
      ## simulate data
      my.methyldata <- dataSim_custom(replicates = replicates,
                                      sites = Nsites,
                                      treatment = treatment,
                                      batch = batch,
                                      batch.effect = batch_effect,
                                      disp.batch.effect = disp_batch_effect,
                                      add.info = TRUE)
      
      ## calculate beta values
      coverage.dat <- as.matrix(
        methylKit::getData(my.methyldata[[1]])[, paste0("coverage", 1:replicates)]
      )
      numCs.dat <- as.matrix(
        methylKit::getData(my.methyldata[[1]])[, paste0("numCs", 1:replicates)]
      )
      numCs.dat[is.na(numCs.dat)] <- 0
      beta.dat <- numCs.dat / coverage.dat
      colnames(beta.dat) <- paste0("beta", 1:replicates)
      
      ##########################
      #### batch adjustment ####
      ##########################
      # M-value ComBat
      cat("Running M-value ComBat...\n")
      beta.dat.mcb <- Mvalue_ComBat(beta.dat,
                                    batch = batch, 
                                    group = treatment, 
                                    full_mod = TRUE,
                                    pseudo_beta = 1e-4)
      
      # ComBat-met
      cat("Running Combat-met...\n")
      beta.dat.cbm <- ComBat_met(beta.dat,
                                 batch = batch,
                                 group = treatment,
                                 full_mod = TRUE,
                                 pseudo_beta = 1e-4,
                                 shrink = FALSE)
      
      # ComBat-biseq
      cat("Running ComBat-biseq...\n")
      numCs.dat.cbbs <- ComBat_biseq(coverage = coverage.dat, 
                                     numCs = numCs.dat, 
                                     batch = batch, 
                                     group = treatment, 
                                     full_mod = TRUE)
      
      ###########################################
      #### differential methylation analysis ####
      ###########################################
      ## no adjustment + t-test
      cat("Running t-test on results from no adjustment...\n")
      raw.de.res <- diff_met_wf(beta.dat, 
                                group = treatment, 
                                contrast.pair = "group1-group0",
                                pseudo_beta = 1e-4,
                                eb = FALSE)
      
      ## including batch as a covariate in t-test
      cat("Including batch as a covariate in t-test...\n")
      covar.de.res <- diff_met_wf(beta.dat, 
                                  group = treatment, 
                                  contrast.pair = "group1-group0",
                                  pseudo_beta = 1e-4,
                                  eb = FALSE,
                                  batch = batch)
      
      ## SVA + t-test
      cat("Including SVA surrogate variables as covariates in t-test...\n")
      sva.de.res <- diff_met_wf(beta.dat, 
                                group = treatment, 
                                contrast.pair = "group1-group0",
                                pseudo_beta = 1e-4,
                                eb = FALSE,
                                batch = NULL,
                                sva = TRUE)
      
      ## M-value ComBat + t-test
      cat("Running t-test on results from M-value ComBat...\n")
      mcb.de.res <- diff_met_wf(beta.dat.mcb, 
                                group = treatment, 
                                contrast.pair = "group1-group0",
                                pseudo_beta = 1e-4,
                                eb = FALSE)
      
      ## ComBat-met + t-test
      cat("Running t-test on results from ComBat-met...\n")
      cbm.de.res <- diff_met_wf(beta.dat.cbm,
                                group = treatment,
                                contrast.pair = "group1-group0",
                                pseudo_beta = 1e-4,
                                eb = FALSE)
      
      ## ComBat-biseq + likelihood ratio test
      cat("Running MethylSig on results from ComBat-biseq...\n")
      adj.BS.dat.lst <- lapply(1:ncol(numCs.dat.cbbs), function(i) {
        data.frame(chr = methylKit::getData(my.methyldata[[1]])$chr,
                   pos = methylKit::getData(my.methyldata[[1]])$start,
                   N = coverage.dat[, i],
                   X = numCs.dat.cbbs[, i])
      })
      
      adj.BSseq <- makeBSseqData(adj.BS.dat.lst, 1:length(adj.BS.dat.lst))
      pData(adj.BSseq) <- data.frame(Type = treatment,
                                     Batch = batch)
      
      diff_gr <- diff_methylsig(
        bs = adj.BSseq,
        group_column = 'Type',
        comparison_groups = c('case' = '1', 'control' = '0'),
        disp_groups = c('case' = TRUE, 'control' = TRUE),
        local_window_size = 0,
        t_approx = TRUE,
        n_cores = 1)
      
      ## adjust p-values; label ground truth (positive or negative)
      tmp.pval.df <- rbind(raw.de.res[, c("id", "pval", "adj.pval")] %>% 
                             mutate(method = "raw + t-test"),
                           covar.de.res[, c("id", "pval", "adj.pval")] %>% 
                             mutate(method = "covariate in t-test"),
                           sva.de.res[, c("id", "pval", "adj.pval")] %>%
                             mutate(method = "SVA + t-test"),
                           mcb.de.res[, c("id", "pval", "adj.pval")] %>% 
                             mutate(method = "M-value ComBat + t-test"),
                           cbm.de.res[, c("id", "pval", "adj.pval")] %>% 
                             mutate(method = "ComBat-Met + t-test"),
                           data.frame(id = 1:nrow(beta.dat), pval = diff_gr$pvalue) %>% 
                             mutate(adj.pval = p.adjust(pval, method = "BH"),
                                    method = "ComBat-biseq + lrt")) %>%
        inner_join(data.frame(id = 1:nrow(beta.dat),
                              TP = ifelse(1:nrow(beta.dat) %in% my.methyldata[[2]], 
                                          TRUE, FALSE)),
                   by = "id")
      
      pval.df <- rbind(pval.df, tmp.pval.df %>%
                         mutate(disp_batch_effect = disp_batch_effect,
                                batch_effect = batch_effect))
    }
  }
  
  # add simulation index
  pval.df$sim <- j
  pval.df.all <- rbind(pval.df.all, pval.df)
}

# save results
saveRDS(pval.df.all, 
        file = 'path/to/sim_all_DE_pval_data.rds')
