
## Author: Junmin Wang
## Date: January 8th, 2025
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, methylKit, DSS, minfi, BEclear, and missMethyl R packages installed.
## Make sure to change the path of the output to where you want to save it.
## Also, please make sure to provide the correct paths to "dataSim_custom.R", "diff_met_wf_minfi.R", and "apply_BEclear.R"
## This script can be used to simulate bisulfite sequencing data with varying mean and dispersion (precision) batch effects.
## The simulated data are adjusted for batch effects using different methods followed by diff. met. analysis.
## Beware that running 1000 simulations requires lots of computing power, so using a high performance computing environment is strongly recommended.
## To run fewer simulations, change "Nsims" from 1000 to 1.

## load libraries
library(dplyr)
library(methylKit)
library(DSS)
library(minfi)
library(ComBatMet)
library(BEclear)
library(missMethyl)

## load simulation functions
source("path/to/inst/simulation/dataSim_custom.R")
source("path/to/inst/simulation/diff_met_wf_minfi.R")
source("path/to/inst/simulation/helper/apply_BEclear.R")

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
      
      # BEclear
      cat("Running BEclear...\n")
      beta.dat.beclear <- apply_BEclear(beta.dat,
                                        batch = batch,
                                        group = treatment)
      
      # ComBat-met with shrink
      cat("Running Combat-met with shrink...\n")
      beta.dat.cbm.shrink <- ComBat_met(beta.dat,
                                        batch = batch,
                                        group = treatment,
                                        full_mod = TRUE,
                                        pseudo_beta = 1e-4,
                                        shrink = TRUE)
      
      ###########################################
      #### differential methylation analysis ####
      ###########################################
      ## no adjustment + t-test
      cat("Running t-test on results from no adjustment...\n")
      raw.de.res <- diff_met_wf_minfi(beta.dat, 
                                      group = treatment,
                                      pseudo_beta = 1e-4)
      
      ## M-value ComBat + t-test
      cat("Running t-test on results from M-value ComBat...\n")
      mcb.de.res <- diff_met_wf_minfi(beta.dat.mcb, 
                                      group = treatment,
                                      pseudo_beta = 1e-4)
      
      ## ComBat-met + t-test
      cat("Running t-test on results from ComBat-met...\n")
      cbm.de.res <- diff_met_wf_minfi(beta.dat.cbm,
                                      group = treatment,
                                      pseudo_beta = 1e-4)
      
      ## BEclear + t-test
      cat("Running t-test on results from BEclear...\n")
      beclear.de.res <- diff_met_wf_minfi(beta.dat.beclear,
                                          group = treatment,
                                          pseudo_beta = 1e-4)
      
      ## ComBat-met with shrink + t-test
      cat("Running t-test on results from ComBat-met with shrink...\n")
      cbm.shrink.de.res <- diff_met_wf_minfi(beta.dat.cbm.shrink,
                                             group = treatment,
                                             pseudo_beta = 1e-4)
      
      ## adjust p-values; label ground truth (positive or negative)
      tmp.pval.df <- rbind(raw.de.res[, c("id", "pval", "adj.pval")] %>% 
                             mutate(method = "raw + t-test"),
                           mcb.de.res[, c("id", "pval", "adj.pval")] %>% 
                             mutate(method = "M-value ComBat + t-test"),
                           cbm.de.res[, c("id", "pval", "adj.pval")] %>% 
                             mutate(method = "ComBat-Met + t-test"),
                           beclear.de.res[, c("id", "pval", "adj.pval")] %>%
                             mutate(method = "BEclear + t-test"),
                           cbm.shrink.de.res[, c("id", "pval", "adj.pval")] %>%
                             mutate(method = "ComBat-Met (shrink) + t-test")) %>%
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
        file = 'path/to/sim_all_DE_pval_data_minfi.rds')
