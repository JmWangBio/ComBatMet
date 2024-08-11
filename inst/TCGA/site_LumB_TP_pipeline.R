
## Author: Junmin Wang
## Date: July 25th, 2024
## IMPORTANT NOTE: To run this script successfully, you need to have the TCGAbiolinks, dplyr, and stringr R packages installed.
## Make sure to change the paths of the outputs to where you want to save them (lines 48, 61, and 74).
## Also, please make sure to provide the correct path to "BRCA_met.RData", the output from download_TCGA_data.R (line 16).
## This script uses multiple methods to adjust the luminal-B breast cancer site-level methylation site data for batch effects.

## load libraries
library(TCGAbiolinks)
library(dplyr)
library(stringr)
library(ComBatMet)

## load data
load("path/to/BRCA_met.RData")

## extract assay data
BRCA.met.mat <- SummarizedExperiment::assay(data)

## select luminal-B tumor samples
samplesTP <- data$barcode[data$paper_BRCA_Subtype_PAM50 %in% c("LumB")]

## get all batches containing more than one sample
samplesTP <- samplesTP[-grep(
  "A17F|A18O|A244|A29O|A32T|A33F|A36K",  # one-sample batches
  sapply(samplesTP, function(x) {
    str_split(x, "-")[[1]][6]
  }))]

## remove one-sample batches
samplesTP_subset_ordered <- get_IDs(data) %>% 
  filter(barcode %in% samplesTP) %>%
  arrange(plate) %>%
  pull(barcode)

BRCA.met.mat.subset <- BRCA.met.mat[, samplesTP_subset_ordered]

## remove rows with missing values
BRCA.met.mat.subset.complete <- BRCA.met.mat.subset[
  complete.cases(BRCA.met.mat.subset), 
]

###################
## No adjustment ##
###################
saveRDS(BRCA.met.mat.subset.complete, 
        file = "path/to/TP_raw_site_dat.rds")

############################################
## M-value ComBat (reference batch: A12R) ##
############################################
BRCA.met.mat.subset.complete.mcb <- Mvalue_ComBat(bv = BRCA.met.mat.subset.complete, 
                                                  batch = factor(
                                                    sapply(samplesTP_subset_ordered, 
                                                           function(x) {
                                                             str_split(x, "-")[[1]][6]
                                                           })),
                                                  ref.batch = "A12R")
saveRDS(BRCA.met.mat.subset.complete.mcb, 
        file = "path/to/TP_Mvalue_ComBat_site_dat.rds")

########################################
## ComBat-met (reference batch: A12R) ##
########################################
BRCA.met.mat.subset.complete.cbm <- ComBat_met(bv = BRCA.met.mat.subset.complete, 
                                               batch = factor(
                                                 sapply(samplesTP_subset_ordered, 
                                                        function(x) {
                                                          str_split(x, "-")[[1]][6]
                                                        })), 
                                               ref.batch = "A12R")
saveRDS(BRCA.met.mat.subset.complete.cbm, 
        file = "path/to/TP_ComBat_met_site_dat.rds")
