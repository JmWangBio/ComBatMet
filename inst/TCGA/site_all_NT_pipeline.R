
## Author: Junmin Wang
## Date: July 25th, 2024
## IMPORTANT NOTE: To run this script successfully, you need to have the TCGAbiolinks, dplyr, and stringr R packages installed.
## Make sure to change the paths of the outputs to where you want to save them (lines 51, 64, and 77).
## Also, please make sure to provide the correct path to "BRCA_met.RData", the output from download_TCGA_data.R (line 16).
## This script uses multiple methods to adjust the non-tumor site-level methylation data for batch effects.

## load libraries
library(TCGAbiolinks)
library(dplyr)
library(stringr)
library(ComBatMet)

## load data
load("path/to/BRCA_met.RData")

## extract assay data
BRCA.met.mat <- SummarizedExperiment::assay(data)

## select normal samples, i.e., "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(BRCA.met.mat),
  typesample = c("NT")
)

## get all batches containing more than one sample
samplesNT <- samplesNT[-grep(
  "A19Z",  # one-sample batch
  sapply(samplesNT, function(x) {
    str_split(x, "-")[[1]][6]
  }))]

## remove one-sample batches
samplesNT_subset_ordered <- get_IDs(data) %>% 
  filter(barcode %in% samplesNT) %>%
  arrange(plate) %>%
  pull(barcode)

BRCA.met.mat.subset <- BRCA.met.mat[, samplesNT_subset_ordered]

## remove rows with missing values
BRCA.met.mat.subset.complete <- BRCA.met.mat.subset[
  complete.cases(BRCA.met.mat.subset), 
]

###################
## No adjustment ##
###################
saveRDS(BRCA.met.mat.subset.complete, 
        file = "path/to/NT_raw_site_dat.rds")

############################################
## M-value ComBat (reference batch: A12R) ##
############################################
BRCA.met.mat.subset.complete.mcb <- Mvalue_ComBat(bv = BRCA.met.mat.subset.complete, 
                                                  batch = factor(
                                                    sapply(samplesNT_subset_ordered, 
                                                           function(x) {
                                                             str_split(x, "-")[[1]][6]
                                                           })),
                                                  ref.batch = "A12R")
saveRDS(BRCA.met.mat.subset.complete.mcb, 
        file = "path/to/NT_Mvalue_ComBat_site_dat.rds")

########################################
## ComBat-met (reference batch: A12R) ##
########################################
BRCA.met.mat.subset.complete.cbm <- ComBat_met(bv = BRCA.met.mat.subset.complete, 
                                               batch = factor(
                                                 sapply(samplesNT_subset_ordered, 
                                                        function(x) {
                                                          str_split(x, "-")[[1]][6]
                                                        })), 
                                               ref.batch = "A12R")
saveRDS(BRCA.met.mat.subset.complete.cbm, 
        file = "path/to/NT_ComBat_met_site_dat.rds")
