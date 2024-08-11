
## Author: Junmin Wang
## Date: July 25th, 2024
## IMPORTANT NOTE: To run this script successfully, you need to have the TCGAbiolinks, dplyr, tibble, tidyr, 
## stringr, IlluminaHumanMethylation450kprobe, edgeR, and reshape2 R packages installed.
## Make sure to change the paths of the outputs to where you want to save them (lines 82, 95, and 108).
## Also, please make sure to provide the correct path to "BRCA_met.RData", the output from download_TCGA_data.R (line 22).
## This script uses multiple methods to adjust the non-tumor gene-level methylation data for batch effects.

## load libraries
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(IlluminaHumanMethylation450kprobe)
library(edgeR)
library(reshape2)
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

## find nearest TSS for all HM450 probes
data("IlluminaHumanMethylation450kprobe")
HM450_TSS <- edgeR::nearestTSS(chr = IlluminaHumanMethylation450kprobe$chr, 
                               locus = IlluminaHumanMethylation450kprobe$site, 
                               species = "Hs")

## find HM450 sites only in the promoter regions
HM450_probe_TSS <- cbind(IlluminaHumanMethylation450kprobe[, c("Probe_ID", "chr", "site")], 
                         HM450_TSS)
HM450_probe_TSS_prreg <- HM450_probe_TSS[(HM450_probe_TSS$distance > 0) & 
                                           (HM450_probe_TSS$distance < 2000), ]

## keep data for sites only in the promoter regions; merge sites to genes
BRCA.met.mat.subset.complete.long <- reshape2::melt(BRCA.met.mat.subset.complete)
colnames(BRCA.met.mat.subset.complete.long) <- c("Probe_ID", "sample", "value")

BRCA.gene.mat.subset.complete <- BRCA.met.mat.subset.complete.long %>% 
  inner_join(HM450_probe_TSS_prreg, by = "Probe_ID") %>%
  group_by(sample, symbol) %>%
  summarise(value = mean(value)) %>%
  pivot_wider(names_from = "sample", 
              values_from = "value") %>%
  column_to_rownames(var = "symbol") %>%
  as.matrix()

###################
## No adjustment ##
###################
saveRDS(BRCA.gene.mat.subset.complete, 
        file = "path/to/NT_raw_gene_dat.rds")

############################################
## M-value ComBat (reference batch: A12R) ##
############################################
BRCA.gene.mat.subset.complete.mcb <- Mvalue_ComBat(bv = BRCA.gene.mat.subset.complete, 
                                                   batch = factor(
                                                     sapply(samplesNT_subset_ordered, 
                                                            function(x) {
                                                              str_split(x, "-")[[1]][6]
                                                            })),
                                                   ref.batch = "A12R")
saveRDS(BRCA.gene.mat.subset.complete.mcb, 
        file = "path/to/NT_Mvalue_ComBat_gene_dat.rds")

########################################
## ComBat-met (reference batch: A12R) ##
########################################
BRCA.gene.mat.subset.complete.cbm <- ComBat_met(bv = BRCA.gene.mat.subset.complete, 
                                                batch = factor(
                                                  sapply(samplesNT_subset_ordered, 
                                                         function(x) {
                                                           str_split(x, "-")[[1]][6]
                                                         })), 
                                                ref.batch = "A12R")
saveRDS(BRCA.gene.mat.subset.complete.cbm, 
        file = "path/to/NT_ComBat_met_gene_dat.rds")
