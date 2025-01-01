
## Author: Junmin Wang
## Date: December 31st, 2024

# load libraries
library(TCGAbiolinks)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(sva)

## load data
load("/path/to/BRCA_met.RData")

# extract assay data
BRCA.met.mat <- SummarizedExperiment::assay(data)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(
  barcode = colnames(BRCA.met.mat),
  typesample = c("NT")
)

samplesNT <- samplesNT[-grep("A19Z",  # one-sample batches
                             sapply(samplesNT, function(x) {
                               str_split(x, "-")[[1]][6]
                             }))]

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(
  barcode = colnames(BRCA.met.mat), 
  typesample = c("TP")
)

samplesTP <- data$barcode[data$paper_BRCA_Subtype_PAM50 %in% c("LumB")]
samplesTP <- samplesTP[-grep("A17F|A18O|A244|A29O|A32T|A33F|A36K",  # one-sample batches
                             sapply(samplesTP, function(x) {
                               str_split(x, "-")[[1]][6]
                             }))]

# remove one-sample batches
samplesNT_subset_ordered <- get_IDs(data) %>% 
  filter(barcode %in% samplesNT) %>%
  arrange(plate) %>%
  pull(barcode)

samplesTP_subset_ordered <- get_IDs(data) %>% 
  filter(barcode %in% samplesTP) %>%
  arrange(plate) %>%
  pull(barcode)

# make metadata data frame
metadata <- data.frame(barcode = c(samplesNT_subset_ordered, samplesTP_subset_ordered),
                       group = c(rep("NT", length(samplesNT_subset_ordered)),
                                 rep("TP", length(samplesTP_subset_ordered))))

# save metadata
write.csv(metadata, 
          file = "/path/to/metadata.csv", 
          row.names = FALSE)
