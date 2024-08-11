
## Author: Junmin Wang
## Date: August 6th, 2024
## This script calculates and plots % variation explained by batch in non-tumor and tumor data.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, stringr, and ggplot2
## packages installed. Also, make sure to provide the correct paths for "calc_explained_variation.R" (line 18) and "meta_data.rds" (line 30).
## "meta_data.rds" specifies the batch information and which samples are tumor samples. Ensure that you also provide the correct paths for "NT_raw_gene_dat.rds", 
## "NT_Mvalue_ComBat_gene_dat.rds", "NT_ComBat_met_gene_dat.rds", "TP_raw_gene_dat.rds", "TP_Mvalue_ComBat_gene_dat.rds",
## and "TP_ComBat_met_gene_dat.rds" (lines 21 - 27), which are the intermediate output files produced from 
## "gene_all_NT_pipeline.R" and "gene_LumB_TP_pipeline.R".

## load libraries
library(dplyr)
library(stringr)
library(ggplot2)

## load functions
source("path/to/calc_explained_variation.R")

## load beta-value data
NT_raw_gene_dat <- readRDS("path/to/NT_raw_gene_dat.rds")
NT_Mvalue_ComBat_gene_dat <- readRDS("path/to/NT_Mvalue_ComBat_gene_dat.rds")
NT_ComBat_met_gene_dat <- readRDS("path/to/NT_ComBat_met_gene_dat.rds")

TP_raw_gene_dat <- readRDS("path/to/TP_raw_gene_dat.rds")
TP_Mvalue_ComBat_gene_dat <- readRDS("path/to/TP_Mvalue_ComBat_gene_dat.rds")
TP_ComBat_met_gene_dat <- readRDS("path/to/TP_ComBat_met_gene_dat.rds")

## load meta-data
meta_data <- readRDS("path/to/meta_data.rds")
samplesNT <- meta_data %>% 
  filter(type == "Solid Tissue Normal") %>% 
  pull(ID)
samplesNT <- samplesNT[!grepl("A19Z", samplesNT)]  # A19Z is a one-sample batch.

samplesTP <- meta_data %>% 
  filter(type == "Primary Tumor",
         subtype == "LumB") %>%
  pull(ID)
samplesTP <- samplesTP[!grepl("A17F|A18O|A244|A29O|A32T|A33F|A36K", samplesTP)]  # A17F, A18O, A244, A29O, A32T, A33F, and A36K are one-sample batches.

# order sample by plate
samplesNT_subset_ordered <- meta_data %>% 
  filter(ID %in% samplesNT) %>%
  arrange(plate) %>%
  pull(ID)

samplesTP_subset_ordered <- meta_data %>%
  filter(ID %in% samplesTP) %>%
  arrange(plate) %>%
  pull(ID)

##########################
####### non-tumor ########
##########################
# retrieve batch and condition
NT_batch <- as.factor(sapply(samplesNT_subset_ordered, function(x) {
  str_split(x, "-")[[1]][6]
}))
NT_condition <- as.factor(c(rep(1, length(samplesNT_subset_ordered))))

# calculate % variance explained
NT_raw_gene_pve <- calc_explained_variation(NT_raw_gene_dat,
                                            condition = NT_condition,
                                            batch = NT_batch)

NT_Mvalue_ComBat_gene_pve <- calc_explained_variation(NT_Mvalue_ComBat_gene_dat, 
                                                      condition = NT_condition,
                                                      batch = NT_batch)

NT_ComBat_met_gene_pve <- calc_explained_variation(NT_ComBat_met_gene_dat, 
                                                   condition = NT_condition,
                                                   batch = NT_batch)

NT_gene_pve_df <- rbind(data.frame(type = "No adjustment", val = NT_raw_gene_pve$Batch), 
                        data.frame(type = "M-value ComBat", val = NT_Mvalue_ComBat_gene_pve$Batch),
                        data.frame(type = "ComBat-met", val = NT_ComBat_met_gene_pve$Batch)) %>%
  mutate(type = factor(type, levels = c("No adjustment",
                                        "M-value ComBat",
                                        "ComBat-met")))

# plot violin plot
p1 <- ggplot(NT_gene_pve_df, aes(x = type, y = val * 100)) +
  geom_violin(aes(fill = type)) +
  geom_boxplot(outlier.size = 0.1, width = 0.06) +
  scale_y_continuous(trans = "log10") +
  labs(y = "% Explained variation") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  ggtitle("Non-Tumor")

######################
####### tumor ########
######################
# retrieve batch and condition
TP_batch <- as.factor(sapply(samplesTP_subset_ordered, function(x) {
  str_split(x, "-")[[1]][6]
}))
TP_condition <- as.factor(c(rep(1, length(samplesTP_subset_ordered))))

# calculate % variance explained
TP_raw_gene_pve <- calc_explained_variation(TP_raw_gene_dat,
                                            condition = TP_condition,
                                            batch = TP_batch)

TP_Mvalue_ComBat_gene_pve <- calc_explained_variation(TP_Mvalue_ComBat_gene_dat, 
                                                      condition = TP_condition,
                                                      batch = TP_batch)

TP_ComBat_met_gene_pve <- calc_explained_variation(TP_ComBat_met_gene_dat, 
                                                   condition = TP_condition,
                                                   batch = TP_batch)

TP_gene_pve_df <- rbind(data.frame(type = "No adjustment", val = TP_raw_gene_pve$Batch), 
                        data.frame(type = "M-value ComBat", val = TP_Mvalue_ComBat_gene_pve$Batch),
                        data.frame(type = "ComBat-met", val = TP_ComBat_met_gene_pve$Batch)) %>%
  mutate(type = factor(type, levels = c("No adjustment",
                                        "M-value ComBat",
                                        "ComBat-met")))

# plot violin plot
p2 <- ggplot(TP_gene_pve_df, aes(x = type, y = val * 100)) +
  geom_violin(aes(fill = type)) +
  geom_boxplot(outlier.size = 0.1, width = 0.06) +
  scale_y_continuous(trans = "log10") +
  labs(y = "% Explained variation") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  ggtitle("Tumor")
