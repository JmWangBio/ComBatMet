
## Author: Junmin Wang
## Date: January 8th, 2025
## This script makes box plots of beta values in non-tumor and tumor data.
## IMPORTANT NOTE: To run this script successfully, you need to have the reshape2, stringr, and ggplot2 packages installed.
## Ensure that you also provide the correct paths for "NT_raw_site_dat.rds" and "TP_raw_site_dat.rds" (lines 15 - 16), 
## which are the intermediate output files produced from "site_all_NT_pipeline.R" and "site_LumB_TP_pipeline.R".

## load libraries
library(reshape2)
library(stringr)
library(ggplot2)

## load beta-value data
NT_raw_site_dat <- readRDS("path/to/NT_raw_site_dat.rds")
TP_raw_site_dat <- readRDS("path/to/TP_raw_site_dat.rds")

#########################################
######## box plot of beta values ########
#########################################
## non-tumor
NT_raw_site_dat_long <- melt(NT_raw_site_dat)
colnames(NT_raw_site_dat_long) <- c("site", "sample", "value")
palette <- c("salmon", "lightblue")
NT_col_vec <- palette[as.numeric(as.factor(sapply(unique(NT_raw_site_dat_long$sample), function(x) {
  str_split(x, "-")[[1]][6]
}))) %% 2 + 1]

# box plot
p1 <- ggplot(NT_raw_site_dat_long, 
             aes(x = sample, y = value, 
                 fill = sample)) +
  geom_boxplot(linetype = "dashed") +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), 
               outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black")) +
  labs(y = "% Methylated") +
  ggtitle("Non-Tumor") +
  scale_fill_manual(breaks = unique(NT_raw_site_dat_long$sample),
                    values = NT_col_vec)

## tumor
TP_raw_site_dat_long <- melt(TP_raw_site_dat)
colnames(TP_raw_site_dat_long) <- c("site", "sample", "value")
palette <- c("salmon", "lightblue")
TP_col_vec <- palette[as.numeric(as.factor(sapply(unique(TP_raw_site_dat_long$sample), function(x) {
  str_split(x, "-")[[1]][6]
}))) %% 2 + 1]

# box plot
p2 <- ggplot(TP_raw_site_dat_long, 
             aes(x = sample, y = value, 
                 fill = sample)) +
  geom_boxplot(linetype = "dashed") +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), 
               outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.title.y = element_text(color = "black")) +
  labs(y = "% Methylated") +
  ggtitle("Tumor") +
  scale_fill_manual(breaks = unique(TP_raw_site_dat_long$sample),
                    values = TP_col_vec) 
