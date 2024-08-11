
## Author: Junmin Wang
## Date: August 6th, 2024
## This script plots the PCAs in non-tumor and tumor data.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tibble, stringr, and ggplot2 packages installed.
## Ensure that you also provide the correct paths for "NT_raw_gene_dat.rds", "NT_Mvalue_ComBat_gene_dat.rds", 
## "NT_ComBat_met_gene_dat.rds", "TP_raw_gene_dat.rds", "TP_Mvalue_ComBat_gene_dat.rds",
## and "TP_ComBat_met_gene_dat.rds" (lines 18 - 24), which are the intermediate output files produced from 
## "gene_all_NT_pipeline.R" and "gene_LumB_TP_pipeline.R".

## load libraries
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)

## load beta-value data
NT_raw_gene_dat <- readRDS("path/to/NT_raw_gene_dat.rds")
NT_Mvalue_ComBat_gene_dat <- readRDS("path/to/NT_Mvalue_ComBat_gene_dat.rds")
NT_ComBat_met_gene_dat <- readRDS("path/to/NT_ComBat_met_gene_dat.rds")

TP_raw_gene_dat <- readRDS("path/to/TP_raw_gene_dat.rds")
TP_Mvalue_ComBat_gene_dat <- readRDS("path/to/TP_Mvalue_ComBat_gene_dat.rds")
TP_ComBat_met_gene_dat <- readRDS("path/to/TP_ComBat_met_gene_dat.rds")

######################
######## PCA #########
######################
## non-tumor
# NT raw data
pca.res <- prcomp(t(NT_raw_gene_dat))
p1 <- ggplot(as.data.frame(pca.res$x) %>%
               rownames_to_column(var = "ID") %>%
               mutate(plate = sapply(ID, function(x) {
                 str_split(x, "-")[[1]][6]
               })),
             aes(x = PC1, y = PC2, color = plate)) +
  geom_point(size = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.margin = unit(c(0.2, 0.6, 0.2, 0.6), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  labs(shape = "Group", color = "Plate",
       x = sprintf("PC1 (%s%%)", round(pca.res$sdev[1]^2 / sum(pca.res$sdev^2) * 100, 2)),
       y = sprintf("PC2 (%s%%)", round(pca.res$sdev[2]^2 / sum(pca.res$sdev^2) * 100, 2))) +
  scale_color_manual(breaks = c("A093", "A10Q", "A12R", "A138", "A13T", "A145", 
                                "A14H", "A14N", "A161", "A16A", "A17F"),
                     values = c("black", "red", "yellow", "green", "blue", "purple",
                                "pink", "gray", "brown", "orange", "cyan")) +
  ggtitle("No adjustment")

# NT M-value ComBat
pca.res <- prcomp(t(NT_Mvalue_ComBat_gene_dat))
p2 <- ggplot(as.data.frame(pca.res$x) %>%
               rownames_to_column(var = "ID") %>%
               mutate(plate = sapply(ID, function(x) {
                 str_split(x, "-")[[1]][6]
               })),
             aes(x = PC1, y = PC2, color = plate)) +
  geom_point(size = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.margin = unit(c(0.2, 0.6, 0.2, 0.6), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  labs(shape = "Group", color = "Plate",
       x = sprintf("PC1 (%s%%)", round(pca.res$sdev[1]^2 / sum(pca.res$sdev^2) * 100, 2)),
       y = sprintf("PC2 (%s%%)", round(pca.res$sdev[2]^2 / sum(pca.res$sdev^2) * 100, 2)))+
  scale_color_manual(breaks = c("A093", "A10Q", "A12R", "A138", "A13T", "A145", 
                                "A14H", "A14N", "A161", "A16A", "A17F"),
                     values = c("black", "red", "yellow", "green", "blue", "purple",
                                "pink", "gray", "brown", "orange", "cyan")) +
  ggtitle("M-value ComBat")

# NT ComBat-Met
pca.res <- prcomp(t(NT_ComBat_met_gene_dat))
p3 <- ggplot(as.data.frame(pca.res$x) %>%
               rownames_to_column(var = "ID") %>%
               mutate(plate = sapply(ID, function(x) {
                 str_split(x, "-")[[1]][6]
               })),
             aes(x = PC1, y = PC2, color = plate)) +
  geom_point(size = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.margin = unit(c(0.2, 0.6, 0.2, 0.6), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  labs(shape = "Group", color = "Plate",
       x = sprintf("PC1 (%s%%)", round(pca.res$sdev[1]^2 / sum(pca.res$sdev^2) * 100, 2)),
       y = sprintf("PC2 (%s%%)", round(pca.res$sdev[2]^2 / sum(pca.res$sdev^2) * 100, 2)))+
  scale_color_manual(breaks = c("A093", "A10Q", "A12R", "A138", "A13T", "A145", 
                                "A14H", "A14N", "A161", "A16A", "A17F"),
                     values = c("black", "red", "yellow", "green", "blue", "purple",
                                "pink", "gray", "brown", "orange", "cyan")) +
  ggtitle("ComBat-met")

## tumor
# TP raw data
pca.res <- prcomp(t(TP_raw_gene_dat))
p4 <- ggplot(as.data.frame(pca.res$x) %>%
               rownames_to_column(var = "ID") %>%
               mutate(plate = sapply(ID, function(x) {
                 str_split(x, "-")[[1]][6]
               })),
             aes(x = PC1, y = PC2, color = plate)) +
  geom_point(size = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.margin = unit(c(0.2, 0.6, 0.2, 0.6), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  labs(shape = "Group", color = "Plate",
       x = sprintf("PC1 (%s%%)", round(pca.res$sdev[1]^2 / sum(pca.res$sdev^2) * 100, 2)),
       y = sprintf("PC2 (%s%%)", round(pca.res$sdev[2]^2 / sum(pca.res$sdev^2) * 100, 2)))+
  scale_color_manual(breaks = c("A10A", "A10N", "A10P", "A12R", "A138", "A13K", 
                                "A145", "A14H", "A14N", "A161", "A16A", "A16G",
                                "A17Z", "A212", "A21R", "A268", "A28C", "A357",
                                "A41Q"),
                     values = c("black", "red", "yellow", "green", "blue", "purple",
                                "pink", "gray", "brown", "orange", "cyan", "salmon",
                                "darkgreen", "violet", "lightblue", "gold", "navy", "yellowgreen",
                                "firebrick")) +
  ggtitle("No adjustment")

# TP M-value ComBat
pca.res <- prcomp(t(TP_Mvalue_ComBat_gene_dat))
p5 <- ggplot(as.data.frame(pca.res$x) %>%
               rownames_to_column(var = "ID") %>%
               mutate(plate = sapply(ID, function(x) {
                 str_split(x, "-")[[1]][6]
               })),
             aes(x = PC1, y = PC2, color = plate)) +
  geom_point(size = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.margin = unit(c(0.2, 0.6, 0.2, 0.6), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  labs(shape = "Group", color = "Plate",
       x = sprintf("PC1 (%s%%)", round(pca.res$sdev[1]^2 / sum(pca.res$sdev^2) * 100, 2)),
       y = sprintf("PC2 (%s%%)", round(pca.res$sdev[2]^2 / sum(pca.res$sdev^2) * 100, 2)))+
  scale_color_manual(breaks = c("A10A", "A10N", "A10P", "A12R", "A138", "A13K", 
                                "A145", "A14H", "A14N", "A161", "A16A", "A16G",
                                "A17Z", "A212", "A21R", "A268", "A28C", "A357",
                                "A41Q"),
                     values = c("black", "red", "yellow", "green", "blue", "purple",
                                "pink", "gray", "brown", "orange", "cyan", "salmon",
                                "darkgreen", "violet", "lightblue", "gold", "navy", "yellowgreen",
                                "firebrick")) +
  ggtitle("M-value ComBat")

# TP ComBat-Met
pca.res <- prcomp(t(TP_ComBat_met_gene_dat))
p6 <- ggplot(as.data.frame(pca.res$x) %>%
               rownames_to_column(var = "ID") %>%
               mutate(plate = sapply(ID, function(x) {
                 str_split(x, "-")[[1]][6]
               })),
             aes(x = PC1, y = PC2, color = plate)) +
  geom_point(size = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.margin = unit(c(0.2, 0.6, 0.2, 0.6), "cm"),
        plot.title = element_text(hjust = 0.5)) +
  labs(shape = "Group", color = "Plate",
       x = sprintf("PC1 (%s%%)", round(pca.res$sdev[1]^2 / sum(pca.res$sdev^2) * 100, 2)),
       y = sprintf("PC2 (%s%%)", round(pca.res$sdev[2]^2 / sum(pca.res$sdev^2) * 100, 2)))+
  scale_color_manual(breaks = c("A10A", "A10N", "A10P", "A12R", "A138", "A13K", 
                                "A145", "A14H", "A14N", "A161", "A16A", "A16G",
                                "A17Z", "A212", "A21R", "A268", "A28C", "A357",
                                "A41Q"),
                     values = c("black", "red", "yellow", "green", "blue", "purple",
                                "pink", "gray", "brown", "orange", "cyan", "salmon",
                                "darkgreen", "violet", "lightblue", "gold", "navy", "yellowgreen",
                                "firebrick")) +
  guides(color = guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("ComBat-met")
