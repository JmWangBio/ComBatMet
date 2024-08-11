
## Author: Junmin Wang
## Date: July 25th, 2024
## This script queries, downloads, and prepares the TCGA breast cancer methylation data.
## IMPORTANT NOTE: To run this script successfully, you need to have the TCGAbiolinks R package installed.
## Make sure to change the paths of outputs to where you want to save them (lines 18, 22, and 24).

## load libraries
library(TCGAbiolinks)

## query BRCA methylation data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "DNA Methylation",
                  platform = "Illumina Human Methylation 450",
                  data.type = "Methylation Beta Value")

## download data
GDCdownload(query, directory = "path/to/output")

## prepare data in the right format
data <- GDCprepare(query, save = TRUE,
                   save.filename = "path/to/BRCA_met.RData",
                   summarizedExperiment = TRUE,
                   directory = "path/to/output",
                   remove.files.prepared = FALSE)
