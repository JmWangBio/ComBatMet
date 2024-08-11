
## Author: Junmin Wang
## Date: August 6th, 2024
## This script conducts differential expression analysis followed by pathway analysis.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr, tibble, limma, AnnotationDbi, and org.Hs.eg.db packages installed.
## Ensure that you also provide the correct paths for "NT_raw_gene_dat.rds", "NT_Mvalue_ComBat_gene_dat.rds", 
## "NT_ComBat_met_gene_dat.rds", "TP_raw_gene_dat.rds", "TP_Mvalue_ComBat_gene_dat.rds",
## and "TP_ComBat_met_gene_dat.rds" (lines 19 - 25), which are the intermediate output files produced from 
## "gene_all_NT_pipeline.R" and "gene_LumB_TP_pipeline.R".

## load libraries
library(dplyr)
library(tibble)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)

## load beta-value data
NT_raw_gene_dat <- readRDS("path/to/NT_raw_gene_dat.rds")
NT_Mvalue_ComBat_gene_dat <- readRDS("path/to/NT_Mvalue_ComBat_gene_dat.rds")
NT_ComBat_met_gene_dat <- readRDS("path/to/NT_ComBat_met_gene_dat.rds")

TP_raw_gene_dat <- readRDS("path/to/TP_raw_gene_dat.rds")
TP_Mvalue_ComBat_gene_dat <- readRDS("path/to/TP_Mvalue_ComBat_gene_dat.rds")
TP_ComBat_met_gene_dat <- readRDS("path/to/TP_ComBat_met_gene_dat.rds")

## combine tumor and non-tumor data
raw_gene_dat <- NT_raw_gene_dat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>%
  inner_join(TP_raw_gene_dat %>% 
               as.data.frame() %>% 
               rownames_to_column(var = "gene"), 
             by = "gene") %>%
  column_to_rownames(var = "gene") %>%
  as.matrix()

Mvalue_ComBat_gene_dat <- NT_Mvalue_ComBat_gene_dat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>%
  inner_join(TP_Mvalue_ComBat_gene_dat %>% 
               as.data.frame() %>% 
               rownames_to_column(var = "gene"), 
             by = "gene") %>%
  column_to_rownames(var = "gene") %>%
  as.matrix()

ComBat_met_gene_dat <- NT_ComBat_met_gene_dat %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>%
  inner_join(TP_ComBat_met_gene_dat %>% 
               as.data.frame() %>% 
               rownames_to_column(var = "gene"), 
             by = "gene") %>%
  column_to_rownames(var = "gene") %>%
  as.matrix()

## convert beta values to M values
raw_gene_mv <- log(raw_gene_dat / (1 - raw_gene_dat))
Mvalue_ComBat_gene_mv <- log(Mvalue_ComBat_gene_dat / (1 - Mvalue_ComBat_gene_dat))
ComBat_met_gene_mv <- log(ComBat_met_gene_dat / (1 - ComBat_met_gene_dat))

## estimate weights
group <- rep(c("normal", "cancer"), 
             times = c(ncol(NT_raw_gene_dat), ncol(TP_raw_gene_dat)))
design <- model.matrix(~ -1 + group)
var.group <- group

raw_gene_arrayw <- arrayWeights(raw_gene_mv,
                                design = design, 
                                var.group = var.group)

Mvalue_ComBat_gene_arrayw <- arrayWeights(Mvalue_ComBat_gene_mv,
                                          design = design, 
                                          var.group = var.group)

ComBat_met_gene_arrayw <- arrayWeights(ComBat_met_gene_mv,
                                       design = design, 
                                       var.group = var.group)

## diff. exp. analysis (linear model fitting)
raw_fit <- lmFit(raw_gene_mv,
                 design = design,
                 weights = raw_gene_arrayw)

Mvalue_ComBat_fit <- lmFit(Mvalue_ComBat_gene_mv,
                           design = design,
                           weights = Mvalue_ComBat_gene_arrayw)

ComBat_met_fit <- lmFit(ComBat_met_gene_mv,
                        design = design,
                        weights = ComBat_met_gene_arrayw)

## make contrasts
contrasts <- makeContrasts(contrasts = "groupcancer-groupnormal",
                           levels = design)

## eBayes
raw_contrasts_fit <- eBayes(contrasts.fit(raw_fit, contrasts))
Mvalue_ComBat_contrasts_fit <- eBayes(contrasts.fit(Mvalue_ComBat_fit, contrasts))
ComBat_met_contrasts_fit <- eBayes(contrasts.fit(ComBat_met_fit, contrasts))

## extract results
raw_res <- topTable(raw_contrasts_fit, number = Inf)
Mvalue_ComBat_res <- topTable(Mvalue_ComBat_contrasts_fit, number = Inf)
ComBat_met_res <- topTable(ComBat_met_contrasts_fit, number = Inf)

## unique genes recovered by ComBat-met only
ComBat_met_uniq_genes <- setdiff(ComBat_met_res %>% 
                                   filter(adj.P.Val < 5*1e-8) %>% 
                                   rownames(), 
                                 union(Mvalue_ComBat_res %>% 
                                         filter(adj.P.Val < 5*1e-8) %>% 
                                         rownames(), 
                                       raw_res %>% 
                                         filter(adj.P.Val < 5*1e-8) %>% 
                                         rownames()))

## pathway analysis
symbol2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                       key = ComBat_met_uniq_genes,
                                       columns = "ENTREZID",
                                       keytype = "SYMBOL")

pw_res <- limma::kegga(symbol2entrez$ENTREZID) %>%
  arrange(P.DE)
