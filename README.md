
# ComBatMet

<!-- badges: start -->
<!-- badges: end -->

ComBatMet (ComBat-met) is a regression framework to adjust batch effects in DNA methylation studies. 
This package includes three batch correction methods. ComBat-met fits beta regression models to the 
beta-values, calculates batch-free distributions, and maps the quantiles of the estimated distributions 
to their batch-free counterparts. M-value ComBat log-transforms beta values to M-values prior to applying batch 
correction via ComBat. ComBat-biseq fits beta-binomial regression models to the bisulfite sequencing data 
followed by a batch correction procedure similar to ComBat-met.

## Installation

You can install ComBatMet like so:

``` r
library(devtools)
install_github("JmWangBio/ComBatMet")
```

Next, I will provide an example to demonstrate how to get started with this package.

## Example 

In this example, I will show how ComBatMet can be applied to batch correct the DNA methylation data. 
First, let’s simulate some DNA methylation data. Suppose that in a hypothetical methylation study of cancer, 
50 probes/sites are quantified. The diseased group (i.e., D) and healthy group (i.e. H) has four replicates each. 
Two replicates belong to batch 1, and the other two replicates, batch 2.

``` r
# Load library
library(ComBatMet)

# Generate a random beta-value matrix
bv_dat <- matrix(runif(n = 400, min = 0, max = 1), nrow = 50, ncol = 8)
batch <- c(rep(1, 4), rep(2, 4))
group <- rep(c(0, 1), 4)
```

Now let's adjust for batch effects.

``` r
# Adjust for batch effects including biological conditions
adj_bv_mat <- ComBat_met(bv_mat, batch = batch, group = group, full_mod = TRUE)

# Adjust for batch effects without including biological conditions
adj_bv_mat <- ComBat_met(bv_mat, batch = batch, group = group, full_mod = FALSE)
```

## Data Analysis

Scripts used for generating the simulated data and analysing the data 
from the TCGA data as shown in the manuscript are stored in the “inst” folder.

Scripts inside the "simulation" subfolder are used to generate and analyze the simulated data. 
To compare different batch correction workflows, run "dataSim_all_DE_pipeline.R" followed by 
"analyze_all_DE_data.R". To understand how sample size affects the true and false positive rates 
of ComBat-biseq, run "dataSim_biseq_lrt_pipeline.R" followed by "analyze_biseq_lrt_data.R".

Scripts inside the "TCGA" subfolder are used to download and analyze the TCGA data. 
Run "download_TCGA_data.R" to download the breast cancer subset data. "gene_all_NT_pipeline.R" and 
"gene_LumB_TP_pipeline.R" adjust the gene-level data for batch effects, while "site_all_NT_pipeline.R" 
and "site_LumB_TP_pipeline.R", the site-level data. "dea_and_pwa.R" conducts the differential expression 
analysis followed by pathway analysis. "plot_gene_box.R" makes the box plots of beta-values. "plot_gene_PCA.R" 
makes the PCA plots. "plot_gene_perc_explained_variation.R" and "plot_site_perc_explained_variation.R" make 
violin plots for % variation explained by batch in the gene-level and site-level data, respectively.
