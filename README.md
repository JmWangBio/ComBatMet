<!-- badges: start -->
[![R-CMD-check](https://github.com/JmWangBio/ComBatMet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JmWangBio/ComBatMet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

ComBatMet
================
Junmin Wang
05/25/2026

## Why Use ComBat-met?

- **Built for methylation data.** Beta-values are bounded between 0 and 1 and have skewed properties that Gaussian-based methods cannot properly handle. ComBat-met uses a beta regression framework that respects the true distributional structure of methylation data.

- **Outperforms existing methods.** ComBat-met consistently removes more batch-associated variation while preserving biological signals, achieving higher statistical power without inflating false positive rates.

- **Fast.** A dataset with 900,000 methylation sites (e.g., Illumina EPIC v2 array) across 8 samples finishes in under one minute on a single core.

- **Scalable.** ComBat-met handles 100 samples (900,000 sites) in just 4 minutes using 9 cores.

- **Works with both array and sequencing data.** ComBat-met supports beta-values from methylation arrays (e.g., Illumina EPIC) via `ComBat_met()`, and cytosine count data from NGS via `ComBat_biseq()`.

ComBat-met was published in NAR Genomics and Bioinformatics. Whenever using ComBat-met, please cite:

> Wang J. ComBat-met: Adjusting Batch Effects in DNA Methylation Data. *NAR Genom Bioinform*. 2025 May 19;7(2): lqaf062.
> DOI: [https://doi.org/10.1093/nargab/lqaf062](https://doi.org/10.1093/nargab/lqaf062)

## Installation and Usage

You can install ComBatMet like so:

``` r
library(devtools)
install_github("JmWangBio/ComBatMet")
```

Next, I will provide an example to demonstrate how to get started with this package.
First, let’s simulate some DNA methylation data. Suppose that in a hypothetical methylation 
study of cancer, 50 probes/sites are quantified. The diseased group (i.e., D) and healthy 
group (i.e., H) has four replicates each. Two replicates belong to batch 1, and the other 
two replicates, batch 2.

``` r
# Load library
library(ComBatMet)

# Generate a random beta-value matrix
bv_mat <- matrix(runif(n = 400, min = 0, max = 1), nrow = 50, ncol = 8)
batch <- c(rep(1, 4), rep(2, 4))
group <- rep(c(0, 1), 4)
```

Now let's adjust for batch effects.

``` r
# Adjust for batch effects including biological conditions
adj_bv_mat <- ComBat_met(bv_mat, batch = batch, group = group, full_mod = TRUE)
```

The example above provides a simple introduction to get you started. For more detailed examples, 
including reference-batch correction, parallelization, and other advanced features, please refer 
to the vignette.

## Introduction

Integration of genomics data is frequently impeded by technical artifacts, commonly 
referred to as batch effects. While many batch correction methods are available, they 
often fail to adequately account for the unique properties of DNA methylation data. 

To address this limitation, we developed **ComBat-met**, a novel batch adjustment framework 
specifically tailored for methylation studies. This method employs a beta regression model 
to estimate batch-free distributions, and realigns quantiles of the data to their corrected 
counterparts, providing robust adjustment for technical variations (see Fig. 1 for details of 
the workflow).

<div align="center">
  <img src="inst/images/Fig1.png" alt="Wf" width="800"/>
</div>

**Fig. 1.** Overview of the ComBat-met workflow.

## Benchmarking

We compared ComBat-met to existing batch adjustment methods using simulated methylation 
datasets with known ground truth. Each method was applied under identical conditions, and 
their ability to remove batch effects while preserving biological variability was assessed. 
Below is a description of the methods included in the comparison besides ComBat-met: 

- **M-value ComBat**: Beta values are log-transformed into M-values, and batch correction 
is performed using ComBat.
- **SVA (Surrogate Variable Analysis)**: Beta values are log-transformed into M-values, 
followed by batch correction using SVA.
- **Including Batch as a Covariate in DE**: Batch is explicitly included as a covariate 
during differential expression analysis to account for technical variation.
- **ComBat-biseq**: Uses beta-binomial regression to estimate batch-free distributions of 
bisulfite sequencing data and correct batch effects.
- **BEclear**: Uses latent factor models to detect and correct batch effects.
- **RUVm (Removing Unwanted Variation)**: Uses the two-stage RUVm approach to correct 
batch effects.

The results demonstrated that ComBat-met consistently outperformed others in reducing 
batch-induced variability while preserving biological variability (Fig. 2).

<div align="center">
  <img src="inst/images/Fig2.png" alt="benchmark" width="700"/>
</div>

**Fig. 2.** Median true positive rates and false positive rates calculated based on 
simulated data.

## Application to TCGA Data

We applied ComBat-met to methylation data from TCGA and compared its performance to other 
batch correction methods by analyzing the percentage of variation explained by batch effects. 
Results showed that ComBat-met consistently achieved the smallest percentage of batch-associated 
variation in both normal and tumor samples (Fig. 3).

<div align="center">
  <img src="inst/images/Fig3.png" alt="TCGA" width="800"/>
</div>
<br>

**Fig. 3.** Percentage of variation explained by batch.

To demonstrate the importance of batch adjustment, we also evaluated the performance of a 
**neural network** classifier both **before** and **after** applying ComBat-met. For this 
analysis:

- We selected three random methylation probes in each iteration to simulate a minimal 
and unbiased probe set, avoiding cherry-picking variables that might artificially boost 
performance.
- A feed-forward, fully connected neural network with two hidden layers was trained to 
classify normal and cancerous samples.
- Accuracy was calculated for models trained on unadjusted data and batch-adjusted data.
- By repeating this process across iterations, we observed that batch adjustment by ComBat-met 
consistently improved classification accuracy, emphasizing the importance of effective batch 
correction in methylation studies (Fig. 4).

<div align="center">
  <img src="inst/images/Fig4.png" alt="ML" width="800"/>
</div>
<br>

**Fig. 4.** Architecture of the neural network used for classification; box plot comparing 
classification accuracy before and after batch adjustment using ComBat-met.

## Code Descriptions

Scripts used for generating the simulated data and analysing the data 
from the TCGA data as shown in the manuscript are stored in the “inst” folder.

<table>
  <thead>
    <tr>
      <th>Purpose</th>
      <th>Folder</th>
      <th>File name</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan="2">
        <div>Benchmarking with simulation data</div>
      </td>
      <td rowspan="2">
        <div>simulation/</div>
      </td>
      <td>dataSim_all_DE_pipeline.R<br>dataSim_biseq_lrt_pipeline.R<br>dataSim_all_DE_pipeline_minfi.R</td>
      <td>simulation of methylation or bisulfite sequencing data followed by differential methylation analysis (limma or minfi)</td>
    </tr>
    <tr>
      <td>analyze_all_DE_data.R<br>analyze_biseq_lrt_data.R<br>analyze_all_DE_data_minfi.R</td>
      <td>calculating the true and false positive rates of each workflow</td>
    </tr>
    <tr>
      <td rowspan="8">
        <div>Benchmarking with TCGA data</div>
      </td>
      <td rowspan="5">
        <div>TCGA/</div>
      </td>
      <td>download_TCGA_data.R</td>
      <td>downloading the breast cancer subtype data</td>
    </tr>
    <tr>
      <td>site_all_NT_pipeline.R<br>site_LumB_TP_pipeline.R</td>
      <td>adjusting the site-level data for batch effects in normal and tumor samples</td>
    </tr>
    <tr>
      <td>plot_site_box.R</td>
      <td>making box plots of beta-values</td>
    </tr>
    <tr>
      <td>plot_site_pca.R</td>
      <td>making PCA plots for the site-level data</td>
    </tr>
    <tr>
      <td>plot_site_perc_explained_variation.R</td>
      <td>making violin plots for % variation explained by batch in the site-level data</td>
    </tr>
    <tr>
      <td rowspan="3">
        <div>ML/</div>
      </td>
      <td>TCGA_extract_metadata.R</td>
      <td>extracting metadata for samples (normal or cancerous)</td>
    </tr>
    <tr>
      <td>ML.py</td>
      <td>evaluating the impact of batch adjustment on classification by repeatedly selecting random probes, training a neural network, and comparing accuracies before and after batch correction</td>
    </tr>
    <tr>
      <td>visualization.R</td>
      <td>making a box plot to visualize differences in accuracy before and after batch adjustment</td>
    </tr>
  </tbody>
</table>

## Common Q & A

**Q: Can I provide M-values as input?**

A: Yes. ComBat-met accepts both beta-values and M-values. Set `dtype = "M-value"` when calling `ComBat_met()`, and the output will also be M-values.

**Q: Does ComBat-met work on unbalanced experimental designs?**

A: Yes. ComBat-met handles unbalanced designs where batches differ in sample size.

**Q: How long does ComBat-met take to run?**

A: Runtime depends on the size of your dataset and the computing resources available. A dataset with 900,000 methylation sites (e.g., Illumina EPIC v2 array) across 8 samples completes in under one minute on a single core. You can further reduce runtime via parallelisation by setting `ncores` to the number of available cores — 900,000 sites across 100 samples finishes in approximately 4 minutes using 9 cores.

**Q: Can I use ComBat-met to adjust sequencing read counts directly?**

A: Yes, via `ComBat_biseq()`, which uses a beta-binomial regression framework that explicitly models the uncertainty from variable read depth. A dataset with 1 million methylation sites across 8 samples finishes in approximately 4 minutes on a single core; 1 million sites across 100 samples (10 batches) finishes in approximately 11 minutes using 9 cores.

**Q: Is ComBat-met compatible with RNA methylation data?**

A: Yes, `ComBat_met()` can be applied to RNA methylation data. However, keep in mind that unlike DNA methylation, RNA methylation coverage varies widely because it is tied to transcript abundance. In such cases, `ComBat_biseq()` may be better suited.

**Q: How does ComBat-met handle beta regression?**

A: The current implementation uses C++ routines with LAPACK/BLAS kernels, inspired by the GLM regression framework in edgeR and ComBat-seq. Precision is estimated within each batch via a fast grid-based search, replacing the original `betareg`-based engine. As demonstrated in the vignettes, adjusted values from the old and new implementations are 99% correlated, confirming that the speedup comes with no meaningful loss of accuracy. If you prefer to use the original `betareg` implementation for any reason, you can do so by setting `use_cpp = FALSE`. `ComBat_biseq()` handles beta-binomial regression in a similar way.

**Q: Does ComBat-met work across operating systems?**

A: Yes. ComBat-met is tested on macOS, Windows, and Linux via GitHub Actions CI, and passes checks on all three platforms.

## Contribute

We welcome your feedback and contributions! If you have suggestions, find an issue, or want to add new features, feel free to open an issue or submit a pull request. Together, we can make this project even better!
