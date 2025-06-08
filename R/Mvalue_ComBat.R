#' Adjust for batch effects in DNA methylation data by converting beta-values to M-values followed by ComBat
#' 
#' Mvalue_ComBat converts beta-values to M-values through logit-transformation, 
#' adjusts M-values for batch effects using ComBat, and converts the adjusted M-values back to 
#' adjusted beta-values through reverse logit-transformation.
#' Forward and reverse logit-transformation are omitted if M-values are provided.
#' 
#' @inheritParams ComBat_met
#'
#' @return \code{Mvalue_ComBat} returns a feature x sample matrix with the same data type as input, adjusted for batch effects.
#' @export
#'
#' @examples
#' # Generate a random beta-value matrix
#' bv_mat <- matrix(runif(n = 400, min = 0, max = 1), nrow = 50, ncol = 8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#'
#' # Adjust for batch effects including biological conditions
#' adj_bv_mat <- Mvalue_ComBat(bv_mat, dtype = "b-value", batch = batch, 
#' group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_bv_mat <- Mvalue_ComBat(bv_mat, dtype = "b-value", batch = batch, 
#' group = group, full_mod = FALSE)
#' 
#' # Generate a random M-value matrix
#' mv_mat <- matrix(rnorm(n = 400, mean = 0, sd = 1), nrow = 50, ncol = 8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#'
#' # Adjust for batch effects including biological conditions
#' adj_mv_mat <- Mvalue_ComBat(mv_mat, dtype = "M-value", batch = batch, 
#' group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_mv_mat <- Mvalue_ComBat(mv_mat, dtype = "M-value", batch = batch, 
#' group = group, full_mod = FALSE)
#' 
#' # Store data as a SummarizedExperiment object
#' library(SummarizedExperiment)
#' colData <- data.frame(batch = batch, group = group)
#' dat.nn <- bv_mat
#' dat.se <- SummarizedExperiment(assays = list(val = dat.nn), colData = colData)
#' adj_bv_mat <- Mvalue_ComBat(dat.se, dtype = "b-value", full_mod = TRUE)
#' 

Mvalue_ComBat <- function(vmat, dtype = "b-value", 
                          batch, group = NULL, 
                          covar_mod = NULL, full_mod = TRUE, 
                          mean.only = FALSE, pseudo_beta = 1e-4,
                          ref.batch = NULL) {
  ## check if vmat is a SummarizedExperiment object
  if (inherits(vmat, "SummarizedExperiment")) {
    ## Get the metadata (column data)
    col_data <- SummarizedExperiment::colData(vmat)
    
    ## Check for the batch column
    if ("batch" %in% colnames(col_data)) {
      batch <- col_data[, "batch", drop = TRUE]
    } else {
      stop("Metadata must contain a batch column.")
    }
    
    ## Check for the group column
    if ("group" %in% colnames(col_data)) {
      group <- col_data[, "group", drop = TRUE]
    } else {
      cat("Metadata does not have a group column.\n")
      group <- NULL
    }
    
    ## Check for covariates
    all_cols <- colnames(col_data)
    covar_cols <- setdiff(all_cols, c("batch", "group"))
    if (length(covar_cols) > 0) {
      covar_mod <- col_data[, covar_cols, drop = FALSE]
    } else {
      cat("Metadata does not have other covariate columns besides batch or group.\n")
      covar_mod <- NULL
    }    
    
    ## Extract assay data
    if (length(SummarizedExperiment::assays(vmat)) == 0) {
      stop("No assay data found in the SummarizedExperiment object.")
    }
    vmat <- SummarizedExperiment::assay(vmat)
  }
  
  ## convert extreme 0 or 1 values to pseudo-beta
  if (dtype == "b-value") {
    if (pseudo_beta <= 0 | pseudo_beta >= 0.5) {
      stop("Invalid pseudo beta-values.")
    }
    vmat[vmat == 0] <- pseudo_beta
    vmat[vmat == 1] <- 1 - pseudo_beta
  }
  
  ## convert beta values to M values if needed
  if (dtype == "b-value") {
    mv <- log(vmat / (1 - vmat))
  } else {
    mv <- vmat
  }
  
  ## construct the model matrix
  mod <- covar_mod
  if (full_mod) {
    mod <- cbind(group, mod)
  }
  
  ## run ComBat
  if (mean.only) {
    adj_mv <- sva::ComBat(mv, batch, mod = mod, mean.only = TRUE, 
                     par.prior = TRUE, ref.batch = ref.batch)
  } else {
    adj_mv <- sva::ComBat(mv, batch, mod = mod, mean.only = FALSE, 
                     par.prior = TRUE, ref.batch = ref.batch)
  }
  
  ## convert adjusted M values back to beta values
  if (dtype == "b-value") {
    adj_vmat <- exp(adj_mv) / (1 + exp(adj_mv))
  } else {
    adj_vmat <- adj_mv
  }
  return(adj_vmat)
}
