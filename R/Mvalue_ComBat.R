#' Adjust for batch effects in DNA methylation data by converting beta-values to M-values followed by ComBat
#' 
#' Mvalue_ComBat converts beta-values to M-values through logit-transformation, 
#' adjusts M-values for batch effects using ComBat, and converts the adjusted M-values back to 
#' adjusted beta-values through reverse logit-transformation.
#' 
#' @inheritParams ComBat_met
#'
#' @return \code{Mvalue_ComBat} returns a feature x sample beta-value matrix, adjusted for batch effects.
#' @export
#'
#' @examples
#' # Generate a random beta-value matrix
#' bv_mat <- matrix(runif(n = 400, min = 0, max = 1), nrow = 50, ncol = 8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#'
#' # Adjust for batch effects including biological conditions
#' adj_bv_mat <- Mvalue_ComBat(bv_mat, batch = batch, group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_bv_mat <- Mvalue_ComBat(bv_mat, batch = batch, group = group, full_mod = FALSE)
#' 

Mvalue_ComBat <- function(bv, batch, group = NULL, 
                          covar_mod = NULL, full_mod = TRUE, 
                          mean.only = FALSE, pseudo_beta = 1e-4,
                          ref.batch = NULL) {
  ## convert extreme 0 or 1 values to pseudo-beta
  if (pseudo_beta <= 0 | pseudo_beta >= 0.5) {
    stop("Invalid pseudo beta-values.")
  }
  bv[bv == 0] <- pseudo_beta
  bv[bv == 1] <- 1 - pseudo_beta
  
  ## convert beta values to M values
  mv <- log(bv / (1 - bv))
  
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
  adj_bv <- exp(adj_mv) / (1 + exp(adj_mv))
  return(adj_bv)
}
