
## Author: Junmin Wang
## Date: July 18th, 2024
## This script can be used for differential methylation analysis.

diff_met_wf <- function(bv, group, contrast.pair, 
                        pseudo_beta = 1e-4, eb = TRUE,
                        batch = NULL, sva = FALSE) {
  ## convert extreme 0 or 1 values to pseudo-beta
  if (pseudo_beta <= 0 | pseudo_beta >= 0.5) {
    stop("Invalid pseudo beta-values.")
  }
  bv[bv <= 0] <- pseudo_beta
  bv[bv >= 1] <- 1 - pseudo_beta
  
  ## convert beta values to M values
  mv <- log(bv / (1 - bv))
  
  ## create design matrix
  group <- factor(group)
  if (is.null(batch)) {
    design <- model.matrix(~ 0 + group)
  } else {
    batch = factor(batch)
    design <- model.matrix(~ 0 + group + batch)
  }
  
  ## run sva if requested
  if (sva & is.null(batch)) {
    svobj <- sva::sva(mv, model.matrix(~group))
    if (svobj$n.sv > 0) {
      sv <- svobj$sv
      colnames(sv) <- paste0("sv", 1:ncol(sv))
      design <- cbind(design, sv)
    }
  }
  
  ## fit linear model
  fit <- limma::lmFit(mv, design)
  
  ## create contrasts
  contrasts <- limma::makeContrasts(contrasts = contrast.pair,
                                    levels = design)
  
  ## empirical Bayes shrinkage
  contrasts.fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts))
  
  ## calculate moderated or ordinary t-statistics, p-values, and adj. p-values
  if (eb) {
    contrasts.fit.res <- limma::topTable(contrasts.fit, coef = 1,
                                         number = Inf, sort.by = "none")
    t_stat <- contrasts.fit.res$t
    pval <- contrasts.fit.res$P.Value
    adj.pval <- contrasts.fit.res$adj.P.Val
  } else {
    t_stat <- contrasts.fit$coef[, 1] / contrasts.fit$stdev.unscaled[, 1] / 
      contrasts.fit$sigma
    pval <- 2 * pt(abs(t_stat), df = fit$df.residual, lower.tail = FALSE)
    adj.pval <- p.adjust(pval, method = "BH")
  }
  
  ## create output
  res <- data.frame(id = 1:nrow(mv),
                    t = t_stat,
                    pval = pval,
                    adj.pval = adj.pval)
  return(res)
}
