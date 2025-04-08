#' Adjust for batch effects using a beta regression framework in DNA methylation data
#'
#' ComBat-met fits beta regression models to the beta-values or M-values, 
#' calculates batch-free distributions, and maps the quantiles of the estimated distributions 
#' to their batch-free counterparts.
#'
#' @param vmat matrix of beta-values or M-values
#' @param dtype data type: b-value or M-value; note that the input and output have the same data type.
#' @param batch vector for batch
#' @param group optional vector for biological condition of interest
#' @param covar_mod optional model matrix representing co-variates to be included in the model
#' @param full_mod Boolean variable indicating whether to include biological condition of interest in the model
#' @param shrink Boolean variable indicating whether to apply EB-shrinkage on parameter estimation
#' @param mean.only Boolean variable indicating whether to apply EB-shrinkage on the estimation of precision effects
#' @param feature.subset.n number of features to use in non-parametric EB estimation, only useful when shrink equals TRUE
#' @param pseudo_beta pseudo beta-values to be used for replacing extreme 0 and 1 beta-values. 
#' Value needs to be between 0 and 0.5. Only active when dtype equals b-value.
#' @param ref.batch NULL by default. If given, that batch will be selected as a reference for batch correction.
#' @param ncores number of cores to be used for parallel computing. By default, ncores is set to one.
#' 
#' @return \code{ComBat_met} returns a feature x sample matrix with the same data type as input, adjusted for batch effects.
#' @export
#' 
#' @examples
#' # Generate a random beta-value matrix
#' bv_mat <- matrix(runif(n = 400, min = 0, max = 1), nrow = 50, ncol = 8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#'
#' # Adjust for batch effects including biological conditions
#' adj_bv_mat <- ComBat_met(bv_mat, dtype = "b-value", batch = batch, group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_bv_mat <- ComBat_met(bv_mat, dtype = "b-value", batch = batch, group = group, full_mod = FALSE)
#'
#' # Generate a random M-value matrix
#' mv_mat <- matrix(rnorm(n = 400, mean = 0, sd = 1), nrow = 50, ncol = 8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#'
#' # Adjust for batch effects including biological conditions
#' adj_mv_mat <- ComBat_met(mv_mat, dtype = "M-value", batch = batch, group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_mv_mat <- ComBat_met(mv_mat, dtype = "M-value", batch = batch, group = group, full_mod = FALSE)
#' 
#' # Adjust for batch effects including biological conditions (multiple threads)
#' adj_mv_mat <- ComBat_met(mv_mat, dtype = "M-value", batch = batch, group = group, full_mod = TRUE, 
#' ncores = 2)
#' 

ComBat_met <- function(vmat, dtype = "b-value", 
                       batch, group = NULL, covar_mod = NULL, full_mod = TRUE,
                       shrink = FALSE, mean.only = FALSE, feature.subset.n = NULL,
                       pseudo_beta = 1e-4, ref.batch = NULL, ncores = 1) {
  ########  Preparation  ########
  ## check if vmat has the correct format
  if (!(is.matrix(vmat) && is.numeric(vmat))) {
    if (is.data.frame(vmat) && all(sapply(vmat, is.numeric))) {
      vmat <- as.matrix(vmat)
    } else {
      stop("vmat must be a matrix of beta-values or M-values.")
    }
  }
  
  ## convert extreme 0 or 1 values to pseudo-beta
  if (dtype == "b-value") {
    if (pseudo_beta <= 0 | pseudo_beta >= 0.5) {
      stop("Invalid pseudo beta-values.")
    }
    vmat[vmat == 0] <- pseudo_beta
    vmat[vmat == 1] <- 1 - pseudo_beta    
  }
  
  ## check if coverage matrix and batch have matching sizes
  if (ncol(vmat) != length(batch)) {
    stop("Coverage matrix and batch vector do not have matching sizes.")
  }
  
  ## check if coverage matrix and group have matching sizes
  if (!is.null(group)) {
    if (ncol(vmat) != length(group)) {
      stop("Coverage matrix and group vector do not have matching sizes.")
    }
  }
  
  ## Does not allow beta values outside the (0, 1) range
  if (dtype == "b-value") {
    if (sum(vmat >= 1, na.rm = TRUE) > 0 | sum(vmat <= 0, na.rm = TRUE) > 0) {
      stop("All beta values must be between 0 and 1.")
    }    
  }
  
  ## Does not support more than 1 batch variable
  if (length(dim(batch)) > 1) {
    stop("ComBat-met does not allow more than one batch variable!")
  }
  
  ## Does not support 1 sample across all batches
  batch <- as.factor(batch)
  if (all(table(batch) <= 1)) {
    stop("ComBat-met doesn't support only 1 sample across all batches!")
  }
  
  ## Does not support only 1 batch level 
  if (length(levels(batch)) <= 1) {
    stop("Found only one batch. No need to adjust for batch effects!")
  }
  
  ## Correct for mean batch effects only if any batch has only 1 sample
  if (any(table(batch) == 1)) {
    cat("At least one batch contains only 1 sample. Only mean batch effects will be corrected.\n")
    mean.only <- TRUE
  }
  
  ## Remove features with zero variance across all batches
  zero.var.rows.lst <- lapply(levels(batch), function(b) {
    which(apply(vmat[, batch == b], 1, function(x) {stats::var(x, na.rm = TRUE) == 0}))
  })
  all.zero.var.rows <- Reduce(intersect, zero.var.rows.lst)
  
  if (length(all.zero.var.rows) > 0) {
    cat(sprintf("Found %s features with uniform values across all batches; 
                these features will not be adjusted for batch effects.\n",
                length(all.zero.var.rows)))
  }
  
  keep <- setdiff(1:nrow(vmat), all.zero.var.rows)
  vmatOri <- vmat
  vmat <- vmatOri[keep, ]
  
  ## Create a vector for correction types
  if (mean.only) {
    mean.only.vec <- rep(TRUE, length(keep))
  } else {
    mean.only.vec <- rep(FALSE, length(keep))
  }
  
  ## Calculate batch-related statistics
  # number of batches
  n_batch <- nlevels(batch)
  # list of samples in each batch
  batches_ind <- lapply(1:n_batch, function(i) {
    which(batch == levels(batch)[i])
  })
  # number of samples in each batch
  n_batches <- sapply(batches_ind, length)
  # total number of samples
  n_sample <- sum(n_batches)
  cat(sprintf("Found %s batches and %s samples.\n", n_batch, n_sample))
  
  ## Make design matrix
  # biological condition matrix
  group <- as.factor(group)
  if (full_mod & nlevels(group) > 1) {
    cat("Using full model in ComBat-met.\n")
    mod <- stats::model.matrix(~group)
  } else {
    cat("Using null model in ComBat-met.\n")
    mod <- stats::model.matrix(~1, data = as.data.frame(t(vmat)))
  }
  # covariate matrix
  if (!is.null(covar_mod)) {
    if (is.data.frame(covar_mod)) {
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i) {
        stats::model.matrix(~covar_mod[, i])
      }))
    }
  }
  # combine covariate matrix with biological condition matrix
  mod <- cbind(mod, covar_mod)
  # combine with batch matrix
  batchmod <- stats::model.matrix(~-1 + batch)
  if (!is.null(ref.batch)) {
    ## check for reference batch and make appropriate changes
    if (!(ref.batch %in% levels(batch))) {
      stop("Reference level ref. batch is not one of the levels of the batch variable.")
    }
    cat(sprintf("Using batch %s as the reference batch\n", ref.batch))
    ref <- which(levels(batch) == ref.batch)
  } else {
    ref <- NULL
  }
  
  design <- cbind(batchmod, mod)
  
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  cat("Adjusting for", ncol(design) - ncol(batchmod), 'covariate(s) or covariate level(s).\n')
  
  ## Check if the design is confounded
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n_batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-met.")
    }
    if (ncol(design) > (n_batch+1)) {
      if ((qr(design[, -c(1:n_batch)])$rank < ncol(design[, -c(1:n_batch)]))) {
        stop('The covariates are confounded! Please remove one or more of the covariates 
             so the design is not confounded.')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded 
             covariates and rerun ComBat-met.")
      }
    }
  }
  
  ## Check for missing values in count matrix
  NAs = any(is.na(vmat))
  if (NAs) {cat(c('Found', sum(is.na(vmat)), 'missing data values\n'), sep = ' ')}
  
  ## convert M-values to beta-values if needed
  if (dtype == "b-value") {
    bv <- vmat
  } else {
    bv <- exp(vmat) / (1 + exp(vmat))
  }
  
  ########  Estimate parameters from beta GLM  ########
  cat("Fitting the GLM model\n")
  
  if (!is.numeric(ncores) || ncores != as.integer(ncores) || ncores <= 0) {
    stop("ncores must be a positive integer.")
  }
  num_cores <- max(1, parallel::detectCores() - 1)
  num_cores <- min(ncores, num_cores)
  cl <- parallel::makeCluster(num_cores)
  
  # define the function to be applied in parallel
  param_estim <- function(k) {
    result<- list(gamma_hat = NULL, mu_hat = NULL, phi_hat = NULL, delta_hat = NULL,
                  zero_modvar = 0, zero_modvar_batch = 0, moderr = 0)
    
    # mark rows with NA values
    full_mat <- cbind(design, bv[k, ])
    nona <- which(stats::complete.cases(full_mat))
    
    # check if the data are all NAs
    if (length(nona) == 0) {
      result$zero_modvar <- 1
      return(result)
    }
    
    # check if the model has zero model variance
    if (qr(full_mat[nona, ])$rank < ncol(full_mat)) {
      result$zero_modvar <- 1
      return(result)
    }
    
    # if precision correction enabled, check whether the model has zero model 
    # variance within any batch
    if (!mean.only.vec[k]) {
      for (i in 1:length(batches_ind)) {
        if (qr(full_mat[intersect(batches_ind[[i]], nona), c(i, (n_batch+1):ncol(full_mat))])$rank < 
            ncol(full_mat) - n_batch + 1) {
          result$zero_modvar_batch <- 1
          return(result)
        }
      }
    }
    
    # model fit
    if (mean.only.vec[k]) {
      glm_f <- tryCatch({
        betareg::betareg.fit(x = design[nona, ], y = bv[k, ][nona])
      }, error = function(e) {
        e
      })
    } else {
      glm_f <- tryCatch({
        betareg::betareg.fit(x = design[nona, ], y = bv[k, ][nona], 
                             z = batchmod[nona, ])
      }, error = function(e) {
        e
      })
    }
    
    # if error with model fitting
    if (inherits(glm_f, "error")) {
      result$moderr <- 1
      return(result)
    }
    
    # compute mean and precision intercepts as batch-size-weighted average from batches
    if (!is.null(ref.batch)) {
      alpha_x <- glm_f$coefficients$mean[ref]
    } else {
      alpha_x <- glm_f$coefficients$mean[1:n_batch] %*% 
        as.matrix(colSums(batchmod[nona, ]) / length(nona))
    }
    
    if (mean.only.vec[k]) {
      alpha_z <- glm_f$coefficients$precision
    } else {
      if (!is.null(ref.batch)) {
        alpha_z <- glm_f$coefficients$precision[ref]
      } else {
        alpha_z <- glm_f$coefficients$precision %*%
          as.matrix(colSums(batchmod[nona, ]) / length(nona))
      }
    }
    
    # estimate parameters
    gamma_hat <- glm_f$coefficients$mean[1:n_batch] - as.numeric(alpha_x)
    mu_hat <- rep(NA, nrow(full_mat))
    mu_hat[nona] <- glm_f$fitted.values
    phi_hat <- as.numeric(exp(alpha_z)) * rep(1, nrow(full_mat))
    if (mean.only.vec[k]) {
      delta_hat <- rep(0, n_batch)
    } else {
      delta_hat <- glm_f$coefficients$precision - as.numeric(alpha_z)
    }
    
    # store result
    result$gamma_hat <- gamma_hat
    result$mu_hat <- mu_hat
    result$phi_hat <- phi_hat
    result$delta_hat <- delta_hat
    return(result)
  }
  
  # run in parallel
  result_lst <- parallel::parLapply(cl, 1:nrow(bv), param_estim)
  parallel::stopCluster(cl)
  
  gamma_hat_lst <- lapply(result_lst, function(x) x$gamma_hat)
  mu_hat_lst <- lapply(result_lst, function(x) x$mu_hat)
  phi_hat_lst <- lapply(result_lst, function(x) x$phi_hat)
  delta_hat_lst <- lapply(result_lst, function(x) x$delta_hat)
  n_zero_modvar <- sum(unlist(lapply(result_lst, function(x) x$zero_modvar)))
  n_zero_modvar_batch <- sum(unlist(lapply(result_lst, function(x) x$zero_modvar_batch)))
  n_moderr <- sum(unlist(lapply(result_lst, function(x) x$moderr)))
  
  cat(sprintf("Found %s features with zero model variance; 
              these features won't be adjusted for batch effects.\n",
              n_zero_modvar))
  cat(sprintf("Errors encountered in %s features with model fitting; 
              these features won't be adjusted for batch effects.\n",
              n_moderr))
  if (!mean.only) {
    cat(sprintf("Found %s features with zero model variance within at least one batch; 
                these features won't be adjusted for batch effects.\n",
                n_zero_modvar_batch))
  }
  
  # convert NULLs to NAs
  gamma_hat_lst[sapply(gamma_hat_lst, is.null)] <- NA
  mu_hat_lst[sapply(mu_hat_lst, is.null)] <- NA
  phi_hat_lst[sapply(phi_hat_lst, is.null)] <- NA
  delta_hat_lst[sapply(delta_hat_lst, is.null)] <- NA
  
  # reformat lists as matrices
  gamma_hat_mat <- do.call('rbind', gamma_hat_lst)
  mu_hat_mat <- do.call('rbind', mu_hat_lst)
  phi_hat_mat <- do.call('rbind', phi_hat_lst)
  delta_hat_mat <- do.call('rbind', delta_hat_lst)
  
  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########
  if (shrink) {
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_beta
    monte_carlo_res <- lapply(1:n_batch, function(ii) {
      if (ii == 1) {
        mcres <- mcint_fun(dat = bv[, batches_ind[[ii]]], 
                           mu = mu_hat_mat[, batches_ind[[ii]]], 
                           gamma = gamma_hat_mat[, ii], 
                           phi = phi_hat_mat[, batches_ind[[ii]]],
                           delta = delta_hat_mat[, ii],
                           feature.subset.n = feature.subset.n,
                           ncores = num_cores)
      } else {
        invisible(utils::capture.output(mcres <- mcint_fun(dat = bv[, batches_ind[[ii]]], 
                                                           mu = mu_hat_mat[, batches_ind[[ii]]], 
                                                           gamma = gamma_hat_mat[, ii], 
                                                           phi = phi_hat_mat[, batches_ind[[ii]]], 
                                                           delta = delta_hat_mat[, ii],
                                                           feature.subset.n = feature.subset.n,
                                                           ncores = num_cores)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0('batch', levels(batch))
    
    gamma_star_mat <- lapply(monte_carlo_res, function(res) {res$gamma_star})
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    delta_star_mat <- lapply(monte_carlo_res, function(res) {res$delta_star})
    delta_star_mat <- do.call(cbind, delta_star_mat)
    
    ## set gamma and delta equal to 0 for reference batch (probably unnecessary, 
    ## but just to make sure)
    if (!is.null(ref.batch)) {
      gamma_star_mat[, ref] <- 0
      delta_star_mat[, ref] <- 0
    }
    
    if (mean.only) {
      cat("Apply shrinkage to mean only\n")
      delta_star_mat <- delta_hat_mat
    }
  } else {
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat_mat
    delta_star_mat <- delta_hat_mat
  }
  
  ########  Obtain adjusted batch-free distribution  ########
  mu_star_mat <- matrix(NA, nrow = nrow(bv), ncol = ncol(bv))
  phi_star_mat <- phi_hat_mat
  
  for (jj in 1:n_batch) {
    logit_mu_star_subset <- log(mu_hat_mat[, batches_ind[[jj]]] / 
                                  (1 - mu_hat_mat[, batches_ind[[jj]]])) - 
      vec2mat(gamma_star_mat[, jj], n_batches[jj])
    mu_star_mat[, batches_ind[[jj]]] <- exp(logit_mu_star_subset) / 
      (1 + exp(logit_mu_star_subset))
    if (!mean.only) {
      log_phi_star_subset <- log(phi_hat_mat[, batches_ind[[jj]]]) +
        vec2mat(delta_hat_mat[, jj], n_batches[jj]) -
        vec2mat(delta_star_mat[, jj], n_batches[jj])
      phi_star_mat[, batches_ind[[jj]]] <- exp(log_phi_star_subset)
    }
  }
  
  ########  Adjust the data  ########
  cat("Adjusting the data\n")
  adj_bv_raw <- matrix(NA, nrow = nrow(bv), ncol = ncol(bv))
  for (kk in 1:n_batch) {
    bv_sub <- bv[, batches_ind[[kk]]]
    old_mu <- mu_hat_mat[, batches_ind[[kk]]]
    old_phi <- phi_hat_mat[, batches_ind[[kk]]] * exp(delta_hat_mat[, kk])
    new_mu <- mu_star_mat[, batches_ind[[kk]]]
    new_phi <- phi_star_mat[, batches_ind[[kk]]]
    adj_bv_raw[, batches_ind[[kk]]] <- match_quantiles_beta(bv_sub = bv_sub,
                                                            old_mu = old_mu, 
                                                            old_phi = old_phi, 
                                                            new_mu = new_mu, 
                                                            new_phi = new_phi)
  }
  
  ## Add back features with values unqualified for beta regression  
  ## so that dimensions won't change
  if (dtype == "b-value") {
    adj_vmat <- vmatOri
    adj_vmat[keep, ] <- adj_bv_raw
  } else {
    adj_vmat <- vmatOri
    adj_vmat[keep, ] <- log(adj_bv_raw / (1 - adj_bv_raw))
  }
  return(adj_vmat)
}
