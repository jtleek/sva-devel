# Helper functions for ComBat_met
# 
# When using ComBat_met, Mvalue_ComBat, and all functions listed in this 
# script, please cite: 
# Wang J (2025) ComBat-met: adjusting batch effects in DNA methylation data. 
# NAR Genomics and Bioinformatics, 7 (2), lqaf062. doi: 10.1093/nargab/lqaf062
# Please reach out to the authors with any questions about these functions.

#' Validate function inputs
#' 
#' Checks all input arguments for correct format, type, and consistency
#' 
#' @param vmat Matrix of beta/M values
#' @param dtype Data type: beta value or M value
#' @param batch Batch covariate
#' @param group Group covariate
#' @param pseudo_beta Value of extreme beta value replacement
#' @return List of validated inputs:
#' \itemize{
#'   \item vmat: Validated matrix
#'   \item batch: Factor-converted batch covariate
#'   \item group: Factor-converted group covariate
#' }
#' @noRd
validate_inputs <- function(vmat, dtype, batch, group, pseudo_beta) {
    if (!(is.matrix(vmat) && is.numeric(vmat))) {
        if (is.data.frame(vmat) && all(vapply(vmat, is.numeric, logical(1)))) {
            vmat <- as.matrix(vmat)
            message("Input is numeric.")
        } else {
            stop(
                "vmat must be a matrix of beta-values or M-values.")
        }
    }
    
    if (dtype == "b-value" && (pseudo_beta <= 0 | pseudo_beta >= 0.5)) {
        stop("Invalid pseudo beta-values.")
    }
    
    if (ncol(vmat) != length(batch)) {
        stop("Coverage matrix and batch vector do not have matching sizes.")
    }
    
    if (!is.null(group) && ncol(vmat) != length(group)) {
        stop("Coverage matrix and group vector do not have matching sizes.")
    }
    
    if (dtype == "b-value" && 
        (sum(vmat >= 1, na.rm = TRUE) > 0 | sum(vmat <= 0, na.rm = TRUE) > 0)
    ) {
        stop("All beta values must be between 0 and 1.")
    }
    
    if (length(dim(batch)) > 1) {
        stop("ComBat-met does not allow more than one batch variable!")
    }
    
    batch <- as.factor(batch)
    if (all(table(batch) <= 1)) {
        stop("ComBat-met doesn't support only 1 sample across all batches!")
    }
    
    if (length(levels(batch)) <= 1) {
        stop("Found only one batch. No need to adjust for batch effects!")
    }
    
    if (!is.null(group)) {
        group <- as.factor(group)
    }
    
    return(list(vmat = vmat, batch = batch, group = group))
}

#' Preprocess data
#' 
#' Handles extreme values and special cases
#' 
#' @param vmat Matrix of values
#' @param dtype Data type: beta value or M value
#' @param batch Batch covariate
#' @param mean.only Logical indicating whether to apply mean-only correction
#' @param pseudo_beta Value of extreme beta value replacement
#' @return List containing:
#' \itemize{
#'   \item vmatOri: Original data matrix
#'   \item vmat: Processed data matrix
#'   \item keep: Indices of retained features
#'   \item mean.only.vec: Logical vector indicating whether to apply 
#'   mean-only correction for each feature
#' }
#' @noRd
preprocess_data <- function(vmat, dtype, batch, mean.only, pseudo_beta) {
    ## convert extreme 0 or 1 values to pseudo-beta
    if (dtype == "b-value") {
        vmat[vmat == 0] <- pseudo_beta
        vmat[vmat == 1] <- 1 - pseudo_beta    
    }
    
    ## Correct for mean batch effects only if any batch has only 1 sample
    if (any(table(batch) == 1)) {
        message(
            "At least one batch contains only 1 sample. Only mean batch 
            effects will be corrected.")
        mean.only <- TRUE
    }
    
    ## Remove features with zero variance across all batches
    zero.var.rows.lst <- lapply(levels(batch), function(b) {
        which(
            apply(
                vmat[, batch == b], 1, 
                function(x) {stats::var(x, na.rm = TRUE) == 0}))
    })
    all.zero.var.rows <- Reduce(intersect, zero.var.rows.lst)
    
    if (length(all.zero.var.rows) > 0) {
        message(sprintf(
            "Found %s features with uniform values across all batches; 
            these features will not be adjusted for batch effects.",
            length(all.zero.var.rows)))
    }
    
    keep <- setdiff(seq_len(nrow(vmat)), all.zero.var.rows)
    vmatOri <- vmat
    vmat <- vmatOri[keep, ]
    
    ## Create a vector for correction types
    mean.only.vec <- if (mean.only) 
        rep(TRUE, length(keep)) else rep(FALSE, length(keep))
    
    return(list(
        vmatOri = vmatOri, vmat = vmat, keep = keep, 
        mean.only.vec = mean.only.vec))
}

#' Prepare design matrices
#' 
#' Creates design matrices for batch, group, and biological variables
#' 
#' @param batch Batch covariate
#' @param group Group covariate
#' @param covar_mod Covariate matrix excluding batch and group
#' @param full_mod Logical indicating whether to account for the group 
#' covariate
#' @param ref.batch Name of reference batch
#' @return List containing
#' \itemize{
#'   \item design: Design matrix
#'   \item batchmod: Batch design matrix
#'   \item n_batch: Number of batches
#'   \item batches_ind: List of samples in each batch
#'   \item ref: Integer index of reference batch
#' }
#' @noRd
prepare_design_matrices <- function(
        batch, group, covar_mod, full_mod, ref.batch) {
    n_batch <- nlevels(batch)
    batches_ind <- lapply(seq_len(n_batch), function(i) {
        which(batch == levels(batch)[i])
    })
    if (full_mod & !is.null(group) & nlevels(group) > 1) {
        message("Using full model in ComBat-met.")
        mod <- stats::model.matrix(~group)
    } else {
        message("Using null model in ComBat-met.")
        mod <- stats::model.matrix(~1, data = as.data.frame(batch))
    }
    if (!is.null(covar_mod)) {
        if (is.data.frame(covar_mod)) {
            covar_mod <- do.call(
                cbind, lapply(seq_len(ncol(covar_mod)), function(i) {
                    stats::model.matrix(~covar_mod[, i])
                }))
        }
    }
    mod <- cbind(mod, covar_mod)
    batchmod <- stats::model.matrix(~-1 + batch)
    if (!is.null(ref.batch)) {
        if (!(ref.batch %in% levels(batch))) {
            stop(
                "Reference level ref. batch is not one of the levels of the 
                batch variable.")
        }
        message(sprintf("Using batch %s as the reference batch", ref.batch))
        ref <- which(levels(batch) == ref.batch)
    } else {
        ref <- NULL
    }
    
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    message(
        "Adjusting for ", ncol(design) - ncol(batchmod), 
        ' covariate(s) or covariate level(s).')
    if (qr(design)$rank < ncol(design)) {
        handle_confounded_design(design, n_batch)
    }
    
    return(list(
        design = design, batchmod = batchmod, n_batch = n_batch, 
        batches_ind = batches_ind, ref = ref))
}

#' Check for confounded designs
#' 
#' Check if the design matrix is confounded
#' 
#' @param design Design matrix
#' @param n_batch Number of batches
#' @return No return value; throws errors if confounded
#' @noRd
handle_confounded_design <- function(design, n_batch) {
    if (ncol(design) == (n_batch+1)) {
        stop(
            "The covariate is confounded with batch! Remove the covariate and 
            rerun ComBat-met.")
    }
    if (ncol(design) > (n_batch+1)) {
        if ((
            qr(design[, -c(seq_len(n_batch))])$rank < 
            ncol(design[, -c(seq_len(n_batch))]))) {
            stop(
                'The covariates are confounded! Please remove one or more of 
                the covariates so the design is not confounded.')
        } else {
            stop(
                "At least one covariate is confounded with batch! Please 
                remove confounded covariates and rerun ComBat-met.")
        }
    }
}

#' Convert M values to beta values
#' 
#' Transforms M values to beta values when needed
#' 
#' @param vmat Matrix of values
#' @param dtype Data type: beta value or M value
#' @return Matrix of beta values
#' @noRd 
convert_to_beta_values <- function(vmat, dtype) {
    ## convert M-values to beta-values if needed
    if (dtype == "b-value") {
        return(vmat)
    } else {
        return(exp(vmat) / (1 + exp(vmat)))
    }
}

#' Estimate GLM parameters
#' 
#' Fits beta regression models to estimate batch effects and other parameters
#' 
#' @param bv Matrix of beta values
#' @param design Design matrix
#' @param batchmod Batch design matrix
#' @param n_batch Number of batches
#' @param batches_ind List of samples in each batch
#' @param mean.only.vec Logical vector indicating whether to apply mean-only 
#' correction for each feature
#' @param ref.batch Name of reference batch
#' @param ref Integer index of reference batch
#' @param ncores Number of cores for parallel processing
#' @return List of processed estimation results
#' @noRd
estimate_parameters <- function(bv, design, batchmod, n_batch, batches_ind,
                                mean.only.vec, ref.batch, ref, ncores) {
    message("Fitting the GLM model")
    
    # Validate ncores
    validate_ncores(ncores)
    
    # Set up parallel processing
    cl <- setup_parallel_processing(ncores)
    on.exit(parallel::stopCluster(cl))
    
    # Run parameter estimation in parallel
    result_lst <- run_parallel_estimation(
        cl, bv, design, batchmod, n_batch, batches_ind, 
        mean.only.vec, ref.batch, ref)
    
    # Process and summarize results
    processed_results <- process_estimation_results(result_lst, n_batch)
    
    # Report any issues found
    report_estimation_issues(
        processed_results$n_zero_modvar, 
        processed_results$n_zero_modvar_batch, 
        processed_results$n_moderr,
        mean.only.vec)
    
    return(processed_results)
}

#' Validate number of cores
#' 
#' Checks if the ncores parameter is valid (positive integer)
#' 
#' @param ncores Number of cores to use for parallel processing
#' @return No return value; throws error if validation fails.
#' @noRd
validate_ncores <- function(ncores) {
    if (!is.numeric(ncores) || ncores != as.integer(ncores) || ncores <= 0) {
        stop("ncores must be a positive integer.")
    }
}

#' Set up parallel processing
#' 
#' Initializes parallel processing with specified number of cores
#' 
#' @param ncores Number of cores to use for parallel processing
#' @return A parallel cluster object
#' @noRd
setup_parallel_processing <- function(ncores) {
    num_cores <- max(1, parallel::detectCores() - 1)
    num_cores <- min(ncores, num_cores)
    parallel::makeCluster(num_cores)
}

#' Run parameter estimation in parallel
#' 
#' Executes parameter estimation across multiple cores
#' 
#' @param cl Parallel cluster object
#' @param bv Matrix of beta values
#' @param design Design matrix
#' @param batchmod Batch design matrix
#' @param n_batch Number of batches
#' @param batches_ind List of samples in each batch
#' @param mean.only.vec Logical vector to indicate whether to apply mean-only 
#' correction for each feature
#' @param ref.batch Name of reference batch
#' @param ref Integer index of reference batch
#' @return List of estimation results for each feature
#' @noRd
run_parallel_estimation <- function(cl, bv, design, batchmod, n_batch, 
                                    batches_ind,mean.only.vec, ref.batch, 
                                    ref) {
    parallel::clusterExport(cl, varlist = c("bv", "design", "batchmod", 
                                            "n_batch", "batches_ind",
                                            "mean.only.vec", "ref.batch", 
                                            "ref",
                                            "estimate_single_feature",
                                            "check_na_values",
                                            "check_model_variance",
                                            "fit_glm_model",
                                            "calculate_parameters"),
                            envir = environment())
    
    # run in parallel
    parallel::parLapply(cl, seq_len(nrow(bv)), function(k) {
        estimate_single_feature(k, bv, design, batchmod, n_batch, batches_ind,
                                mean.only.vec, ref.batch, ref)
    })
}

#' Estimate parameters for a single feature
#' 
#' Performs parameter estimation for one feature (row) in the methylation 
#' data matrix
#' 
#' @param k Row index of feature to estimate
#' @param bv Matrix of beta values
#' @param design Design matrix
#' @param batchmod Batch design matrix
#' @param n_batch Number of batches
#' @param batches_ind List of samples in each batch
#' @param mean.only.vec Logical vector to indicate whether to apply mean-only 
#' correction for each feature
#' @param ref.batch Name of reference batch
#' @param ref Integer index of reference batch
#' @return List containing estimated parameters and any issues encountered
#' @noRd
estimate_single_feature <- function(k, bv, design, batchmod, n_batch, 
                                    batches_ind, mean.only.vec, ref.batch, 
                                    ref) {
    result<- list(
        gamma_hat = NULL, mu_hat = NULL, phi_hat = NULL, delta_hat = NULL,
        zero_modvar = 0, zero_modvar_batch = 0, moderr = 0)
    
    # Check for NA values
    na_check <- check_na_values(k, design, bv)
    if (na_check$has_issue) return(na_check$result)
    
    # Check model variance
    variance_check <- check_model_variance(
        k, design, bv, batchmod, n_batch, batches_ind,
        mean.only.vec, na_check$nona)
    if (variance_check$has_issue) return(variance_check$result)
    
    # Fit GLM
    glm_fit <- fit_glm_model(
        k, bv, design, batchmod, mean.only.vec, na_check$nona)
    # if error with model fitting
    if (inherits(glm_fit, "error")) {
        result$moderr <- 1
        return(result)
    }
    
    # Calculate parameters
    calculate_parameters(
        glm_fit, n_batch, batchmod, na_check$nona, ref.batch, ref, 
        mean.only.vec[k])
}

#' Check whether the methylation data matrix contains NAs
#' 
#' Identifies and handles NA values in the methylation data matrix
#' 
#' @param k Row index of feature to check
#' @param design Design matrix
#' @param bv Matrix of beta values
#' @return List containing:
#' \itemize{
#'   \item has_issue: Logical indicating whether NA issues were found
#'   \item result: If has_issue=TRUE, contains error result
#'   \item nona: If has_issue=FALSE, contains non-NA indices
#' }
#' @noRd
check_na_values <- function(k, design, bv) {
    # mark rows with NA values
    full_mat <- cbind(design, bv[k, ])
    nona <- which(stats::complete.cases(full_mat))
    
    # check if the data are all NAs
    if (length(nona) == 0) {
        return(list(has_issue = TRUE,
                    result = list(
                        gamma_hat = NULL, mu_hat = NULL, 
                        phi_hat = NULL, delta_hat = NULL,
                        zero_modvar = 1, zero_modvar_batch = 0, moderr = 0)))
    }
    
    list(has_issue = FALSE, nona = nona)
}

#' Check model variance
#' 
#' Verifies that the model has sufficient variance for parameter estimation
#' 
#' @param k Row index of feature to check
#' @param design Design matrix
#' @param bv Matrix of beta values
#' @param batchmod Batch design matrix
#' @param n_batch Number of batches
#' @param batches_ind List of samples in each batch
#' @param mean.only.vec Logical vector to indicate whether to apply mean-only 
#' correction for each feature
#' @param nona Indices of non-NA values
#' @return List containing:
#' \itemize{
#'   \item has_issue: Logical indicating whether variance issues were found
#'   \item result: if has_issue=TRUE, contains error result
#' }
#' @noRd
check_model_variance <- function(
        k, design, bv, batchmod, n_batch, batches_ind, mean.only.vec, nona) {
    full_mat <- cbind(design, bv[k, ])
    
    # check if the model has zero model variance
    if (qr(full_mat[nona, ])$rank < ncol(full_mat)) {
        return(list(has_issue = TRUE,
                    result = list(
                        gamma_hat = NULL, mu_hat = NULL, 
                        phi_hat = NULL, delta_hat = NULL,
                        zero_modvar = 1, zero_modvar_batch = 0, moderr = 0)))
    }
    
    # if precision correction enabled, check whether the model has zero model 
    # variance within any batch
    if (!mean.only.vec[k]) {
        for (i in seq_len(length(batches_ind))) {
            if (qr(full_mat[intersect(batches_ind[[i]], nona), 
                            c(i, (n_batch+1):ncol(full_mat))])$rank < 
                ncol(full_mat) - n_batch + 1) {
                return(list(has_issue = TRUE,
                            result = list(
                                gamma_hat = NULL, mu_hat = NULL, 
                                phi_hat = NULL, delta_hat = NULL,
                                zero_modvar = 0, zero_modvar_batch = 1, 
                                moderr = 0)))
            }
        }
    }
    
    list(has_issue = FALSE)
}

#' Fit GLM for a feature
#' 
#' Fits beta regression model for a single feature
#' 
#' @param k Row index of feature to fit
#' @param bv Matrix of beta values
#' @param design Design matrix
#' @param batchmod Batch design matrix
#' @param mean.only.vec Logical vector indicating whether to apply mean-only 
#' correction for each feature
#' @param nona Indices of non-NA values
#' @return Fitted model object or error object if fitting failed
#' @noRd
fit_glm_model <- function(k, bv, design, batchmod, mean.only.vec, nona) {
    # model fit
    if (mean.only.vec[k]) {
        tryCatch({
            betareg::betareg.fit(x = design[nona, ], y = bv[k, ][nona])
        }, error = function(e) {
            e
        })
    } else {
        tryCatch({
            betareg::betareg.fit(
                x = design[nona, ], y = bv[k, ][nona], z = batchmod[nona, ])
        }, error = function(e) {
            e
        })
    }
}

#' Calculate parameters from GLM fit
#' 
#' Extracts and calculates parameters from fitted GLM
#' 
#' @param glm_f Fitted GLM object
#' @param n_batch Number of batches
#' @param batchmod Batch design matrix
#' @param nona Indices of non-NA values
#' @param ref.batch Name of reference batch
#' @param ref Integer index of reference batch
#' @param tmp.mean.only Logical indicating whether to apply mean-only 
#' correction for the selected feature
#' @return List containing estimated parameters:
#' \itemize{
#'   \item gamma_hat: Estimated additive mean batch effects within each batch
#'   \item mu_hat: Estimated average mean across batches
#'   \item phi_hat: Estimated average precision across batches
#'   \item delta_hat: Estimated precision batch effects within each batch
#'   \item zero_modvar: Flag for zero model variance
#'   \item zero_modvar_batch: Flag for zero variance within batch
#'   \item moderr: Flag for model fitting error
#' }
#' @noRd
calculate_parameters <- function(
        glm_f, n_batch, batchmod, nona, ref.batch, ref, tmp.mean.only) {
    result <- list(
        gamma_hat = NULL, mu_hat = NULL, phi_hat = NULL, delta_hat = NULL, 
        zero_modvar = 0, zero_modvar_batch = 0, moderr = 0)
    
    # compute mean and precision intercepts as batch-size-weighted average 
    # from batches
    if (!is.null(ref.batch)) {
        alpha_x <- glm_f$coefficients$mean[ref]
    } else {
        alpha_x <- glm_f$coefficients$mean[seq_len(n_batch)] %*% 
            as.matrix(colSums(batchmod[nona, ]) / length(nona))
    }
    
    if (tmp.mean.only) {
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
    result$gamma_hat <- glm_f$coefficients$mean[seq_len(n_batch)] - 
        as.numeric(alpha_x)
    result$mu_hat <- rep(NA, nrow(batchmod))
    result$mu_hat[nona] <- glm_f$fitted.values
    result$phi_hat <- as.numeric(exp(alpha_z)) * rep(1, nrow(batchmod))
    
    if (tmp.mean.only) {
        result$delta_hat <- rep(0, n_batch)
    } else {
        result$delta_hat <- glm_f$coefficients$precision - as.numeric(alpha_z)
    }
    
    result
}

#' Process estimation results
#' 
#' Converts list of results into matrices and counts issues
#' 
#' @param result_lst List of parameter estimation results for all features
#' @param n_batch Number of batches
#' @return List containing:
#' \itemize{
#'   \item gamma_hat_mat: Matrix of additive mean batch effects within each 
#'   batch
#'   \item mu_hat_mat: Matrix of estimated average mean across batches
#'   \item phi_hat_mat: Matrix of estimated average precision across batches
#'   \item delta_hat_mat: Matrix of precision batch effects within each batch
#'   \item n_zero_modvar: Count of features with zero model variance
#'   \item n_zero_modvar_batch: Count of features with zero variance within 
#'   batch
#'   \item n_moderr: Count of model fitting errors
#' }
#' @noRd
process_estimation_results <- function(result_lst, n_batch) {
    gamma_hat_lst <- lapply(result_lst, function(x) x$gamma_hat)
    mu_hat_lst <- lapply(result_lst, function(x) x$mu_hat)
    phi_hat_lst <- lapply(result_lst, function(x) x$phi_hat)
    delta_hat_lst <- lapply(result_lst, function(x) x$delta_hat)
    
    n_zero_modvar <- sum(unlist(lapply(result_lst, function(x) x$zero_modvar)))
    n_zero_modvar_batch <- sum(
        unlist(lapply(result_lst, function(x) x$zero_modvar_batch)))
    n_moderr <- sum(unlist(lapply(result_lst, function(x) x$moderr)))
    
    # convert NULLs to NAs
    gamma_hat_lst[vapply(gamma_hat_lst, is.null, logical(1))] <- NA
    mu_hat_lst[vapply(mu_hat_lst, is.null, logical(1))] <- NA
    phi_hat_lst[vapply(phi_hat_lst, is.null, logical(1))] <- NA
    delta_hat_lst[vapply(delta_hat_lst, is.null, logical(1))] <- NA
    
    # reformat lists as matrices
    list(
        gamma_hat_mat = do.call('rbind', gamma_hat_lst),
        mu_hat_mat = do.call('rbind', mu_hat_lst),
        phi_hat_mat = do.call('rbind', phi_hat_lst),
        delta_hat_mat = do.call('rbind', delta_hat_lst),
        n_zero_modvar = n_zero_modvar,
        n_zero_modvar_batch = n_zero_modvar_batch,
        n_moderr = n_moderr
    )
}

#' Report estimation issues
#' 
#' Prints summary of any issues encountered during parameter estimation
#' 
#' @param n_zero_modvar Count of features with zero model variance
#' @param n_zero_modvar_batch Count of features with zero variance within batch
#' @param n_moderr Count of model fitting errors
#' @param mean.only.vec Logical vector indicating whether to apply mean-only 
#' correction for each feature
#' @return No return value; prints messages to console
#' @noRd
report_estimation_issues <- function(
        n_zero_modvar, n_zero_modvar_batch, n_moderr, mean.only.vec) {
    message(sprintf(
        "Found %s features with zero model variance; these features won't be 
        adjusted for batch effects. 
        Issues encountered in %s features with model fitting; these features 
        won't be adjusted for batch effects.",
        n_zero_modvar, n_moderr))
    if (!all(mean.only.vec)) {
        message(sprintf(
            "Found %s features with zero model variance within at least one 
            batch; these features won't be adjusted for batch effects.",
            n_zero_modvar_batch))
    }
}

#' Applies shrinkage to batch effect estimation
#' 
#' Computes posterior parameter estimates using Monte Carlo integration when 
#' shrink=TRUE; otherwise returns original GLM estimates
#' 
#' @param shrink Logical indicating whether to apply shrinkage correction
#' @param bv Matrix of beta values
#' @param params List containing GLM estimates
#' @param batch Batch covariate
#' @param n_batch Number of batches
#' @param batches_ind List of samples in each batch
#' @param feature.subset.n Number of features to subset for faster computation
#' @param ncores Number of cores for parallel processing
#' @param mean.only Logical indicating whether to apply mean batch effects only
#' @param ref.batch Name of reference batch
#' @param ref Integer index of reference batch
#' @return List containing:
#' \itemize{
#'   \item gamma_star_mat: Matrix of shrunken additive mean batch effects 
#'   within each batch
#'   \item delta_star_mat: Matrix of shrunken precision batch effects within 
#'   each batch
#' }
#' @noRd
apply_shrinkage <- function(shrink, bv, params, batch, n_batch, batches_ind,
                            feature.subset.n, ncores, mean.only, ref.batch, 
                            ref) {
    if (shrink) {
        message("Apply shrinkage - computing posterior estimates")
        mcint_fun <- monte_carlo_int_beta
        monte_carlo_res <- lapply(seq_len(n_batch), function(ii) {
            if (ii == 1) {
                mcres <- mcint_fun(
                    dat = bv[, batches_ind[[ii]]], 
                    mu = params$mu_hat_mat[, batches_ind[[ii]]], 
                    gamma = params$gamma_hat_mat[, ii], 
                    phi = params$phi_hat_mat[, batches_ind[[ii]]],
                    delta = params$delta_hat_mat[, ii],
                    feature.subset.n = feature.subset.n, ncores = ncores)
            } else {
                invisible(
                    utils::capture.output(
                        mcres <- mcint_fun(
                            dat = bv[, batches_ind[[ii]]], 
                            mu = params$mu_hat_mat[, batches_ind[[ii]]], 
                            gamma = params$gamma_hat_mat[, ii], 
                            phi = params$phi_hat_mat[, batches_ind[[ii]]], 
                            delta = params$delta_hat_mat[, ii], 
                            feature.subset.n = feature.subset.n, 
                            ncores = ncores)))
            }
            return(mcres)
        })
        names(monte_carlo_res) <- paste0('batch', levels(batch))
        gamma_star_mat <- lapply(monte_carlo_res, function(res) res$gamma_star)
        gamma_star_mat <- do.call(cbind, gamma_star_mat)
        delta_star_mat <- lapply(monte_carlo_res, function(res) res$delta_star)
        delta_star_mat <- do.call(cbind, delta_star_mat)
        if (!is.null(ref.batch)) {
            gamma_star_mat[, ref] <- 0
            delta_star_mat[, ref] <- 0
        }
        if (mean.only) {
            message("Apply shrinkage to mean only")
            delta_star_mat <- params$delta_hat_mat
        }
    } else {
        message("Shrinkage off - using GLM estimates for parameters")
        gamma_star_mat <- params$gamma_hat_mat
        delta_star_mat <- params$delta_hat_mat
    }
    return(list(
        gamma_star_mat = gamma_star_mat, delta_star_mat = delta_star_mat))
}

#' Replicate vector to matrix
#' 
#' Expands a vector into a matrix by replicating it column-wise n_times
#' 
#' @param vec Numeric vector to replicate
#' @param n_times Number of times to replicate
#' @return Matrix with n_times columns, each identical to vec
#' @noRd
vec2mat_met <- function(vec, n_times) {
    return(matrix(rep(vec, n_times), ncol = n_times, byrow = FALSE))
}

#' Monte Carlo integration for beta distribution
#' 
#' Computes posterior parameter estimates through Monte Carlo integration of 
#' beta likelihoods
#' 
#' @param dat Matrix of beta values for one batch
#' @param mu Matrix of GLM fitted average means across batches, for the 
#' selected batch
#' @param gamma Vector of additive mean batch effects, for the selected batch
#' @param phi Matrix of average precisions across batches, for the selected 
#' batch
#' @param delta Vector of precision batch effects, for the selected batch
#' @param feature.subset.n Subset size for faster computation
#' @param ncores Number of cores for parallel processing
#' @return List containing:
#' \itemize{
#'   \item gamma_star: Posterior estimates of additive mean batch effects, 
#'   for the selected batch
#'   \item delta_star: Posterior estimates of precision batch effects, for 
#'   the selected batch
#' }
#' @noRd
monte_carlo_int_beta <- function(
        dat, mu, gamma, phi, delta, feature.subset.n, ncores) {
    cl <- parallel::makeCluster(ncores)
    MC_int <- function(i) {
        pos_res <- c(gamma.star = NA, delta.star = NA)
        m <- mu[-i, !is.na(dat[i, ])]
        p <- phi[-i, !is.na(dat[i, ])]
        x <- dat[i, !is.na(dat[i, ])]
        gamma_sub <- gamma[-i]
        delta_sub <- delta[-i]
        # take a subset of features to do integration - save time
        if (!is.null(feature.subset.n) & is.numeric(feature.subset.n) & 
            length(feature.subset.n) == 1) {
            mcint_ind <- sample(
                seq_len(nrow(dat) - 1), feature.subset.n, replace = FALSE)
            m <- m[mcint_ind, ]
            p <- p[mcint_ind, ]
            gamma_sub <- gamma_sub[mcint_ind]
            delta_sub <- delta_sub[mcint_ind]
            G_sub <- feature.subset.n
        } else {
            G_sub <- nrow(dat) - 1
        }
        LH <- vapply(seq_len(G_sub), function(j) {
            prod(stats::dbeta(
                x, shape1 = m[j, ] * (p[j, ] * exp(delta_sub[j])), 
                shape2 = (1 - m[j, ]) * (p[j, ] * exp(delta_sub[j]))))
        }, numeric(1))
        LH[is.na(LH)] <- 0
        if (sum(LH) == 0 | is.na(sum(LH))){
            pos_res <- c(
                gamma.star = as.numeric(gamma[i]),
                delta.star = as.numeric(delta[i]))
        } else {
            pos_res["gamma.star"] <- sum(gamma_sub * LH, na.rm = TRUE) / 
                sum(LH[!is.na(gamma_sub)])
            pos_res["delta.star"] <- sum(delta_sub * LH, na.rm = TRUE) / 
                sum(LH[!is.na(delta_sub)])
        }
        return(pos_res)
    }
    # run in parallel
    pos_res_lst <- parallel::parLapply(cl, seq_len(nrow(dat)), MC_int)
    parallel::stopCluster(cl)
    pos_res_mat <- do.call(rbind, pos_res_lst)
    res <- list(gamma_star = pos_res_mat[, "gamma.star"], 
                delta_star = pos_res_mat[, "delta.star"]) 
    return(res)
}

#' Adjust data for batch effects
#' 
#' Removes batch effects from beta/M values
#' 
#' @param vmatOri Original data matrix
#' @param dtype Data type: beta value or M value
#' @param keep Indices of features being adjusted
#' @param bv Matrix of beta values
#' @param params List of GLM estimates
#' @param adjusted_params List of shrunken estimates
#' @param batches_ind List of samples in each batch
#' @param mean.only Logical indicating whether to apply mean-only correction
#' @return Batch-adjusted matrix matching input dimensions
#' @noRd
adjust_data <- function(vmatOri, dtype, keep, bv, params, adjusted_params, 
                        batches_ind, mean.only) {
    gamma_star_mat <- adjusted_params$gamma_star_mat
    delta_star_mat <- adjusted_params$delta_star_mat
    mu_hat_mat <- params$mu_hat_mat
    phi_hat_mat <- params$phi_hat_mat
    delta_hat_mat <- params$delta_hat_mat
    n_batch <- length(batches_ind)
    n_batches <- vapply(batches_ind, length, integer(1))
    mu_star_mat <- matrix(NA, nrow = nrow(bv), ncol = ncol(bv))
    phi_star_mat <- phi_hat_mat
    for (jj in seq_len(n_batch)) {
        logit_mu_star_subset <- log(
            mu_hat_mat[, batches_ind[[jj]]] / 
                (1 - mu_hat_mat[, batches_ind[[jj]]])) - 
            vec2mat_met(gamma_star_mat[, jj], n_batches[jj])
        mu_star_mat[, batches_ind[[jj]]] <- exp(logit_mu_star_subset) / 
            (1 + exp(logit_mu_star_subset))
        if (!mean.only) {
            log_phi_star_subset <- log(phi_hat_mat[, batches_ind[[jj]]]) +
                vec2mat_met(delta_hat_mat[, jj], n_batches[jj]) -
                vec2mat_met(delta_star_mat[, jj], n_batches[jj])
            phi_star_mat[, batches_ind[[jj]]] <- exp(log_phi_star_subset)
        }
    }
    message("Adjusting the data")
    adj_bv_raw <- matrix(NA, nrow = nrow(bv), ncol = ncol(bv))
    for (kk in seq_len(n_batch)) {
        bv_sub <- bv[, batches_ind[[kk]]]
        old_mu <- mu_hat_mat[, batches_ind[[kk]]]
        old_phi <- phi_hat_mat[, batches_ind[[kk]]] * exp(delta_hat_mat[, kk])
        new_mu <- mu_star_mat[, batches_ind[[kk]]]
        new_phi <- phi_star_mat[, batches_ind[[kk]]]
        adj_bv_raw[, batches_ind[[kk]]] <- match_quantiles_beta(
            bv_sub = bv_sub,
            old_mu = old_mu, 
            old_phi = old_phi, 
            new_mu = new_mu, 
            new_phi = new_phi)
    }
    if (dtype == "b-value") {
        adj_vmat <- vmatOri
        adj_vmat[keep, ] <- adj_bv_raw
    } else {
        adj_vmat <- vmatOri
        adj_vmat[keep, ] <- log(adj_bv_raw / (1 - adj_bv_raw))
    }
    return(adj_vmat)
}

#' Match quantiles for beta-distributed data
#' 
#' Adjusts beta values to match target distribution quantiles
#' 
#' @param bv_sub Matrix of beta values for one batch
#' @param old_mu Original means
#' @param old_phi Original precisions
#' @param new_mu Target means
#' @param new_phi Target precisions
#' @return Matrix of adjusted beta values
#' @noRd
match_quantiles_beta <- function(bv_sub, old_mu, old_phi, new_mu, new_phi) {
    new_bv_sub <- matrix(NA, nrow = nrow(bv_sub), ncol = ncol(bv_sub))
    for (a in seq_len(nrow(bv_sub))) {
        for (b in seq_len(ncol(bv_sub))) {
            tmp_p <- stats::pbeta(
                bv_sub[a, b], shape1 = old_mu[a, b] * old_phi[a, b], 
                shape2 = (1 - old_mu[a, b]) * old_phi[a, b])
            if (is.na(tmp_p)) {
                new_bv_sub[a, b] <- bv_sub[a, b]  
            } else {
                new_bv_sub[a, b] <- stats::qbeta(
                    tmp_p, shape1 = new_mu[a, b] * new_phi[a, b], 
                    shape2 = (1 - new_mu[a, b]) * new_phi[a, b])
            }
        }
    }
    return(new_bv_sub)
}
