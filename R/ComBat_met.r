#' Adjust for batch effects using a beta regression framework in DNA 
#' methylation data
#'
#' ComBat-met fits beta regression models to the beta-values or M-values, 
#' calculates batch-free distributions, and maps the quantiles of the 
#' estimated distributions to their batch-free counterparts.
#'
#' @param vmat matrix of beta-values or M-values
#' @param dtype data type: b-value or M-value; note that the input and output 
#' have the same data type.
#' @param batch vector for batch
#' @param group optional vector for biological condition of interest
#' @param covar_mod optional model matrix representing covariates to be 
#' included in the model
#' @param full_mod Boolean variable indicating whether to include biological 
#' condition of interest in the model
#' @param shrink Boolean variable indicating whether to apply EB-shrinkage on 
#' parameter estimation
#' @param mean.only Boolean variable indicating whether to apply EB-shrinkage 
#' on the estimation of precision effects
#' @param feature.subset.n number of features to use in non-parametric EB 
#' estimation, only useful when shrink equals TRUE
#' @param pseudo_beta pseudo beta-values to be used for replacing extreme 0 
#' and 1 beta-values.
#' Value needs to be between 0 and 0.5. Only active when dtype equals b-value.
#' @param ref.batch NULL by default. If given, that batch will be selected as 
#' a reference for batch correction.
#' @param ncores number of cores to be used for parallel computing. 
#' By default, ncores is set to one.
#'
#' @return \code{ComBat_met} returns a feature x sample matrix with the same 
#' data type as input, adjusted for batch effects.
#' @references
#' Wang J (2025) ComBat-met: adjusting batch effects in DNA methylation data. 
#' NAR Genomics and Bioinformatics, 7 (2), lqaf062. doi: 10.1093/nargab/lqaf062
#' @export
#'
#' @examples
#' # Generate a random beta-value matrix
#' bv_mat <- matrix(runif(n = 400, min = 0, max = 1), nrow = 50, ncol = 8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#'
#' # Adjust for batch effects including biological conditions
#' adj_bv_mat <- ComBat_met(bv_mat, dtype = "b-value", batch = batch, 
#' group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_bv_mat <- ComBat_met(bv_mat, dtype = "b-value", batch = batch, 
#' group = group, full_mod = FALSE)
#'
#' # Generate a random M-value matrix
#' mv_mat <- matrix(rnorm(n = 400, mean = 0, sd = 1), nrow = 50, ncol = 8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#'
#' # Adjust for batch effects including biological conditions
#' adj_mv_mat <- ComBat_met(mv_mat, dtype = "M-value", batch = batch, 
#' group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_mv_mat <- ComBat_met(mv_mat, dtype = "M-value", batch = batch, 
#' group = group, full_mod = FALSE)
#'
#' # Adjust for batch effects including biological conditions (multi-threads)
#' adj_mv_mat <- ComBat_met(mv_mat, dtype = "M-value", batch = batch, 
#' group = group, full_mod = TRUE, ncores = 2)
#'

ComBat_met <- function(
        vmat, dtype = "b-value", batch, group = NULL, covar_mod = NULL,
        full_mod = TRUE, shrink = FALSE, mean.only = FALSE,
        feature.subset.n = NULL, pseudo_beta = 1e-4, ref.batch = NULL,
        ncores = 1) {
    # Preparation
    validated <- validate_inputs(vmat, dtype, batch, group, pseudo_beta)
    vmat <- validated$vmat
    batch <- validated$batch
    group <- validated$group
    
    # Handle extreme values and uniform features
    processed <- preprocess_data(vmat, dtype, batch, mean.only, pseudo_beta)
    vmatOri <- processed$vmatOri
    vmat <- processed$vmat
    keep <- processed$keep
    mean.only.vec <- processed$mean.only.vec
    
    # Prepare design matrices
    design_data <- prepare_design_matrices(
        batch, group, covar_mod, full_mod, ref.batch)
    design <- design_data$design
    batchmod <- design_data$batchmod
    n_batch <- design_data$n_batch
    batches_ind <- design_data$batches_ind
    ref <- design_data$ref
    
    # Convert to beta values if needed
    bv <- convert_to_beta_values(vmat, dtype)
    
    # Estimate parameters
    params <- estimate_parameters(
        bv, design, batchmod, n_batch, batches_ind, mean.only.vec,
        ref.batch, ref, ncores)
    
    # Apply shrinkage if requested
    adjusted_params <- apply_shrinkage(
        shrink, bv, params, batch, n_batch, batches_ind, feature.subset.n,
        ncores, mean.only, ref.batch, ref)
    
    # Adjust the data
    adj_vmat <- adjust_data(vmatOri, dtype, keep, bv, params, adjusted_params,
                            batches_ind, mean.only)
    
    return(adj_vmat)
}
