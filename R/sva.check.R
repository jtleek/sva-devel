#' A function for post-hoc checking of an sva object to check for
#' degenerate cases. 
#' 
#' This function is designed to check for degenerate cases in the
#' sva fit and fix the sva object where possible. 
#' 
#' \code{\link{empirical.controls}} for a direct estimate of the empirical controls. 
#' 
#' @param svaobj The transformed data matrix with the variables in rows and samples in columns
#' @param dat The data set that was used to build the surrogate variables
#' @param mod The model matrix being used to fit the data
#' @param mod0 The null model matrix being used to fit the data
#' 
#' @return sv The estimated surrogate variables, one in each column
#' @return pprob.gam: A vector of the posterior probabilities each gene is affected by heterogeneity
#' @return pprob.b A vector of the posterior probabilities each gene is affected by mod
#' @return n.sv The number of significant surrogate variables
#' 
#' @examples 
#' library(bladderbatch)
#' data(bladderdata)
#' #dat <- bladderEset
#' dat <- bladderEset[1:5000,]
#' 
#' pheno = pData(dat)
#' edata = exprs(dat)
#' mod = model.matrix(~as.factor(cancer), data=pheno)
#' mod0 = model.matrix(~1,data=pheno)
#' 
#' n.sv = num.sv(edata,mod,method="leek")
#' svobj = sva(edata,mod,mod0,n.sv=n.sv)
#' svacheckobj = sva.check(svobj,edata,mod,mod0)
#' 
#' @export
#' 
#' 

sva.check <- function(svaobj,dat,mod,mod0){
  if(mean(svaobj$pprob.gam) > 0.95){
    message("Nearly all genes identified as batch associated, estimates of batch may be compromised. Please check your data carefully before using sva. Reverting to 2-step sva algorithm.")
    svaobj = sva(dat,mod,mod0,method="two-step")
  }
  return(svaobj)
}