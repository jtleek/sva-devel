#' A function for estimating surrogate variables using a supervised approach
#' 
#' This function implements a supervised surrogate variable analysis approach
#' where genes/probes known to be affected by artifacts but not by the biological
#' variables of interest are assumed to be known in advance. This supervised sva
#' approach can be called through the \code{\link{sva}} and \code{\link{svaseq}} functions
#' by specifying controls. 
#' 
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param controls A vector of probabilities (between 0 and 1, inclusive) that each gene is a control. A value of 1 means the gene is certainly a control and a value of 0 means the gene is certainly not a control. 
#' @param n.sv The number of surogate variables to estimate
#'
#' @return sv The estimated surrogate variables, one in each column
#' @return pprob.gam: A vector of the posterior probabilities each gene is affected by heterogeneity (exactly equal to controls for ssva)
#' @return pprob.b A vector of the posterior probabilities each gene is affected by mod (always null for ssva)
#' @return n.sv The number of significant surrogate variables
#' 
#' @examples 
#' library(bladderbatch)
#' data(bladderdata)
#' dat <- bladderEset[1:5000,]
#' 
#' pheno = pData(dat)
#' edata = exprs(dat)
#' mod = model.matrix(~as.factor(cancer), data=pheno)
#' 
#' n.sv = num.sv(edata,mod,method="leek")
#' set.seed(1234)
#' controls <- runif(nrow(edata))
#' ssva_res <- ssva(edata,controls,n.sv)
#' 
#' @export
#' 
#' 

ssva <- function(dat,controls,n.sv){
  if(is.null(n.sv)){stop("ssva error: You must specify the number of surrogate variables")}
  if(dim(dat)[1] != length(controls)){stop("ssva error: You must specify a control vector the same length as the number of genes.")}
  if(any(controls > 1) | any(controls < 0)){stop("ssva error: Control probabilities must be between 0 and 1.")}

  if (n.sv == 0) {
    warning("Returning zero surrogate variables as requested")
    return(list(sv=matrix(nrow=ncol(dat), ncol=0),
                pprob.gam = controls, pprob.b=NULL, n.sv=0))
  }

  dats <- dat*controls
  allZero = rowMeans(dats==0) == 1
  dats = dats[!allZero,]
  ss = svd((dats - rowMeans(dats)))
  sv = ss$v[,1:n.sv, drop=FALSE]
  return(list(sv=sv,pprob.gam = controls, pprob.b=NULL,n.sv=n.sv))
}
