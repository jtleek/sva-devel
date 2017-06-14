#' A function for estimating surrogate variables by estimating empirical control probes
#' 
#' This function is the implementation of the iteratively re-weighted least squares
#' approach for estimating surrogate variables. As a by product, this function
#' produces estimates of the probability of being an empirical control. See the function
#' \code{\link{empirical.controls}} for a direct estimate of the empirical controls. 
#' 
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param mod The model matrix being used to fit the data
#' @param mod0 The null model being compared when fitting the data
#' @param n.sv The number of surogate variables to estimate
#' @param controls A vector of probabilities (between 0 and 1, inclusive) that each gene is a control. A value of 1 means the gene is certainly a control and a value of 0 means the gene is certainly not a control.
#' @param method For empirical estimation of control probes use "irw". If control probes are known use "supervised"
#' @param vfilter You may choose to filter to the vfilter most variable rows before performing the analysis. vfilter must be NULL if method is "supervised"
#' @param numSVmethod If n.sv is NULL, sva will attempt to estimate the number of needed surrogate variables. This should not be adapted by the user unless they are an expert. 
#' 
#' @return sv The estimated surrogate variables, one in each column
#' @return pprob.gam: A vector of the posterior probabilities each gene is affected by heterogeneity
#' @return pprob.b A vector of the posterior probabilities each gene is affected by mod
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
#' mod0 = model.matrix(~1,data=pheno)
#' 
#' n.sv = num.sv(edata,mod,method="leek")
#' svobj = sva(edata,mod,mod0,n.sv=n.sv)
#' 
#' @export
#' 

sva <- function(dat, mod, mod0 = NULL,n.sv=NULL,controls=NULL,method=c("irw","two-step","supervised"),
                vfilter=NULL, numSVmethod = "be") {
  method <- match.arg(method)
  if(!is.null(controls) & !is.null(vfilter)){stop("sva error: if controls is provided vfilter must be NULL.\n")}
  if((method=="supervised") & is.null(controls)){stop("sva error: for a supervised analysis you must provide a vector of controls.\n")}
  if(!is.null(controls) & (method!="supervised")){method = "supervised"; cat("sva warning: controls provided so supervised sva is being performed.\n")}
  
  if(!is.null(vfilter)){
    if(vfilter < 100 | vfilter > dim(dat)[1]){
      stop(paste("sva error: the number of genes used in the analysis must be between 100 and",dim(dat)[1],"\n"))
    }
    tmpv = rowVars(dat)
    ind = which(rank(-tmpv) <= vfilter)
    dat = dat[ind,]
  }
  if (!is.null(n.sv) && n.sv == 0) {
    warning("Returning zero surrogate variables as requested")
    return(list(sv=matrix(nrow=ncol(dat), ncol=0),
                pprob.gam = rep(0, nrow(dat)), pprob.b=NULL, n.sv=0))
  }
  if(is.null(n.sv)){
    n.sv = num.sv(dat,mod,method=numSVmethod,vfilter=vfilter)
  }
  
  if(n.sv > 0){
    cat(paste("Number of significant surrogate variables is: ",n.sv,"\n"))
    if(method=="two-step"){
      return(twostepsva.build(dat=dat, mod=mod,n.sv=n.sv))
    }
    if(method=="irw"){
      return(irwsva.build(dat=dat, mod=mod, mod0 = mod0,n.sv=n.sv))
    }
    if(method=="supervised"){
      return(ssva(dat,controls,n.sv))
    }
  }else{
    cat("No significant surrogate variables\n");
    return(list(sv=matrix(nrow=ncol(dat), ncol=0),
                pprob.gam = rep(0, nrow(dat)), pprob.b=NULL, n.sv=0))
  }
}
