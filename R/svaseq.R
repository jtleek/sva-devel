#' A function for estimating surrogate variables for count based RNA-seq data.
#'
#' This function is the implementation of the iteratively re-weighted least squares
#' approach for estimating surrogate variables. As a by product, this function
#' produces estimates of the probability of being an empirical control. This function first
#' applies a moderated log transform as described in Leek 2014 before calculating the surrogate
#' variables. See the function \code{\link{empirical.controls}} for a direct estimate of the empirical controls.
#'
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param mod The model matrix being used to fit the data
#' @param mod0 The null model being compared when fitting the data
#' @param n.sv The number of surogate variables to estimate
#' @param controls A vector of probabilities (between 0 and 1, inclusive) that each gene is a control. A value of 1 means the gene is certainly a control and a value of 0 means the gene is certainly not a control.
#' @param method For empirical estimation of control probes use "irw". If control probes are known use "supervised"
#' @param vfilter You may choose to filter to the vfilter most variable rows before performing the analysis. vfilter must be NULL if method is "supervised"
#' @param B The number of iterations of the irwsva algorithm to perform
#' @param numSVmethod If n.sv is NULL, sva will attempt to estimate the number of needed surrogate variables. This should not be adapted by the user unless they are an expert.
#' @param constant The function takes log(dat + constant) before performing sva. By default constant = 1, all values of dat + constant should be positive.
#'
#' @return sv The estimated surrogate variables, one in each column
#' @return pprob.gam: A vector of the posterior probabilities each gene is affected by heterogeneity
#' @return pprob.b A vector of the posterior probabilities each gene is affected by mod
#' @return n.sv The number of significant surrogate variables
#'
#' @examples
#' library(zebrafishRNASeq)
#' data(zfGenes)
#' filter = apply(zfGenes, 1, function(x) length(x[x>5])>=2)
#' filtered = zfGenes[filter,]
#' genes = rownames(filtered)[grep("^ENS", rownames(filtered))]
#' controls = grepl("^ERCC", rownames(filtered))
#' group = as.factor(rep(c("Ctl", "Trt"), each=3))
#' dat0 = as.matrix(filtered)
#'
#' mod1 = model.matrix(~group)
#' mod0 = cbind(mod1[,1])
#' svseq = svaseq(dat0,mod1,mod0,n.sv=1)$sv
#' plot(svseq,pch=19,col="blue")
#'
#' @export
#'
svaseq <- function(dat, mod, mod0 = NULL,n.sv=NULL,controls=NULL,method=c("irw","two-step","supervised"),
                vfilter=NULL,B=5, numSVmethod = "be",constant = 1) {
  method <- match.arg(method)
  if(!is.null(controls) & !is.null(vfilter)){stop("sva error: if controls is provided vfilter must be NULL.\n")}
  if((method=="supervised") & is.null(controls)){stop("sva error: for a supervised analysis you must provide a vector of controls.\n")}
  if(!is.null(controls) & (method!="supervised")){method = "supervised"; cat("sva warning: controls provided so supervised sva is being performed.\n")}

  if(any(dat < 0)){stop("svaseq error: counts must be zero or greater")}
  dat = log(dat + constant)


  if(!is.null(vfilter)){
    if(vfilter < 100 | vfilter > dim(dat)[1]){
      stop(paste("sva error: the number of genes used in the analysis must be between 100 and",dim(dat)[1],"\n"))
    }
    tmpv = rowVars(dat)
    ind = which(rank(-tmpv) < vfilter)
    dat = dat[ind,]
  }

  if(is.null(n.sv)){
    n.sv = num.sv(dat=dat,mod=mod,method=numSVmethod,vfilter=vfilter)
  }

  if(n.sv > 0){
    cat(paste("Number of significant surrogate variables is: ",n.sv,"\n"))

    if(method=="two-step"){
      return(twostepsva.build(dat=dat, mod=mod,n.sv=n.sv))
    }
    if(method=="irw"){
      return(irwsva.build(dat=dat, mod=mod, mod0 = mod0,n.sv=n.sv,B=B))
    }
    if(method=="supervised"){
      return(ssva(dat,controls,n.sv=n.sv))
    }
  }else{
    cat("No significant surrogate variables\n"); return(list(sv=0,pprob.gam=0,pprob.b=0,n.sv=0))
  }

}
