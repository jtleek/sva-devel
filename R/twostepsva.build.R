#' A function for estimating surrogate variables with the two step approach of Leek and Storey 2007
#' 
#' This function is the implementation of the two step approach for estimating surrogate
#' variables proposed by Leek and Storey 2007 PLoS Genetics. This function is primarily
#' included for backwards compatibility. Newer versions of the sva algorithm are available
#' through \code{\link{sva}}, \code{\link{svaseq}}, with low level functionality available
#' through \code{\link{irwsva.build}} and \code{\link{ssva}}.
#' 
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param mod The model matrix being used to fit the data
#' @param n.sv The number of surogate variables to estimate
#'
#' @return sv The estimated surrogate variables, one in each column
#' @return pprob.gam: A vector of the posterior probabilities each gene is affected by heterogeneity
#' @return pprob.b A vector of the posterior probabilities each gene is affected by mod (this is always null for the two-step approach)
#' @return n.sv The number of significant surrogate variables
#' 
#' @examples 
#' library(bladderbatch)
#' library(limma)
#' data(bladderdata)
#' dat <- bladderEset
#' 
#' pheno = pData(dat)
#' edata = exprs(dat)
#' mod = model.matrix(~as.factor(cancer), data=pheno)
#' 
#' n.sv = num.sv(edata,mod,method="leek")
#' svatwostep <- twostepsva.build(edata,mod,n.sv)
#' 
#' @export
#' 


twostepsva.build <- function(dat, mod, n.sv){

  n <- ncol(dat)
  m <- nrow(dat)
  H <- mod %*% solve(t(mod) %*% mod) %*% t(mod) 
  res <- dat - t(H %*% t(dat))
  uu <- svd(res)
  ndf <- n - ceiling(sum(diag(H)))
  dstat <-  uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  res.sv <- uu$v[,1:n.sv, drop=FALSE]
  
  use.var <- matrix(rep(FALSE, n.sv*m), ncol=n.sv)
  pp <- matrix(rep(FALSE, n.sv*m), ncol=n.sv)
  
  
  for(i in 1:n.sv) {
    mod <- cbind(rep(1,n),res.sv[,i])
    mod0 <- cbind(rep(1,n))
    pp[,i] <-f.pvalue(dat,mod,mod0)
    use.var[,i] <- edge.lfdr(pp[,i]) < 0.10
  }
  
  for(i in ncol(use.var):1) {
    if(sum(use.var[,i]) < n) {
      use.var <- as.matrix(use.var[,-i])
      n.sv <- n.sv - 1
      if(n.sv <= 0){break}
    }
  }
  if(n.sv >0){
    sv <- matrix(0,nrow=n,ncol=n.sv)
    dat <- t(scale(t(dat),scale=FALSE))
    for(i in 1:n.sv) {
      uu <- svd(dat[use.var[,i],])
      maxcor <- 0
      for(j in 1:(n-1)) {
        if(abs(cor(uu$v[,j], res.sv[,i])) > maxcor) {
          maxcor <- abs(cor(uu$v[,j], res.sv[,i]))
          sv[,i] <- uu$v[,j]
        }
      }
    }
    pprob.gam <- use.var %*% rep(1,n.sv) > 0
    retval <- list(sv=sv,pprob.gam=pprob.gam,pprob.b=NULL,n.sv=n.sv)
    return(retval)
  }
  else{
    sv <- rep(0,n)
    ind <- rep(0,m)
    n.sv <- 0
    retval <- list(sv=sv,pprob.gam=rep(0,m),pprob.b = NULL,n.sv=n.sv)
    return(retval)
  }
}
