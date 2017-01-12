#' A function for estimating surrogate variables by estimating empirical control probes
#' 
#' This function is the implementation of the iteratively re-weighted least squares
#' approach for estimating surrogate variables. As a buy product, this function
#' produces estimates of the probability of being an empirical control. See the function
#' \code{\link{empirical.controls}} for a direct estimate of the empirical controls. 
#' 
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param mod The model matrix being used to fit the data
#' @param mod0 The null model being compared when fitting the data
#' @param n.sv The number of surogate variables to estimate
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
#' 
#' n.sv = num.sv(edata,mod,method="leek")
#' res <- irwsva.build(edata, mod, mod0 = NULL,n.sv) 
#' 
#' @export
#' 

irwsva.build <- function(dat, mod, mod0 = NULL,n.sv) {

  n <- ncol(dat)
  m <- nrow(dat)
  if(is.null(mod0)){mod0 <- mod[,1]}
  Id <- diag(n)
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))  
  uu <- eigen(t(resid)%*%resid)
  vv <- uu$vectors
  ndf <- n - dim(mod)[2]
  
  pprob <- rep(1,m)
  one <- rep(1,n)
  Id <- diag(n)
  df1 <- dim(mod)[2] + n.sv
  df0 <- dim(mod0)[2]  + n.sv
  dats2 <- dat
  
  rm(resid)

  while(num.sv(dats2,mod,method="be") != 0){
    mod.b <- cbind(mod,uu$vectors[,1:n.sv])
    mod0.b <- cbind(mod0,uu$vectors[,1:n.sv])
    ptmp <- f.pvalue(dat,mod.b,mod0.b)
    pprob.b <- (1-edge.lfdr(ptmp))
    
    mod.gam <- cbind(mod0,uu$vectors[,1:n.sv])
    mod0.gam <- cbind(mod0)
    ptmp <- f.pvalue(dat,mod.gam,mod0.gam)
    pprob.gam <- (1-edge.lfdr(ptmp))
    pprob <- pprob.gam*(1-pprob.b)
    dats <- dat*pprob
    dats <- dats - rowMeans(dats)
    dats2 <- dat*pprob.b*(1-pprob.gam)
    dats2 <- dats2 - rowMeans(dats2)
    uu <- eigen(t(dats)%*%dats)
    cat(paste(i," "))
  }
  
  sv = svd(dats)$v[,1:n.sv]
  retval <- list(sv=sv,pprob.gam = pprob.gam, pprob.b=pprob.b,n.sv=n.sv)
  return(retval)
  
}
