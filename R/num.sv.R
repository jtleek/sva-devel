#' A function for calculating the number of surrogate variables to estimate in a model 
#' 
#' This function estimates the number of surrogate variables that should be included
#' in a differential expression model. The default approach is based on a permutation
#' procedure originally prooposed by Buja and Eyuboglu 1992. The function also provides
#' an interface to the asymptotic approach proposed by Leek 2011 Biometrics.
#' 
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param mod The model matrix being used to fit the data
#' @param method One of "be" or "leek" as described in the details section
#' @param vfilter You may choose to filter to the vfilter most variable rows before performing the analysis
#' @param B The number of permutaitons to use if method = "be"
#' @param seed Set a seed when using the permutation approach 
#' 
#' @return n.sv The number of surrogate variables to use in the sva software
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
#' 
#' @export
#' 
num.sv <- function(dat, mod,method=c("be","leek"),vfilter=NULL,B=20,seed=NULL) {  
  if(!is.null(vfilter)){
    if(vfilter < 100 | vfilter > dim(dat)[1]){
      stop(paste("The number of genes used in the analysis must be between 100 and",dim(dat)[1],"\n"))
    }
    tmpv = rowVars(dat)
    ind = which(rank(-tmpv) < vfilter)
    dat = dat[ind,]
  }
  
  method <- match.arg(method)
  if(method=="be"){
    if(!is.null(seed)){set.seed(seed)}
    warn <- NULL
    n <- ncol(dat)
    m <- nrow(dat)
    H <- mod %*% solve(t(mod) %*% mod) %*% t(mod) 
    res <- dat - t(H %*% t(dat))
    uu <- svd(res)
    ndf <- min(m,n) - ceiling(sum(diag(H)))
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    dstat0 <- matrix(0,nrow=B,ncol=ndf)
    
    for(i in 1:B){
      res0 <- t(apply(res, 1, sample, replace=FALSE))
      res0 <- res0 - t(H %*% t(res0))
      uu0 <- svd(res0)
      dstat0[i,] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
    psv <- rep(1,n)
    for(i in 1:ndf){
      psv[i] <- mean(dstat0[,i] >= dstat[i])
    }
    for(i in 2:ndf){
      psv[i] <- max(psv[(i-1)],psv[i]) 
    }
    
    nsv <- sum(psv <= 0.10)
    return(as.numeric(list(n.sv = nsv)))
  }else{
    dat <- as.matrix(dat)
    dims <- dim(dat)
    a <- seq(0,2,length=100)
    n <- floor(dims[1]/10)
    rhat <- matrix(0,nrow=100,ncol=10)
    P <- (diag(dims[2])-mod %*% solve(t(mod) %*% mod) %*% t(mod))  
    for(j in 1:10){
      dats <- dat[1:(j*n),]
      ee <- eigen(t(dats) %*% dats)
      sigbar <- ee$values[dims[2]]/(j*n)
      R <- dats %*% P
      wm <- (1/(j*n))*t(R) %*% R - P*sigbar
      ee <- eigen(wm)
      v <- c(rep(TRUE, 100), rep(FALSE, dims[2]))
      v <- v[order(c(a*(j*n)^(-1/3)*dims[2],ee$values), decreasing = TRUE)]
      u <- 1:length(v)
      w <- 1:100
      rhat[,j] <- rev((u[v==TRUE]-w))
    }
    ss <- rowVars(rhat)
    
    bumpstart <- which.max(ss > (2*ss[1]))
    start <- which.max(c(rep(1e5,bumpstart),ss[(bumpstart+1):100]) < 0.5*ss[1])
    finish <- which.max(ss*c(rep(0,start),rep(1,100-start)) > ss[1])
    if(finish==1){finish <- 100}
    
    n.sv <- modefunc(rhat[start:finish,10])
    return(n.sv)
    print(method)
  }
}

