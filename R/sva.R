num.sv <- function(dat, mod,method=c("leek","be"),vfilter=NULL,B=20, sv.sig=0.10,seed=NULL) {
#Input
#=============================================================================
#dat: A m genes by n arrays matrix of expression data
#mod: A model matrix for the terms included in the analysis
#method: Leek is the method proposed in Leek 2011 Biometrics. be is the method
#         proposed in Buja and 1992
#B: The number of null iterations to perform
#sv.sig: The significance cutoff for the surrogate variables
#Output
#=============================================================================
#n.sv: Number of significant surrogate variables. 
#id: An indicator of the significant surrogate variables
#B: number of null iterations
#sv.sig: significance level for surrogate variables
  
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
    uu <- fast.svd(res,tol=0)
    ndf <- n - ceiling(sum(diag(H)))
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    dstat0 <- matrix(0,nrow=B,ncol=ndf)
    
    for(i in 1:B){
      res0 <- t(apply(res, 1, sample, replace=FALSE))
      res0 <- res0 - t(H %*% t(res0))
      uu0 <- fast.svd(res0,tol=0)
      dstat0[i,] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
    psv <- rep(1,n)
    for(i in 1:ndf){
      psv[i] <- mean(dstat0[,i] >= dstat[i])
    }
    for(i in 2:ndf){
      psv[i] <- max(psv[(i-1)],psv[i]) 
    }
    
    nsv <- sum(psv <= sv.sig)
    return(list(n.sv = nsv))
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
      v <- c(rep(T, 100), rep(F, dims[2]))
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
  


irwsva.build <- function(dat, mod, mod0 = NULL,n.sv,B=5) {
#Input
#=============================================================================
#dat: A m genes by n arrays matrix of expression data
#mod: A model matrix for the terms included in the analysis
#mod0: A null model matrix
#n.sv: The number of surrogate variables to build
#B: The number of iterations to perform

  
#Output
#=============================================================================
#sv: A n by n.sv matrix where each column is a distinct surrogate variable
#pprob.gam: A vector of the posterior probabilities each gene is affected by heterogeneity
#pprob.b: A vector of the posterior probabilities each gene is affected by mod
#n.sv: The number of significant surrogate variables
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
  
  rm(resid)
  
  cat(paste("Iteration (out of", B,"):"))
  for(i in 1:B){
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
    uu <- eigen(t(dats)%*%dats)
    cat(paste(i," "))
  }

 sv = fast.svd(dats,tol=0)$v[,1:n.sv]
 retval <- list(sv=sv,pprob.gam = pprob.gam, pprob.b=pprob.b,n.sv=n.sv)
  return(retval)
  
}

twostepsva.build <- function(dat, mod, n.sv){
#Input
#=============================================================================
#dat: A m genes by n arrays matrix of expression data
#mod: A model matrix for the terms included in the analysis 
#n.sv: The number of surrogate variables to build
#Output
#=============================================================================
#sv: A n by n.sv matrix where each column is a distinct surrogate variable
#ind: An m by n.sv indicator of association between genes and surrogate variables
#v.sv: The percent of variation in the residual matrix explained by each eigengene
#n.sv: The number of significant surrogate variables
  n <- ncol(dat)
  m <- nrow(dat)
  H <- mod %*% solve(t(mod) %*% mod) %*% t(mod) 
  res <- dat - t(H %*% t(dat))
  uu <- fast.svd(res,tol=0)
  ndf <- n - ceiling(sum(diag(H)))
  dstat <-  uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  res.sv <- as.matrix(uu$v[,1:n.sv])
  
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
    dat <- t(scale(t(dat),scale=F))
    for(i in 1:n.sv) {
      uu <- fast.svd(dat[use.var[,i],],tol=0)
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


sva <- function(dat, mod, mod0 = NULL,n.sv=NULL,method=c("irw","two-step"),vfilter=NULL,B=5) {
#Input
#=============================================================================
#dat: A m genes by n arrays matrix of expression data
#mod: A model matrix for the terms included in the analysis
#mod0: A null model matrix
#n.sv: The number of surrogate variables to build
#vfilter: The number of most variable genes/probes to use in the analysis (must be > 100 and < m)
#B: The number of iterations to perform

  
#Output
#=============================================================================
#sv: A n by n.sv matrix where each column is a distinct surrogate variable
#pprob.gam: A vector of the posterior probabilities each gene is affected by heterogeneity
#pprob.b: A vector of the posterior probabilities each gene is affected by mod
#n.sv: The number of significant surrogate variables

  if(is.null(n.sv)){
    n.sv = num.sv(dat,mod,method="leek",vfilter=vfilter)
  }

  if(!is.null(vfilter)){
    if(vfilter < 100 | vfilter > dim(dat)[1]){
      stop(paste("The number of genes used in the analysis must be between 100 and",dim(dat)[1],"\n"))
    }
    tmpv = rowVars(dat)
    ind = which(rank(-tmpv) < vfilter)
    dat = dat[ind,]
  }
  
  if(n.sv > 0){
    cat(paste("Number of significant surrogate variables is: ",n.sv,"\n"))
    method <- match.arg(method)
    if(method=="two-step"){
      return(twostepsva.build(dat=dat, mod=mod,n.sv=n.sv))
    }
    if(method=="irw"){
      return(irwsva.build(dat=dat, mod=mod, mod0 = mod0,n.sv=n.sv,B=B))
    }
  }else{print("No significant surrogate variables"); return(list(sv=0,pprob.gam=0,pprob.b=0,n.sv=0))}
}



fstats <- function(dat,mod,mod0){
  # A function for calculating F-statistics
  # on the rows of dat, comparing the models
  # mod (alternative) and mod0 (null). 
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0,m)
  Id <- diag(n)
  
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))
  
  rss1 <- (resid*resid) %*% rep(1,n)
  rss0 <- (resid0*resid0) %*% rep(1,n)

  fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
  return(fstats)
}





## Borrowed from genefilter
rowVars <- function(x, ...) 
{
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x, ...)), ...)/(n - 1))
}

modefunc <- function(x) {
  return(as.numeric(names(sort(-table(x)))[1]))
}


f.pvalue <- function(dat,mod,mod0){
  # This is a function for performing
  # parametric f-tests on the data matrix
  # dat comparing the null model mod0
  # to the alternative model mod. 
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0,m)
  Id <- diag(n)
  
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  rss1 <- rowSums(resid*resid)
  rm(resid)
  
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))
  rss0 <- rowSums(resid0*resid0)
  rm(resid0)
  
  fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
  p <-  1-pf(fstats,df1=(df1-df0),df2=(n-df1))
  return(p)
}

mono <- function(lfdr){.Call("monotone",as.numeric(lfdr),PACKAGE="sva")}



edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
	pi0 <- mean(p >= lambda)/(1 - lambda)
        pi0 <- min(pi0, 1)

	n = length(p)
	transf = match.arg(transf)
	
	if(transf=="probit") {
		p = pmax(p, eps)
		p = pmin(p, 1-eps)
		x = qnorm(p)
		myd = density(x, adjust=adj)
		mys = smooth.spline(x=myd$x, y=myd$y)
		y = predict(mys, x)$y
		lfdr = pi0*dnorm(x)/y
	}
	
	if(transf=="logit") {
		x = log((p+eps)/(1-p+eps))
		myd = density(x, adjust=adj)
		mys = smooth.spline(x=myd$x, y=myd$y)
		y = predict(mys, x)$y
		dx = exp(x)/(1+exp(x))^2
		lfdr = pi0*dx/y
	}
	
	if(trunc) {lfdr[lfdr > 1] = 1}
	if(monotone) {	
		lfdr = lfdr[order(p)]
                lfdr = mono(lfdr)
		lfdr = lfdr[rank(p)]
	}
	
	return(lfdr)
}

fsva <- function(dbdat,mod,sv,newdat=NULL){

#Input
#=============================================================================
#dbdat: A m genes by n arrays matrix of expression data from the database/training data
#mod: The model matrix for the terms included in the analysis for the training data
#sv: The surrogate variable object created by running sva on dbdat using mod. 
#newdat (optional): A set of test samples to be adjusted using the training database

  
#Output
#=============================================================================
#db: An adjusted version of the training database where the effect of batch/expression heterogeneity has been removed
#new: An adjusted version of the new samples, adjusted one at a time using the fsva methodology. 

  sv$sv = cbind(sv$sv)
  ndb = dim(dbdat)[2]
  nnew = dim(newdat)[2]
  nmod = dim(mod)[2]
  ntot = ndb + nnew
  mod0 = mod
  
  mod = cbind(mod,sv$sv)
  gammahat = (dbdat %*% mod %*% solve(t(mod) %*% mod))[,(nmod+1):(nmod + sv$n.sv)]
  db = dbdat - gammahat %*% t(sv$sv)
  mod = mod0
  adjusted=NA
  cat(paste("Correcting sample (out of", nnew,"):"))
  if(!is.null(newdat)){
    adjusted = newdat*0
    for(i in 1:nnew){
      tmp = cbind(dbdat,newdat[,i])
      tmpd = (1-sv$pprob.b)*sv$pprob.gam*tmp
      ss = fast.svd(t(scale(t(tmpd),scale=F)),tol=0)
      mod = cbind(mod,sv$sv)
      sgn = rep(NA,sv$n.sv)
      for(j in 1:sv$n.sv){
        sgn[j] = sign(cor(ss$v[1:ndb,j],sv$sv[1:ndb,j]))
      }
      gammahat2 = t(t(gammahat)*sgn)
      
      adjusted[,i] = newdat[,i] - gammahat2 %*% ss$v[(ndb+1),1:sv$n.sv]
      mod = mod0
      cat(paste(i," "))
    }
  }
  return(list(db=db,new=adjusted))
}




