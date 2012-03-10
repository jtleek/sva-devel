num.sv <- function(dat, mod,method=c("be","leek"),vfilter=NULL,B=20, sv.sig=0.10,seed=NULL) {
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


sva <- function(dat, mod, mod0 = NULL,n.sv=NULL,method=c("irw","two-step"),vfilter=NULL,B=5, numSVmethod = "be") {
#Input
#=============================================================================
#dat: A m genes by n arrays matrix of expression data
#mod: A model matrix for the terms included in the analysis
#mod0: A null model matrix
#n.sv: The number of surrogate variables to build
#vfilter: The number of most variable genes/probes to use in the analysis (must be > 100 and < m)
#B: The number of iterations to perform
#numSVmethod: The method for determining the number of SVs (see num.sv for more details)

  
#Output
#=============================================================================
#sv: A n by n.sv matrix where each column is a distinct surrogate variable
#pprob.gam: A vector of the posterior probabilities each gene is affected by heterogeneity
#pprob.b: A vector of the posterior probabilities each gene is affected by mod
#n.sv: The number of significant surrogate variables

  if(is.null(n.sv)){
    n.sv = num.sv(dat,mod,method=numSVmethod,vfilter=vfilter)
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

fsva <- function(dbdat,mod,sv,newdat=NULL,method=c("fast","exact")){

#Input
#=============================================================================
#dbdat: A m genes by n arrays matrix of expression data from the database/training data
#mod: The model matrix for the terms included in the analysis for the training data
#sv: The surrogate variable object created by running sva on dbdat using mod. 
#newdat (optional): A set of test samples to be adjusted using the training database
#method: If method ="fast" then the SVD is calculated using an online approach, this may introduce slight bias. If method="exact" the exact SVD is calculated, but will be slower
  
#Output
#=============================================================================
#db: An adjusted version of the training database where the effect of batch/expression heterogeneity has been removed
#new: An adjusted version of the new samples, adjusted one at a time using the fsva methodology. 
#newsv: Surrogate variables for the new samples

  method <- match.arg(method)
  ndb = dim(dbdat)[2]
  nnew = dim(newdat)[2]
  nmod = dim(mod)[2]
  ntot = ndb + nnew
  n.sv<-sv$n.sv

  # Adjust the database
  mod = cbind(mod,sv$sv)
  gammahat = (dbdat %*% mod %*% solve(t(mod) %*% mod))[,(nmod+1):(nmod + sv$n.sv)]
  db = dbdat - gammahat %*% t(sv$sv)
  adjusted=NA

  if(!is.null(newdat)){
  # Calculate the new surrogate variables
    if(method=="fast"){
      wts <- (1-sv$pprob.b)*sv$pprob.gam
      WX <- wts*dbdat
      svd.wx = fast.svd(t(scale(t(WX),scale=F)),tol=0)
      D <- svd.wx$d[1:n.sv]
      U <- svd.wx$u[,1:n.sv]
      P <- t(wts*t(1/D * t(U)))
      newV <- P %*% newdat
      sgn = rep(NA,n.sv)
      for(j in 1:sv$n.sv){
        if(sv$n.sv>1){
          sgn[j] = sign(cor(svd.wx$v[1:ndb,j],sv$sv[1:ndb,j]))
        }
        if(sv$n.sv==1){
          sgn[j] = sign(cor(svd.wx$v[1:ndb,j],sv$sv[1:ndb]))
        }
      }
      newV = newV*sgn

      
    }else if(method=="exact"){
      newV = matrix(nrow=nnew,ncol=n.sv)
      for(i in 1:nnew){
        tmp = cbind(dbdat,newdat[,i])
        tmpd = (1-sv$pprob.b)*sv$pprob.gam*tmp
        ss = fast.svd(t(scale(t(tmpd),scale=F)),tol=0)
        sgn = rep(NA,sv$n.sv)
        for(j in 1:sv$n.sv){
          if(sv$n.sv>1){
            sgn[j] = sign(cor(ss$v[1:ndb,j],sv$sv[1:ndb,j]))
          }
          if(sv$n.sv==1){
            sgn[j] = sign(cor(ss$v[1:ndb,j],sv$sv[1:ndb]))
          }
        }
        newV[i,]<-ss$v[(ndb+1),1:sv$n.sv]*sgn
      }
     newV = t(newV) 
    }

    newV = t(newV)
    newV = scale(newV)/sqrt(dim(newV)[1])
    newV = t(newV)
    adjusted = newdat - gammahat %*% newV
    newV = t(newV)
  }

  return(list(db=db,new=adjusted,newsv = newV))
}



# this assumes that the data is already preprocessed, normalized, and clean
#	but I retained some functions that do this
# ---- Inputs: -----#
#	> dat = probe x sample genomic measure matrix, ie 'p' for methylation data 
#	> batch = the "batch" variable its
#		right now, its just set to be  processing month
#	> mod0 = sample x covariate information, that you want to 'protect' NOT including the batch variable
#	> par.prior = use parametric priors? taken from the original function:
#		'T' uses the parametric adjustments, if 'F' uses the nonparametric adjustments
#		--if you are unsure what to use, try the parametric adjustments (they run faster)
#		'par.prior' if  and check the plots to see if these priors are reasonable; 
#	> prior.plots = plot the priors? also taken from the original function:
#		if true will give prior plots with black as a kernal estimate 
#		of the empirical batch effect density and red as the parametric estimate. 
#        > numCovs = the column numbers of the variables in mod to be treated as continuous

ComBat <- function(dat, batch, mod, numCovs = NULL, par.prior=TRUE,prior.plots=FALSE) {

	mod = cbind(mod,batch)
	
	# check for intercept, and drop if present
	check = apply(mod, 2, function(x) all(x == 1))
	mod = as.matrix(mod[,!check])
	
	colnames(mod)[ncol(mod)] = "Batch"
	
	if(sum(check) > 0 & !is.null(numCovs)) numCovs = numCovs-1
	
	design <- design.mat(mod,numCov = numCovs)	

	batches <- list.batch(mod)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	
	## Check for missing values
	NAs = any(is.na(dat))
	if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
        #print(dat[1:2,])
	##Standardize Data across genes
	cat('Standardizing Data across genes\n')
	if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))}else{B.hat=apply(dat,1,Beta.NA,design)} #Standarization Model
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	if (!NAs){var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=T)}

	stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
	s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

	##Get regression batch effect parameters
	cat("Fitting L/S model and finding priors\n")
	batch.design <- design[,1:n.batch]
	if (!NAs){
		gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
	} else{
		gamma.hat=apply(s.data,1,Beta.NA,batch.design)
	
	}
	delta.hat <- NULL
	for (i in batches){
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
		}

	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean)
	t2 <- apply(gamma.hat, 1, var)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)

	
	##Plot empirical and parametric priors

	if (prior.plots & par.prior){
		par(mfrow=c(2,2))
		tmp <- density(gamma.hat[1,])
		plot(tmp,  type='l', main="Density Plot")
		xx <- seq(min(tmp$x), max(tmp$x), length=100)
		lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
		qqnorm(gamma.hat[1,])	
		qqline(gamma.hat[1,], col=2)	
	
		tmp <- density(delta.hat[1,])
		invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
		tmp1 <- density(invgam)
		plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
		lines(tmp1, col=2)
		qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
		lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
		title('Q-Q Plot')
	}
	
	##Find EB batch adjustments

	gamma.star <- delta.star <- NULL
	if(par.prior){
		cat("Finding parametric adjustments\n")
		for (i in 1:n.batch){
			temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],
				delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
	}else{
		cat("Finding nonparametric adjustments\n")
		for (i in 1:n.batch){
			temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
		}


	### Normalize the Data ###
	cat("Adjusting the Data\n")

	bayesdata <- s.data
	j <- 1
	for (i in batches){
		bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
		}

	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean

	return(bayesdata)

}

# filters data based on presence/absence call
filter.absent <- function(x,pct){
	present <- T
	col <- length(x)/2
	pct.absent <- (sum(x[2*(1:col)]=="A") + sum(x[2*(1:col)]=="M"))/col
	if(pct.absent > pct){present <- F}
	present
	}

# Next two functions make the design matrix (X) from the sample info file 
build.design <- function(vec, des=NULL, start=2){
	tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
	for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
	cbind(des,tmp)
	}

design.mat <- function(mod, numCov){

	tmp <- which(colnames(mod) == 'Batch')
	tmp1 <- as.factor(mod[,tmp])
	cat("Found",nlevels(tmp1),'batches\n')
	design <- build.design(tmp1,start=1)

	if(!is.null(numCov)) {
		theNumCov = as.matrix(mod[,numCov])
		mod0 = as.matrix(mod[,-c(numCov,tmp)])
	} else 	mod0 = as.matrix(mod[,-tmp])

	ncov <- ncol(mod0)
	
	cat("Found",ncov,' categorical covariate(s)\n')
	if(!is.null(numCov)) cat("Found",ncol(theNumCov),' continuous covariate(s)\n')
	if(ncov>0){
		for (j in 1:ncov){
			tmp1 <- as.factor(as.matrix(mod0)[,j])
			design <- build.design(tmp1,des=design)
			}
		}
	if(!is.null(numCov)) design = cbind(design,theNumCov)
	return(design)

}



# Makes a list with elements pointing to which array belongs to which batch
list.batch <- function(mod){
	tmp1 <- as.factor(mod[,which(colnames(mod) == 'Batch')])
	batches <- NULL
	for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
	batches
	}

# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
	tmp <- strsplit(colnames(dat),'\\.')
	tr <- NULL
	for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
	tr
	}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}


# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
	n <- apply(!is.na(sdat),1,sum)
	g.old <- g.hat
	d.old <- d.hat
	change <- 1
	count <- 0
	while(change>conv){
		g.new <- postmean(g.hat,g.bar,n,d.old,t2)
		sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
		d.new <- postvar(sum2,n,a,b)
		change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
		g.old <- g.new
		d.old <- d.new
		count <- count+1
		}
	#cat("This batch took", count, "iterations until convergence\n")
	adjust <- rbind(g.new, d.new)
	rownames(adjust) <- c("g.star","d.star")
	adjust
	}

#likelihood function used below
L <- function(x,g.hat,d.hat){prod(dnorm(x,g.hat,sqrt(d.hat)))}

# Monte Carlo integration functients
int.eprior <- function(sdat,g.hat,d.hat){
	g.star <- d.star <- NULL
	r <- nrow(sdat)
	for(i in 1:r){
		g <- g.hat[-i]
		d <- d.hat[-i]		
		x <- sdat[i,!is.na(sdat[i,])]
		n <- length(x)
		j <- numeric(n)+1
		dat <- matrix(as.numeric(x),length(g),n,byrow=T)
		resid2 <- (dat-g)^2
		sum2 <- resid2%*%j
		LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
		LH[LH=="NaN"]=0
		g.star <- c(g.star,sum(g*LH)/sum(LH))
		d.star <- c(d.star,sum(d*LH)/sum(LH))
		#if(i%%1000==0){cat(i,'\n')}
		}
	adjust <- rbind(g.star,d.star)
	rownames(adjust) <- c("g.star","d.star")
	adjust	
	} 

#fits the L/S model in the presence of missing data values

Beta.NA = function(y,X){
	des=X[!is.na(y),]
	y1=y[!is.na(y)]
	B <- solve(t(des)%*%des)%*%t(des)%*%y1
	B
	}


