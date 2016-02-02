#' Adjust for batch effects using an empirical Bayes framework
#' 
#' ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology 
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for 
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects. The input
#' data are assumed to be cleaned and normalized before batch effect removal. 
#' 
#' @param dat Genomic measure matrix (dimensions probe x sample) - for example, expression matrix
#' @param batch {Batch covariate (only one batch allowed)}
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the mean of the batch effect (no scale adjustment)
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a reference for batch adjustment. 
#'
#' @return data A probe x sample genomic measure matrix, adjusted for batch effects.
#' 
#' @examples 
#' library(bladderbatch)
#' data(bladderdata)
#' dat <- bladderEset[1:50,]
#' 
#' pheno = pData(dat)
#' edata = exprs(dat)
#' batch = pheno$batch
#' mod = model.matrix(~as.factor(cancer), data=pheno)
#' 
#' # parametric adjustment
#' combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
#' 
#' # non-parametric adjustment, mean-only version
#' combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
#' 
#' # reference-batch version, with covariates
#' combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
#' 
#' @export
#' 

ComBat <- function(dat, batch, mod=NULL, par.prior=TRUE,prior.plots=FALSE,mean.only=FALSE,ref.batch=NULL) {
  # make batch a factor and make a set of indicators for batch
  if(mean.only==TRUE){cat("Using the 'mean only' version of ComBat\n")}
  if(length(dim(batch))>1){stop("This version of ComBat only allows one batch variable")}  ## to be updated soon!
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  if (!is.null(ref.batch)){ # check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))){stop("reference level ref.batch is not one of the levels of the batch variable")}
    cat("Using batch =",ref.batch, "as a reference batch (this batch won't change)\n")
    ref = which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref]=1
  }else{ref=NULL}
  cat("Found",nlevels(batch),'batches\n')
  
  # A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch){batches[[i]] <- which(batch == levels(batch)[i])} # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){mean.only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n.array <- sum(n.batches)
  
  #combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  # check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
  # Number of covariates or covariate levels
  cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  # Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")}
    if(ncol(design)>(n.batch+1)){
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")}}
  }
  
  ## Check for missing values
  NAs = any(is.na(dat))
  if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
  #print(dat[1:2,])
  
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){
    B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
  }else{
    B.hat=apply(dat,1,Beta.NA,design)
  }
  
  ######## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  }else{
    grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
  }
  
  ######## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)){
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ]%*%B.hat))^2)%*%rep(1/n.batches[ref],n.batches[ref])
    }else{
      var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
    }
  }else{
    if(!is.null(ref.batch)){
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- apply(ref.dat-t(design[batches[[ref]], ]%*%B.hat),1,var,na.rm=TRUE)
    }else{
      var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=TRUE)
    }
  }
  
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
    if(mean.only==TRUE){delta.hat <- rbind(delta.hat,rep(1,nrow(s.data)))}else{
      delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=TRUE))
    }
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
      if(mean.only){
        gamma.star <- rbind(gamma.star,postmean(gamma.hat[i,],gamma.bar[i],1,1,t2[i]))
        delta.star <- rbind(delta.star,rep(1,nrow(s.data)))
      }else{temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
            gamma.star <- rbind(gamma.star,temp[1,])
            delta.star <- rbind(delta.star,temp[2,])
      }
    }
  }else{
    cat("Finding nonparametric adjustments\n")
    for (i in 1:n.batch){
      if(mean.only){delta.hat[i,]=1}
      temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }
  
  if(!is.null(ref.batch)){
    gamma.star[ref,]=0  ## set reference batch mean equal to 0
    delta.star[ref,]=1  ## set reference batch variance equal to 1
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
  
  ##### tiny change still exist when tested on bladder data
  #### total sum of change within each batch around 1e-15 
  ##### (could be computational system error).  
  ##### Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  
  return(bayesdata)
  
}
