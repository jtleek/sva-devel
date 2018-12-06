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
#' @param BPPARAM (Optional) BiocParallelParam for parallel operation
#'
#' @return data A probe x sample genomic measure matrix, adjusted for batch effects.
#'
#' @importFrom graphics lines par
#' @importFrom stats cor density dnorm model.matrix pf ppoints prcomp predict
#' qgamma qnorm qqline qqnorm qqplot smooth.spline var
#' @importFrom utils read.delim
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

ComBat <- function(dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                   mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam")) {
    ## coerce dat into a matrix
    dat <- as.matrix(dat)
    
    ## find genes with zero variance in any of the batches
    batch <- as.factor(batch)
    zero.rows.lst <- lapply(levels(batch), function(batch_level){
      which(apply(dat[, batch==batch_level], 1, function(x){var(x)==0}))
    })
    zero.rows <- Reduce(union, zero.rows.lst)
    keep.rows <- setdiff(1:nrow(dat), zero.rows)
    
    if (length(zero.rows) > 0) {
      cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n", length(zero.rows)))
      # keep a copy of the original data matrix and remove zero var rows
      dat.orig <- dat
      dat <- dat[keep.rows, ]
    }
  
    ## make batch a factor and make a set of indicators for batch
    if(mean.only==TRUE){
        message("Using the 'mean only' version of ComBat")
    }
    if(length(dim(batch))>1){
        stop("This version of ComBat only allows one batch variable")
    }  ## to be updated soon!
    batchmod <- model.matrix(~-1+batch)  
    if (!is.null(ref.batch)){
        ## check for reference batch, check value, and make appropriate changes
        if (!(ref.batch%in%levels(batch))) {
            stop("reference level ref.batch is not one of the levels of the batch variable")
        }
        cat("Using batch =",ref.batch, "as a reference batch (this batch won't change)\n")
        ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
        batchmod[,ref] <- 1
    } else {
        ref <- NULL
    }
    message("Found", nlevels(batch), "batches")
  
    ## A few other characteristics on the batches
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels(batch)[i])
    } # list of samples in each batch  
    n.batches <- sapply(batches, length)
    if(any(n.batches==1)){
        mean.only=TRUE
        message("Note: one batch has only one sample, setting mean.only=TRUE")
    }
    n.array <- sum(n.batches)
    ## combine batch variable and covariates
    design <- cbind(batchmod,mod)
  
    ## check for intercept in covariates, and drop if present
    check <- apply(design, 2, function(x) all(x == 1))
    if(!is.null(ref)){
        check[ref] <- FALSE
    } ## except don't throw away the reference batch indicator
    design <- as.matrix(design[,!check])
  
    ## Number of covariates or covariate levels
    message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  
    ## Check if the design is confounded
    if(qr(design)$rank < ncol(design)) {
        ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
        if(ncol(design)==(n.batch+1)) {
            stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
        }
        if(ncol(design)>(n.batch+1)) {
            if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
                stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
            } else {
                stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
            }
        }
    }
  
    ## Check for missing values
    NAs <- any(is.na(dat))
    if(NAs){
        message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
    ## print(dat[1:2,])
  
    ##Standardize Data across genes
    cat('Standardizing Data across genes\n')
    if (!NAs){
        B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
    } else { 
        B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
    }
    
    ## change grand.mean for ref batch
    if(!is.null(ref.batch)){
        grand.mean <- t(B.hat[ref, ])
    } else {
        grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
    }
  
    ## change var.pooled for ref batch
    if (!NAs){
        if(!is.null(ref.batch)) {
            ref.dat <- dat[, batches[[ref]]]
            var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
        } else {
            var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
        }
    } else {
        if(!is.null(ref.batch)) {
            ref.dat <- dat[, batches[[ref]]]
            var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
        } else {
            var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
        }
    }
  
    stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
    if(!is.null(design)){
        tmp <- design
        tmp[,c(1:n.batch)] <- 0
        stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
    }  
    s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
    ##Get regression batch effect parameters
    message("Fitting L/S model and finding priors")
    batch.design <- design[, 1:n.batch]
    if (!NAs){
      gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                             as.matrix(s.data)))
    } else{
        gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
    }
    delta.hat <- NULL
    for (i in batches){
        if(mean.only==TRUE) {
            delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
        } else {
            delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
        }
    }
  
    ##Find Priors
    gamma.bar <- rowMeans(gamma.hat)
    t2 <- rowVars(gamma.hat)
    a.prior <- apply(delta.hat, 1, aprior) # FIXME 
    b.prior <- apply(delta.hat, 1, bprior) # FIXME
  
    ## Plot empirical and parametric priors

    if (prior.plots && par.prior) {
        par(mfrow=c(2,2))
        
        ## Top left
        tmp <- density(gamma.hat[1,])
        plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
        xx <- seq(min(tmp$x), max(tmp$x), length=100)
        lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
        
        ## Top Right
        qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
        qqline(gamma.hat[1,], col=2)
        
        ## Bottom Left
        tmp <- density(delta.hat[1,])
        xx <- seq(min(tmp$x), max(tmp$x), length=100)
        tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
        plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
             main=expression(paste("Density Plot of First Batch ", hat(delta))))
        lines(tmp1, col=2)

        ## Bottom Right
        invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
        qqplot(invgam, delta.hat[1,],
               main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
               ylab="Sample Quantiles", xlab="Theoretical Quantiles")
        lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
    }
      
    ## Find EB batch adjustments

    gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
    if (par.prior) {
        message("Finding parametric adjustments")
        results <- bplapply(1:n.batch, function(i) {
            if (mean.only) {
                gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
                delta.star <- rep(1, nrow(s.data))
            }
            else {
                temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                               delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                               b.prior[i])
                gamma.star <- temp[1, ]
                delta.star <- temp[2, ]
            }
            list(gamma.star=gamma.star, delta.star=delta.star)
        }, BPPARAM = BPPARAM)
        for (i in 1:n.batch) {
            gamma.star[i,] <- results[[i]]$gamma.star
            delta.star[i,] <- results[[i]]$delta.star
        }
    }
    else {
        message("Finding nonparametric adjustments")
        results <- bplapply(1:n.batch, function(i) {
            if (mean.only) {
                delta.hat[i, ] = 1
            }
            temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                               gamma.hat[i, ], delta.hat[i, ])
            list(gamma.star=temp[1,], delta.star=temp[2,])
        }, BPPARAM = BPPARAM)
        for (i in 1:n.batch) {
            gamma.star[i,] <- results[[i]]$gamma.star
            delta.star[i,] <- results[[i]]$delta.star
        }
    }
    
    if(!is.null(ref.batch)){
        gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
        delta.star[ref,] <- 1  ## set reference batch variance equal to 1
    }
    
    ## Normalize the Data ###
    message("Adjusting the Data\n")
  
    bayesdata <- s.data
    j <- 1
    for (i in batches){
        bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
        j <- j+1
    }
    
    bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME
  
    ## Do not change ref batch at all in reference version
    if(!is.null(ref.batch)){
        bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
    }
    
    ## put genes with 0 variance in any batch back in data
    if (length(zero.rows) > 0) {
      dat.orig[keep.rows, ] <- bayesdata
      bayesdata <- dat.orig
    }
    
    return(bayesdata)
}
