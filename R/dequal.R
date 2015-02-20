# ' A function to calculate and inspect RNA quality bias in differential expression analysis
#'
#' This function visualizes the effect of RNA Integrity Number (RIN) and
#' the outcome of interest on gene expression levels. Correlation between 
#' these two effects suggests that the differential expression (DE) analysis
#' may be confounded by quality. 
#'
#' @param Y Gene expression matrix, rows = genes, columns = samples
#' @param mod model matrix for the differential expression analysis, e.g. the output of model.matrix(~x)
#' @param RIN vector of RNA integrity numbers, one per sample
#' @param coiIndex covariate of interest index; i.e. which column of mod is
#'		the outcome of interest?
#' @param xlab label of x-axis. Default: "DE Effect"
#' @param verbose return messages? Default: TRUE
#' @param Contour plot null contours? Default: FALSE
#' @param Lines add lines (observed, null, and grid)? Default: TRUE
#' @param seed used for sampling the null data
#' @return DEqual plot and correlation between RIN effect of COI effect
#' @export

DEqualPlot = function(Y, mod, RIN, coiIndex=2, 
	xlab="DE Effect", verbose=TRUE, Contour=FALSE, 
	Lines=TRUE, seed=NULL, ...) {
	
	require(limma)
	require(MASS)
	
	cat("Fitting models.\n")
	# canidate model
	fit = lmFit(Y, mod)
	cf = fit$coef[,coiIndex]
	
	# rin model
	modRin = model.matrix(~RIN)
	fitRin = lmFit(Y, modRin)
			
	## random
	cat("Adding null info.\n")

	if(!is.null(seed)) set.seed(seed)
	ystar = matrix(rnorm(nrow(Y)*ncol(Y)), nc = ncol(Y))
	fitStar = lmFit(ystar, mod)
	fitStarRin = lmFit(ystar, modRin)
	if(Contour) kde = kde2d(fitStar$coef[,coiIndex], fitStarRin$coef[,2], n=50)

	cat("Plotting.\n")
	plot(cf, fitRin$coef[,2],ylab="RIN Effect", xlab=xlab, ...)
	if(Contour) contour(kde,add=TRUE,lwd=2)
	if(Lines) {
		abline(lm( fitRin$coef[,2] ~ cf),	lwd=3)
		abline(lm( fitStarRin$coef[,2] ~ fitStar$coef[,coiIndex]),lwd=3, lty=2)
		abline(h=0,v=0,lty=3)
	}
	return(cor(fitRin$coef[,2],cf))
}
