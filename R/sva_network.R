#' A function to adjust gene expression data before network inference
#' 
#' This function corrects a gene expression matrix prior to network inference by
#' returning the residuals after regressing out the top principal components. 
#' The number of principal components to remove can be determined using a 
#' permutation-based approach using the "num.sv" function with method = "be"

#' @param dat The uncorrected normalized gene expression data matrix with samples in rows and genes in columns
#' @param n.pc The number of principal components to remove 

#' @return dat.adjusted Cleaned gene expression data matrix with the top prinicpal components removed
#' 
#' @examples 
#' library(bladderbatch)
#' data(bladderdata)
#' dat <- bladderEset[1:5000,]
#' 
#' edata = exprs(dat)
#' mod = matrix(1, nrow = dim(dat)[2], ncol = 1)
#' 
#' n.pc = num.sv(edata, mod, method="be")
#' dat.adjusted = sva_network(t(edata), n.pc)
#' 
#' @export
#' 
sva_network=function(dat, n.pc){
	ss<-svd(dat - colMeans(dat))

	dat.adjusted=dat
	for (i in seq_len(dim(dat))[2]){
  		dat.adjusted[,i]<- lm(dat[,i] ~ ss$u[,1:n.pc])$residuals
	}
	return(dat.adjusted)
}
