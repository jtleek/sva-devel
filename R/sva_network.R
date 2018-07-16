#' A function to adjust gene expression data before network inference
#' 
#' This function corrects a gene expression matrix prior to network inference by
#' returning the residuals after regressing out the top principal components. The number 
#' of principal components to remove can be determined using a permutation-based approach
#' using the "num.sv" function with method = "be"

#' @param dat The raw gene expression data matrix (with the variables in rows and samples in columns)
#' @param n.pc The number of principal components to remove 

#' @return dat.adjusted Cleaned gene expression data matrix with the top prinicpal components removed
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
#' n.pc = num.sv(edata, mod, method="be")
#' dat.adjusted = sva_network(edata, n.pc)
#' 
#' @export
#' 
sva_network=function(dat, n.pc){
	ss<-svd(dat - colMeans(dat))

	dat.adjusted=dat
	for (i in 1:dim(dat)[2]){
  		dat.adjusted[,i]<- lm(dat[,i] ~ ss$u[,1:n.pc])$residuals
	}
	return(dat.adjusted)
}