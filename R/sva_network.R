
## sva_network() is a function to adjust expression data before network inference
## The function inputs raw gene expression data (raw) and the number of principal components to remove (n.pc) and returns a cleaned matrix with those prinicpal components removed.
## The number of prinicpal components can be determined using the "num.sv" function: n.pc=num.sv(t(dat), mod,method="be")

sva_network=function(dat, n.pc){

	## singular value decomposition
	ss<-svd(dat - colMeans(dat))

	## use residuals from top n.pc principal components
	dat.adjusted=dat
	for (i in 1:dim(dat)[2]){
  		dat.adjusted[,i]<- lm(dat[,i] ~ ref.ss$u[,1:n.pc])$residuals
	}
	return(dat.adjusted)

}