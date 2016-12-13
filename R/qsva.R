#' A function for computing quality surrogate variables (qSVs)
#'
#' This function computes quality surrogate variables (qSVs)
#' from the library-size- and read-length-normalized
#' degradation matrix for subsequent RNA quality correction
#'
#' @param degradationMatrix the normalized degradation matrix, region by sample
#' @param mod (Optional) statistical model used in DE analysis
#'
#' @return the qSV adjustment variables
#'
#' @examples
#' degCovAdj = read.degradation.matrix(
#'  covFiles = list.files(bwPath,full.names=TRUE),
#'  sampleNames = list.files(bwPath), readLength = 76,
#'  totalMapped = rep(100e6,5),type="bwtool")
#'  head(degCovAdj)
#' qsva(degCovAdj)
#'
#' @export
qsva <- function(degradationMatrix,
		mod = matrix(1,nc=1,nr = ncol(degradationMatrix))) {

	# do PCA
	degPca = prcomp(t(log2(degradationMatrix+1)))

	## how many PCs?
	k = num.sv(log2(degradationMatrix+1), mod)

	# return qSVS
	degPca$x[,1:k]
}
