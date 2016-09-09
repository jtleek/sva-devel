#' A function for reading in coverage data from degradation-susceptible regions 
#'
#' This function reads in degradation regions to form 
#' a library-size- and read-length-normalized 
#' degradation matrix for subsequent RNA quality correction
#' 
#' @param covFiles coverage file(s) for degradation regions
#' @param sampleNames sample names; creates column names of degradation matrix
#' @param totalMapped how many reads per sample (library size normalization) 
#' @param readLength read length in base pairs (read length normalization)
#' @param normFactor common library size to normalize to; 80M reads as default
#' @param type whether input are individual `bwtool` output, `region_matrix` run on individual samples, or `region_matrix` run on all samples together
#'
#' @return the normalized degradation matrix, region by sample
#' 
#' @examples 
#' bwPath = system.file('extdata', 'bwtool', package = 'sva')
#' degCovAdj = read.degradation.matrix(
#'  covFiles = list.files(bwPath,full.names=TRUE),
#'  sampleNames = list.files(bwPath), L = 76, 
#'  totalMapped = rep(100e6,5),type="bwtool")
#'  

read.degradation.matrix <- function(covFiles, sampleNames,
	totalMapped, readLength= 100, normFactor = 80e6,
	type = c("bwtool","region_matrix_single","region_matrix_all"),
	mc.cores=1) {
	
	type <- match.arg(type)
	
	if(length(covFiles) != length(sampleNames) & 
		type %in% c("bwtool", "region_matrix_single")) 
			stop("Must provide one coverage file and sample name per sample for 'bwtool' or 'region_matrix_single' types")
	
	## read in data
	if(type == "bwtool") { # bwtool output
		degCov = mclapply(covFiles, read.delim,
			as.is=TRUE, colClasses = c(rep("NULL",9),
				"numeric"),header=TRUE,mc.cores=mc.cores)
		degCovMat =do.call("cbind",degCov)	
		colnames(degCovMat) = sampleNames
	} else if(type == "region_matrix_single") { # region w/ manifest
		degCov = mclapply(covFiles, read.delim,
			as.is=TRUE,header=TRUE,row.names=1,mc.cores=mc.cores)
		degCovMat =do.call("cbind",degCov)	
		colnames(degCovMat) = sampleNames
	} else if(type == "region_matrix_all") { # region w/ individual file
		degCov = read.delim(covFiles, as.is=TRUE,row.names=1)
		degCov = t(as.matrix(degCov))
		colnames(degCov) = sampleNames
	} else stop("type must be 'bwtool','region_matrix_single','region_matrix_all'")
		
	## normalize for library size and read length
	degCovMat = degCovMat/readLength # read length
	bg = matrix(rep(totalMapped/normFactor), 
		nc = ncol(degCovMat), 
		nr = nrow(degCovMat), byrow=TRUE)
	degCovAdj = degCovMat/bg
	
	# return
	return(degCovAdj)
}