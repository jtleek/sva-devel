#' sva: a package for removing artifacts from microarray and sequencing data 
#'
#' sva has functionality to estimate and remove artifacts from high dimensional data
#' the \code{\link{sva}} function can be used to estimate artifacts from microarray data
#' the \code{\link{svaseq}} function can be used to estimate artifacts from count-based
#' RNA-sequencing (and other sequencing) data. The \code{\link{ComBat}} function can be
#' used to remove known batch effecs from microarray data. The \code{\link{fsva}} function
#' can be used to remove batch effects for prediction problems. 
#' 
#' 
#' A vignette is available by typing \code{browseVignettes("sva")} in the R prompt.
#' 
#' @references For the package: Leek JT, Johnson WE, Parker HS, Jaffe AE, and Storey JD. (2012) The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics DOI:10.1093/bioinformatics/bts034
#' @references For sva: Leek JT and Storey JD. (2008) A general framework for multiple testing dependence. Proceedings of the National Academy of Sciences , 105: 18718-18723. 
#' @references For sva: Leek JT and Storey JD. (2007) Capturing heterogeneity in gene expression studies by `Surrogate Variable Analysis'. PLoS Genetics, 3: e161.
#' @references For Combat: Johnson WE, Li C, Rabinovic A (2007) Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics,  8 (1), 118-127
#' @references For svaseq: Leek JT (2014) svaseq: removing batch and other artifacts from count-based sequencing data. bioRxiv doi: TBD
#' @references For fsva: Parker HS, Bravo HC, Leek JT (2013) Removing batch effects for prediction problems with frozen surrogate variable analysis arXiv:1301.3947
#' 
#' @docType package
#' @author Jeffrey T. Leek, W. Evan Johnson, Hilary S. Parker, Andrew E. Jaffe, John D. Storey
#' @name sva
#' @useDynLib sva monotone
NULL
