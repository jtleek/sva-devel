#' A function for estimating surrogate variables with the two step approach of Leek and Storey 2007
#' 
#' This function is the implementation of the two step approach for estimating surrogate
#' variables proposed by Leek and Storey 2007 PLoS Genetics. This function is primarily
#' included for backwards compatibility. Newer versions of the sva algorithm are available
#' through \code{\link{sva}}, \code{\link{svaseq}}, with low level functionality available
#' through \code{\link{irwsva.build}} and \code{\link{ssva}}.
#' 
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param batch A factor variable giving the known batch levels
#' @param ... Other arguments to the \code{\link{sva}} function. 
#'
#' @return psva.D Data with batch effect removed but biological heterogeneity preserved
#' 
#' @examples 
#' library(bladderbatch)
#' library(limma)
#' data(bladderdata)
#' dat <- bladderEset[1:50,]
#' 
#' pheno = pData(dat)
#' edata = exprs(dat)
#' batch = pheno$batch
#' batch.fac = as.factor(batch)
#' 
#' psva_data <- psva(edata,batch.fac)
#' 
#' @author Elana J. Fertig 
#' 
#' @export
#' 

psva <- function(dat, batch, ...) {
  
  # convert input class / categorical variables to a standard model matrix 
  if (class(batch)=='factor' | class(batch)=='character') {
    mod <- sva.class2Model(batch) 
  } else {
    stop('Invalid batch type for psva (require factor or character):', 
         class(batch))
  } 
  
  # find SV's
  psva.SV <- sva(dat=dat, mod=mod, ...)
  colnames(psva.SV$sv) <- paste('sv',1:ncol(psva.SV$sv))
  
  # fit data
  psva.fit <- lmFit(dat, cbind(mod,psva.SV$sv))
  
  # batch corrected data 
  psva.D <- sweep(psva.fit$coefficients[,paste('sv',1:ncol(psva.SV$sv))]%*%
                    t(psva.SV$sv), 1, psva.fit$coefficients[,"(Intercept)"], FUN="+")
  
  return(psva.D)

}