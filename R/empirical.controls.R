#' A function for estimating the probability that each gene is an empirical control 
#' 
#' This function uses the iteratively reweighted surrogate variable analysis approach
#' to estimate the probability that each gene is an empirical control. 
#' 
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param mod The model matrix being used to fit the data
#' @param mod0 The null model being compared when fitting the data
#' @param n.sv The number of surogate variables to estimate
#' @param B The number of iterations of the irwsva algorithm to perform
#' @param type If type is norm then standard irwsva is applied, if type is counts, then the moderated log transform is applied first
#' 
#' @return pcontrol A vector of probabilites that each gene is a control. 
#' 
#' @export
#' 

empirical.controls <- function(dat, mod, mod0 = NULL,n.sv,B=5,type=c("norm","counts")) {
  type = match.arg(type)
  if((type=="counts") & any(dat < 0)){stop("empirical controls error: counts must be zero or greater")}
  if(type=="counts"){dat = log(dat + 1)}
  svobj = irwsva.build(dat, mod, mod0 = NULL,n.sv,B=5) 
  pcontrol = svobj$pprob.gam
  return(pcontrol)

}