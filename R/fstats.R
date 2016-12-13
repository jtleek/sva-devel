#' A function for quickly calculating f statistics for use in sva
#'
#' This function does simple linear algebra to calculate f-statistics
#' for each row of a data matrix comparing the nested models
#' defined by the design matrices for the alternative (mod) and and null (mod0) cases.
#' The columns of mod0 must be a subset of the columns of mod.
#'
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param mod The model matrix being used to fit the data
#' @param mod0 The null model being compared when fitting the data
#'
#' @return fstats A vector of F-statistics one for each row of dat.
#'
#' @examples
#' library(bladderbatch)
#' data(bladderdata)
#' dat <- bladderEset[1:50,]
#'
#' pheno = pData(dat)
#' edata = exprs(dat)
#' mod = model.matrix(~as.factor(cancer), data=pheno)
#' mod0 = model.matrix(~1,data=pheno)
#'
#' fs <- fstats(edata, mod, mod0)
#'
#' @export
#'
fstats <- function(dat,mod,mod0){
  # A function for calculating F-statistics
  # on the rows of dat, comparing the models
  # mod (alternative) and mod0 (null).
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0,m)
  Id <- diag(n)

  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))

  rss1 <- (resid*resid) %*% rep(1,n)
  rss0 <- (resid0*resid0) %*% rep(1,n)

  fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
  return(fstats)
}
