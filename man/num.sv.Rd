\name{num.sv}
\alias{num.sv}
\title{Estimate number of surrogate variables to include in analysis}
\description{
  Estimate number of important surrogate variables from a gene expression 
  data set. 
}
\usage{
  num.sv(dat, mod,method=c("be","leek"),vfilter=NULL, B=20, sv.sig=0.10,seed=NULL)
}
\arguments{
  \item{dat}{A m genes by n arrays matrix of expression data}
  \item{mod}{A n by k model matrix corresponding to the primary model fit   
  (see model.matrix)}
  \item{method}{The method to use for estimating surrogate variables, for now
  there is only one option (based ib Buja and Eyuboglu 1992). }
   \item{vfilter}{The number of most variable genes to use when building SVs, must be between 100 and m}
  \item{B}{The number of null iterations to perform. Only used when method="be"}
  \item{sv.sig}{The significance cutoff for eigengenes. Only used when method="be"}
  \item{seed}{A numeric seed for reproducible results. Optional, only used when method="be"} 
}
\details{
  The model matrix should include a column for an intercept. num.sv
  estimates the number of surrogate variables to include in the analysis
  as described in Leek and Storey (2007), based on the permutation
  test of Buja and Eyuboglu (1992). 
}


\value{
  A list containing:
  \item{n.sv}{The number of significant surrogate variables}

}
\references{
   Buja A and Eyuboglu N. (1992) Remarks on parrallel analysis. Multivariate
   Behavioral Research, 27(4), 509-540
 
   Leek JT and Storey JD. (2008) A general framework for multiple testing 
   dependence. Proceedings of the National Academy of Sciences, 105: 
   18718-18723. \url{http://www.biostat.jhsph.edu/~jleek/publications.html}

   Leek JT and Storey JD. (2007) Capturing heterogeneity in gene expression 
   studies by surrogate variable analysis. PLoS Genetics, 3: e161. 
   \url{http://www.biostat.jhsph.edu/~jleek/publications.html}
 
  sva Vignette \url{http://www.biostat.jhsph.edu/~jleek/sva/}
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}, John Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{sva}}, \code{\link{irwsva.build}}, \code{\link{twostepsva.build}} }
\examples{
  \dontrun{
  
  ## Load data
  library(bladderbatch)
  data(bladderdata)
  
  ## Obtain phenotypic data
  pheno = pData(bladderEset)
  edata = exprs(bladderEset)
  batch = pheno$batch
  mod = model.matrix(~as.factor(cancer), data=pheno)
  mod0 = model.matrix(~1, data=pheno)

  ## Calculate the number of surrogate variables 
  xx <- num.sv(edata,mod)
  xx
 
 }
}

\keyword{misc}


