\name{irwsva.build}
\alias{irwsva.build}
\title{Build surrogate variables with an iterative algorithm from gene expression and model data}
\description{
  Construct a specified number of surrogate variables from a gene 
  expression data set and a fixed model. 
}
\usage{
  irwsva.build(dat, mod, mod0=NULL,n.sv,B=5)
}
\arguments{
  \item{dat}{A m genes by n arrays matrix of expression data}
  \item{mod}{A n by k model matrix corresponding to the primary model fit   
  (see model.matrix)}
  \item{mod0}{A n by k0 model matrix corresponding to the null model 
  to be compared to  mod.}
  \item{n.sv}{The number of surrogate variables to construct.}
  \item{B}{ The number of iterations of the algorithm to perform.}
}
\details{
  The IRW-SVA estimation algorithm is described in Leek and Storey (2008).The
  basic idea is to estimate surrogate variables based on the subset of rows
  affected by unmodeled dependence, but not affected by the main variable
  parameterized in mod. 
}

\value{
  A list containing:
  \item{sv}{A n by n.sv matrix where each column is a distinct surrogate 
  variable (the main quantity of interest)}
  \item{pprob.gam}{A vector with the posterior probability estimates that
  each row is affected by dependence.} 
  \item{pprob.b}{A vector with the posterior probabiliity estimates that each
  row is affected by the variables in mod, but not in mod0.}
  \item{n.sv}{The number of suggorate variables estimated. }		
}

\references{
   Leek JT and Storey JD. (2008) A general framework for multiple testing 
   dependence. Proceedings of the National Academy of Sciences, 105: 
   18718-18723. \url{http://www.biostat.jhsph.edu/~jleek/publications.html}

   Leek JT and Storey JD. (2007) Capturing heterogeneity in gene expression 
   studies by surrogate variable analysis. PLoS Genetics, 3: e161. 
   \url{http://www.biostat.jhsph.edu/~jleek/publications.html}
 
   sva Vignette \url{http://www.biostat.jhsph.edu/~jleek/sva/}
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}, John Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{sva}}, \code{\link{num.sv}}, \code{\link{twostepsva.build}}, \code{\link{ComBat}}}
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

  ## Construct the surrogate variables 
  svaobj <- irwsva.build(edata,mod,mod0,n.sv=1)

  ## Include them in a downstream analysis

  mod.sv <- cbind(mod,svaobj$sv)
  mod0.sv <- cbind(mod0,svaobj$sv)
  adjusted.pvals <- f.pvalue(edata,mod.sv,mod0.sv)
  
  
 }
}
\keyword{misc}

