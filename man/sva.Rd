\name{sva}
\alias{sva}
\title{Estimate surrogate variables with an iterative algorithm from gene expression and model data}
\description{
  Estimate surrogate variables are estimated using either the iteratively re-weighted
  surrogate variable analysis algorithm of Leek and Storey (2008) or the two-step
  algorithm of Leek and Storey (2007).
}
\usage{
  sva(dat, mod, mod0 = NULL,n.sv=NULL,method=c("irw","two-step"),vfilter=NULL, B=5)
}
\arguments{
  \item{dat}{A m genes by n arrays matrix of expression data}
  \item{mod}{A n by k model matrix corresponding to the primary model fit   
  (see model.matrix)}
  \item{mod0}{A n by k0 model matrix corresponding to the null model 
  to be compared to  mod.}
  \item{n.sv}{Optional. The number of surrogate variables to estimate, can be estimated
using the num.sv function}
  \item{method}{Choose between the iteratively re-weighted or two-step
  surrogate variable estimation algorithms.}
  \item{vfilter}{The number of most variable genes to use when building SVs, must be between 100 and m}
  \item{B}{ The number of iterations of the algorithm to perform.}
}
\details{
  Surrogate variable estimates are formed based on the algorithms
  in Leek and Storey (2007,2008). Surrogate variables can be included
  in a significance analysis to reduce dependence and confounding. 
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
 
  SVA Manual \url{http://www.biostat.jhsph.edu/~jleek/sva/manual.pdf}
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}, John Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{irwsva.build}},\code{\link{twostepsva.build}}, \code{\link{num.sv}}}
\examples{
  \dontrun{

  ## Load data
  data(svadat)
  grp <- rep(c(0,1),each=10)
  mod <- model.matrix(~1+grp)
  mod0 <- cbind(rep(1,20))


  ## Construct the surrogate variables 
  svaobj <- sva(svadat,mod,mod0,method="irw")


  ## Include them in a downstream analysis

  mod.sv <- cbind(mod,svaobj$sv)
  mod0.sv <- cbind(mod0,svaobj$sv)
  adjusted.pvals <- f.pvalue(dat,mod.sv,mod0.sv)
  
 }
}
\keyword{misc}

