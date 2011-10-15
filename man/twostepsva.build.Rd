\name{twostepsva.build}
\alias{twostepsva.build}
\title{Build surrogate variables from gene expression and model data}
\description{
  Construct a specified number of surrogate variables from a gene 
  expression data set based on the two-step algorithm of Leek and Storey (2007).
}
\usage{
  twostepsva.build(dat, mod, n.sv)
}
\arguments{
  \item{dat}{A m genes by n arrays matrix of expression data}
 \item{mod}{A n by k model matrix corresponding to the primary model fit   
  (see model.matrix)}
  \item{n.sv}{The number of surrogate variables to construct.}
}
\details{
  The SVA estimation algorithm is described in Leek and Storey (2007).The
  basic idea is to estimate surrogate variables based on the subset of rows
  affected by unmodeled dependence.
}

\value{
  A list containing:
  \item{sv}{A n by n.sv matrix where each column is a distinct surrogate 
  variable (the main quantity of interest)}
  \item{pprob.gam}{A vector indicating whether each row was used in the
  building of the surrogate variables. 1= row used, 0=not used.} 
  \item{pprob.b}{Null for two-step SVA, see irwsva.build for more info.}
  \item{n.sv}{The number of suggorate variables estimated. }		
}

\references{
   Leek JT and Storey JD. (2007) Capturing heterogeneity in gene expression 
   studies by surrogate variable analysis. PLoS Genetics, 3: e161. 
   \url{http://www.biostat.jhsph.edu/~jleek/publications.html}
 
  SVA Manual \url{http://www.biostat.jhsph.edu/~jleek/sva/manual.pdf}
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}, John Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{sva}}, \code{\link{num.sv}},\code{\link{irwsva.build}}}
\examples{
  \dontrun{
  ## Load data
  data(svadat)
  grp <- rep(c(0,1),each=10)
  mod <- model.matrix(~1+grp)
  mod0 <- cbind(rep(1,20))
  
  ## Construct the surrogate variables 
  svaobj <- twostepsva.build(svadat,mod,1)
  

  ## Include them in a downstream analysis
  mod.sv <- cbind(mod,svaobj$sv)
  mod0.sv <- cbind(mod0,svaobj$sv)
  adjusted.pvals <- f.pvalue(dat,mod.sv,mod0.sv)
  
  
 }
}
\keyword{misc}

