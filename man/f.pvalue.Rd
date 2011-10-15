\name{f.pvalue}
\alias{f.pvalue}
\title{Quickly calculate F-test p-values}
\description{
  Calculate p-values from a parametric F-test comparing the
  models mod and mod0 for each row of the data set dat. 
}
\usage{
  f.pvalue(dat,mod,mod0)
}
\arguments{
  \item{dat}{A m genes by n arrays matrix of expression data}
 \item{mod}{A n by k model matrix corresponding to the primary model fit   
  (see model.matrix)}
  \item{mod0}{A n by k0 model matrix corresponding to the null model 
  to be compared to  mod.}
}
\details{
  The data for test i should be in the ith row of dat, if mod0
  is null, the first column of mod is used as the null model. 
}

\value{
  \item{p}{A vector of p-values one for each row of dat.}
}


\author{Jeffrey T. Leek \email{jleek@jhsph.edu}, John Storey \email{jstorey@princeton.edu}}

\examples{
  \dontrun{
  ## Load data
  data(svadat)
  grp <- rep(c(0,1),each=10)
  mod <- model.matrix(~1+grp)
  mod0 <- cbind(rep(1,20))
  
  #Calculate the p-values
  p <- f.pvalue(dat,mod,mod0)
  hist(p)
  
 }
}
\keyword{misc}

