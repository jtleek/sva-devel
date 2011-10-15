\name{fsva}
\alias{fsva}
\title{Single sample surrogate variable correction for prediciton problems.}
\description{
  fsva corrects training databases and performs "single-sample" correction
on new samples for prediction problems. The effect of surrogate variables
is removed from the training database which can then be used to build a predictor.
The new samples are corrected individually to account for batch effects when
the group status is unknown. 
}
\usage{
  fsva(dbdat,mod,sv,newdat=NULL)
}
\arguments{
  \item{dbdat}{A m genes by n arrays matrix of expression data from the database/training data}
  \item{mod}{The model matrix for the terms included in the analysis for the training data}
  \item{sv}{The surrogate variable object created by running sva on dbdat using mod. }
  \item{newdat}{(Optional)  A set of test samples to be adjusted using the training database}
}
\details{
  Frozen surrogate variable analysis (fsva) can be applied to remove batch effects
for prediction problems. 
}

\value{
  A list containing:
  \item{db}{An adjusted version of the training database where the effect of batch/expression heterogeneity has been removed)}
  \item{new}{An adjusted version of the new samples, adjusted one at a time using the fsva methodology. } 		
}

\references{
  SVA Manual \url{http://www.biostat.jhsph.edu/~jleek/sva/manual.pdf}
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}}
\seealso{\code{\link{sva}},\code{\link{irwsva.build}},\code{\link{twostepsva.build}}, \code{\link{num.sv}}}
\examples{
  \dontrun{

  ## Load data
  data(svadat)
  grp <- rep(c(0,1),each=10)
  mod <- model.matrix(~1+grp)
  mod0 <- cbind(rep(1,20))

  ## Construct the surrogate variables 
  svaobj <- sva(svadat,mod,mod0)

  ## Create some "new samples"
  newdat = svadat + matrix(rnorm(1000*20),nrow=1000)
  
  ## Correct the "new samples"
  data <- fsva(svadat,mod,sv,newdat)

  # Look at the new samples and corrected new samples
  par(mfrow=c(1,2))
  image(newdat)
  image(data$new)
 }
}
\keyword{misc}

