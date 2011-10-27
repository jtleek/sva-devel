\name{ComBat}
\alias{ComBat}
\title{Adjust for batch effects using empirical Bayes framework.}
\description{
ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for batch effects.  Users are returned an expression matrix that has been corrected for batch effects.
}
\usage{
	ComBat(dat, batch, mod, numCovs=NULL, par.prior=TRUE,prior.plots=FALSE)
}
\arguments{
  \item{dat}{Genomic measure matrix (dimensions probe x sample) - for example, expression matrix}
  \item{batch}{Batch covariate (multiple batches allowed)}
  \item{mod}{Model matrix for outcome of interest and other covariates besides batch}
  \item{numCovs}{(Optional) Vector containing the column numbers of the continuous covariates in the model matrix, or NULL if no continuous covariates are used}
  \item{par.prior}{(Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used}
  \item{prior.plots}{(Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric estimate}
}
\details{
  ComBat can be applied to genomic measures when batch is known to remove the effect of batch on the data using an empirical Bayesian framework.  It was described in Johnson et al. 2007.
}

\value{
	A probe x sample genomic measure matrix, adjusted for batch effects.
 }

\references{
  Johnson et al. 2007 \url{http://www.ncbi.nlm.nih.gov/pubmed/16632515}
}

\author{W. Evan Johnson \email{evan@stat.byu.edu}}
\seealso{\code{\link{sva}},\code{\link{irwsva.build}},\code{\link{twostepsva.build}}, \code{\link{num.sv}}}
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
  
  ## Correct for batch using ComBat
  combat_edata = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

 }
}
\keyword{misc}
