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
  fsva(dbdat,mod,sv,newdat=NULL,method=c("fast","exact"))
}
\arguments{
  \item{dbdat}{A m genes by n arrays matrix of expression data from the database/training data}
  \item{mod}{The model matrix for the terms included in the analysis for the training data}
  \item{sv}{The surrogate variable object created by running sva on dbdat using mod. }
  \item{newdat}{(Optional)  A set of test samples to be adjusted using the training database}
  \item{method}{If method = "fast" then the online SVD is used, this may introduce slight bias. If method="exact" the exact SVD is calculated, but will be slower.}
}
\details{
  Frozen surrogate variable analysis (fsva) can be applied to remove batch effects
for prediction problems. 
}

\value{
  A list containing:
  \item{db}{An adjusted version of the training database where the effect of batch/expression heterogeneity has been removed)}
  \item{new}{An adjusted version of the new samples, adjusted one at a time using the fsva methodology. } 
  \item{newsv}{The surrogate variables on the new samples}		
}

\references{
  sva Vignette \url{http://www.biostat.jhsph.edu/~jleek/sva}
}

\author{Jeffrey T. Leek \email{jleek@jhsph.edu}, Hilary S. Parker \email{hiparker@jhsph.edu}}
\seealso{\code{\link{sva}},\code{\link{irwsva.build}},\code{\link{twostepsva.build}}, \code{\link{num.sv}}}
\examples{
  \dontrun{

  ## Load data
  library(bladderbatch)
  library(pamr)
  data(bladderdata)
  
  ## Obtain phenotypic data
  pheno = pData(bladderEset)
  edata = exprs(bladderEset)
  batch = pheno$batch
  mod = model.matrix(~as.factor(cancer), data=pheno)

  ## Build the training and test sets
  set.seed(12354)
  trainIndicator = sample(1:57,size=30,replace=F)
  testIndicator = (1:57)[-trainIndicator]

  trainData = edata[,trainIndicator]
  testData = edata[,testIndicator]

  trainPheno = pheno[trainIndicator,]
  testPheno = pheno[testIndicator,]

  # Fit the sva model to the training set
  trainMod = model.matrix(~cancer,data=trainPheno)
  trainMod0 = model.matrix(~1,data=trainPheno)
  trainSv = sva(trainData,trainMod,trainMod0)  

  #fsva correct and train
  fsvaobj = fsva(dbdat=trainData,mod=trainMod,sv=trainSv,newdat=testData)
  mydataSv = list(x=fsvaobj$db,y=trainPheno$cancer)
  mytrainSv = pamr.train(mydataSv)
  table(pamr.predict(mytrainSv,fsvaobj$new,threshold=1),testPheno$cancer)

 }
}
\keyword{misc}

