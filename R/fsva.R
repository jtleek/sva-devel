#' A function for performing frozen surrogate variable analysis as proposed in Parker, Corrada Bravo and Leek 2013
#'
#' This function performs frozen surrogate variable analysis as described in Parker, Corrada Bravo and Leek 2013.
#' The approach uses a training database to create surrogate variables which are then used
#' to remove batch effects both from the training database and a new data set for prediction purposes.
#' For inferential analysis see \code{\link{sva}}, \code{\link{svaseq}}, with low level functionality available
#' through \code{\link{irwsva.build}} and \code{\link{ssva}}.
#'
#' @param dbdat A m genes by n arrays matrix of expression data from the database/training data
#' @param mod The model matrix for the terms included in the analysis for the training data
#' @param sv The surrogate variable object created by running sva on dbdat using mod.
#' @param newdat (optional) A set of test samples to be adjusted using the training database
#' @param method If method ="fast" then the SVD is calculated using an online approach, this may introduce slight bias. If method="exact" the exact SVD is calculated, but will be slower
#'
#' @return db An adjusted version of the training database where the effect of batch/expression heterogeneity has been removed
#' @return new An adjusted version of the new samples, adjusted one at a time using the fsva methodology.
#' @return newsv Surrogate variables for the new samples
#'
#' @examples
#' library(bladderbatch)
#' library(pamr)
#' data(bladderdata)
#' dat <- bladderEset[1:50,]
#'
#' pheno = pData(dat)
#' edata = exprs(dat)
#'
#' set.seed(1234)
#' trainIndicator = sample(1:57,size=30,replace=FALSE)
#' testIndicator = (1:57)[-trainIndicator]
#' trainData = edata[,trainIndicator]
#' testData = edata[,testIndicator]
#' trainPheno = pheno[trainIndicator,]
#' testPheno = pheno[testIndicator,]
#'
#' mydata = list(x=trainData,y=trainPheno$cancer)
#' mytrain = pamr.train(mydata)
#' table(pamr.predict(mytrain,testData,threshold=2),testPheno$cancer)
#'
#' trainMod = model.matrix(~cancer,data=trainPheno)
#' trainMod0 = model.matrix(~1,data=trainPheno)
#' trainSv = sva(trainData,trainMod,trainMod0)
#'
#' fsvaobj = fsva(trainData,trainMod,trainSv,testData)
#' mydataSv = list(x=fsvaobj$db,y=trainPheno$cancer)
#' mytrainSv = pamr.train(mydataSv)
#' table(pamr.predict(mytrainSv,fsvaobj$new,threshold=1),testPheno$cancer)
#'
#' @export
#'
fsva <- function(dbdat,mod,sv,newdat=NULL,method=c("fast","exact")){

  method <- match.arg(method)
  ndb = dim(dbdat)[2]
  nnew = dim(newdat)[2]
  nmod = dim(mod)[2]
  ntot = ndb + nnew
  n.sv<-sv$n.sv

  # Adjust the database
  mod = cbind(mod,sv$sv)
  gammahat = (dbdat %*% mod %*% solve(t(mod) %*% mod))[,(nmod+1):(nmod + sv$n.sv)]
  db = dbdat - gammahat %*% t(sv$sv)
  adjusted=NA

  if(!is.null(newdat)){
    # Calculate the new surrogate variables
    if(method=="fast"){
      wts <- (1-sv$pprob.b)*sv$pprob.gam
      WX <- wts*dbdat
      svd.wx = svd(t(scale(t(WX),scale=FALSE)))
      D <- svd.wx$d[1:n.sv]
      U <- svd.wx$u[,1:n.sv]
      P <- t(wts*t(1/D * t(U)))
      newV <- P %*% newdat
      sgn = rep(NA,n.sv)
      for(j in 1:sv$n.sv){
        if(sv$n.sv>1){
          sgn[j] = sign(cor(svd.wx$v[1:ndb,j],sv$sv[1:ndb,j]))
        }
        if(sv$n.sv==1){
          sgn[j] = sign(cor(svd.wx$v[1:ndb,j],sv$sv[1:ndb]))
        }
      }
      newV = newV*sgn
      newV = t(newV)
      newV = scale(newV)/sqrt(dim(newV)[1])
      newV = t(newV)


    }else if(method=="exact"){
      newV = matrix(nrow=nnew,ncol=n.sv)
      for(i in 1:nnew){
        tmp = cbind(dbdat,newdat[,i])
        tmpd = (1-sv$pprob.b)*sv$pprob.gam*tmp
        ss = svd(t(scale(t(tmpd),scale=FALSE)))
        sgn = rep(NA,sv$n.sv)
        for(j in 1:sv$n.sv){
          if(sv$n.sv>1){
            sgn[j] = sign(cor(ss$v[1:ndb,j],sv$sv[1:ndb,j]))
          }
          if(sv$n.sv==1){
            sgn[j] = sign(cor(ss$v[1:ndb,j],sv$sv[1:ndb]))
          }
        }
        newV[i,]<-ss$v[(ndb+1),1:sv$n.sv]*sgn
      }
      newV = t(newV)
    }

    adjusted = newdat - gammahat %*% newV
    newV = t(newV)
  }

  return(list(db=db,new=adjusted,newsv = newV))
}
