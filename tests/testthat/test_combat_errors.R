library(sva)
library(testthat)
context("ComBat Error Check")

test_that("Check ComBat errors for too many batch variables and confounded designs",{
  # set up number of batches (nb), total sample size (n), and batch/treatment vaariables 
  nb <- 2  
  ns <- 2 # one of nb or ns needs to be even (see treat variable)
  n <-  ns*nb
  batch <- as.factor(rep(1:nb, each=ns))
  treat <- rep(c(0,1), n/2)
  
  #generate random data
  X <- cbind(model.matrix(~-1+batch), treat) 
  beta <- c(3*(1:nb-1), 2)
  sim <-  100
  y <-  matrix(rnorm(n*sim,X%*%beta),nrow=sim, byrow=TRUE)
  
  ## user gives multiple batch variables
  expect_error(ComBat(y,cbind(batch, batch)),"This version of ComBat only allows one batch variable")
  
  ## test confounded errors:
  expect_error(ComBat(y,batch, batch), "The covariate is confounded with batch! Remove the covariate and rerun ComBat") 
  expect_error(ComBat(y,batch, cbind(treat,batch)), "At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
  expect_error(ComBat(y,batch, cbind(treat,treat)), "The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
  
  ## test ref batch input error
  expect_error(ComBat(y,batch, ref.batch="yellow"), "reference level ref.batch is not one of the levels of the batch variable")
  expect_error(ComBat(y,batch, ref.batch=100), "reference level ref.batch is not one of the levels of the batch variable")
  expect_error(ComBat(y,batch, ref.batch=rep(0,3)), "reference level ref.batch is not one of the levels of the batch variable")
})
