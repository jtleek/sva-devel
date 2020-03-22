library(sva)
library(testthat)
context("ComBat-Seq tests")


### Run ComBat on a simple simulated dataset
test_that("check ComBat-Seq output with several different parameters on simulated genes",{
  count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
  batch <- c(rep(1, 4), rep(2, 4))
  group <- rep(c(0,1), 4)
  
  ### Test choice of parameter
  mod_choice <- capture.output(ComBat_seq(count_matrix, batch, group=rep(1, 8), full_mod=TRUE))[2]
  expect_equal(mod_choice, "Using null model in ComBat-seq.")

  ### Test error message
  expect_error(ComBat_seq(count_matrix, batch, group=batch), "The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")
  expect_error(ComBat_seq(count_matrix, batch=c(1,rep(2,7))), "ComBat-seq doesn't support 1 sample per batch yet")
})
