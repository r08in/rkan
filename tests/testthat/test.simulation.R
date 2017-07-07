library(rkan)
setwd("E:/R/rkan")
#setwd("D:/RProject/RobustCD")
source("Simulation/Simulation.R")
context("test simulation result")

test_that("simulation result of rkan for Case C ",{
  beta1 <- c(1, 0.5, 0, 0)
  beta2 <- c(1, -0.9, -1.3, -0.5)
  exam1 <- list(n=100, p=16, g=rep(4,4), beta=c(beta1,beta2,rep(0,8)))
  #exam1 <- list(n=100, p=400, g=rep(4,100), beta=c(beta1,beta2,rep(0,392)))
  test.res <- simulation(L=50, n=exam1$n, beta=exam1$beta, g=exam1$g, model=c("C"), method="rkan",nk=8)
  expect_equal(test.res[[1]][1:6], list(model="C",CFR=38,CFR2=44,OFR=62, 
                                   PDR=100,FDR=35.9))
})