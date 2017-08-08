set_parameter = function(x, y, nlambda1, nlambda2, nk, lambda1.min=0.05,lambda2.min=1e-03, k.min, beta0, w0) {
  n <- length(y)
  p <- dim(x)[2]
  lambda1Max <- max(abs(t(x) %*% y/n) * abs(beta0))  # max |betaj|*|xj'y/n|
  lambda2Max <- max(y^2 * abs(1 - w0)/n)
  lambda1 <- logSeq2(lambda1Max, lambda1Max * lambda1.min, nlambda1)
  lambda2 <- logSeq2(lambda2Max, lambda2Max * lambda2.min, nlambda2)
  k <- seq(1, k.min, length.out = nk) # generate proportion
  return(list(lambda2 = lambda2, lambda1 = lambda1, k=k))
}

logSeq2 = function(smax, smin, n) {
  if(smin==0){
    c(exp(seq(log(smax), log(1e-10), length.out=n-1)),0)
  }else{
    exp(seq(log(smax), log(smin), length.out=n))
  }
  
}
