#
rkan <- function(x, y, lambda1, lambda2, k, beta0, w0, delta, maxIter){
  
  ## initialize
  n <- length(y)
  p <- dim(x)[2]
  L1 <- length(lambda1)
  L2 <- length(lambda2)
  L3 <- length(k)
  betaPre <- beta0
  wPre <- w0
  
  beta <- array(0,dim=c(L1,L2,L3,p))
  w <- array(0,dim=c(L1,L2,L3,n))
  
  Dmat <- matrix(0,nrow =(2*p+1), ncol=(2*p+1))
  dvec <- rep(0, 2*p+1)
  a <- diag(1,p)
  b <- rep(1,p)
  Amat <- rbind(cbind(-a, a, b),cbind(a, a, b), cbind(0*a, a, 0*b))
  
  ## iteration
  for(l1 in 1:L1){ # for each lambda1
    for(l2 in 1:L2){ # for each lambda2
      for(l3 in 1:L3){ # for each K
        dvec[(p+1):(2*p+1)] <- c(-lambda1[l1]*rep(1,p),k[l3])
        iter <- 0
        while(TRUE){
          iter <- iter + 1
          ## update x,y by wPre
          yy <- y * wPre
          xx <- diag(wPre) %*% x # can use apply
          
          ## prepare for quadprog
          Dmat[(1:p), (1:p)] <- 1/n * t(xx) %*% xx
          dvec[1:p] <- 1/n * t(xx) %*% yy
          sol <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat)
          betaTemp <- sol$solution
          
          ## update w
          rseq <- (y - t(x) %*% betaTemp)^2
          wTemp <- ifelse(rseq > n * lambda2[l2], n * lambda2[l2]/rseq, 1)
          
          ## check convergency
          diff <- c(betaTemp - betaPre, wTemp - wPre)
          betaPre <- betaTemp
          wPre <- wTemp
          if(sum(diff^2)/(p+n) < delta || iter > maxIter){
            break
          }
        } # end inner loop
        beta[l1, l2, l3,] <- betaPre
        w[l1,l2,l3,] <- wPre
      }
    }
  }
}
