#
rkan <- function(x, y, lambda1, lambda2, k, g=p, beta0, w0, delta=1e-6, maxIter=50){
  
  ## initialize
  n <- length(y)
  p <- dim(x)[2]
  if(!is.null(g) & sum(g)!=p)
    stop("group info does not match x.")
  gl <- length(g)
  tp <- 2*p + gl
  L1 <- length(lambda1)
  L2 <- length(lambda2)
  L3 <- length(k)
  betaPre <- beta0
  wPre <- w0
  
  beta <- array(0,dim=c(L1,L2,L3,p))
  w <- array(0,dim=c(L1,L2,L3,n))
  iter <- array(0,dim=c(L1,L2,L3))
  
  Dmat <- matrix(0,nrow =tp, ncol=tp)
  dvec <- rep(0, tp)
  a <- diag(1,p)
  b <- matrix(0, nrow=p, ncol=gl)
  i <- 1
  j <- 0
  for(l in 1 : gl){
    j <- j + g[l]
    b[i:j, l] <- 1
    i <- j + 1
  }
  Amat <- rbind(cbind(-a, a, b),cbind(a, a, b))
  cmax <- 1000
  l <- c(rep(-cmax,p),rep(0,(tp-p)))
  u <- rep(cmax,tp)

  ## iteration
  for(l1 in 1:L1){ # for each lambda1
    for(l2 in 1:L2){ # for each lambda2
      for(l3 in 1:L3){ # for each K
        dvec[(p+1):(2*p)] <- rep(lambda1[l1],p)
        if(gl==1){
          dvec[(2*p+1):tp] <- lambda1[l1]*k[l3]
        } else{
          dvec[(2*p+1):tp] <- lambda1[l1]*k[l3]*g
        }
        while(TRUE){
          iter[l1,l2,l3] <- iter[l1,l2,l3]  + 1
          ## update x,y by wPre
          yy <- y * wPre
          xx <- diag(wPre) %*% x # can use apply
          
          ## prepare for quadprog
          Dmat[(1:p), (1:p)] <- 1/n * t(xx) %*% xx
          dvec[1:p] <- -1/n * t(xx) %*% yy
          sv <- ipop(c=dvec,H=Dmat,A=Amat,b=rep(0,2*p),l=l,u=u,r=rep(cmax,2*p),verb = 1)
          #betaTemp <- primal(sv)[1:p]
          betaTemp <- sv$primal[1:p]
          
          ## update w
          rseq <- (y - x %*% betaTemp)^2
          wTemp <- ifelse(rseq > n * lambda2[l2], n * lambda2[l2]/rseq, 1)
          
          ## check convergency
          diff <- c(betaTemp - betaPre, wTemp - wPre)
          betaPre <- as.vector(betaTemp)
          wPre <- as.vector(wTemp)
          if(sum(diff^2)/(p+n) < delta || iter[l1,l2,l3]  > maxIter){
            break
          }
        } # end inner loop
        beta[l1, l2, l3,] <- betaPre
        w[l1,l2,l3,] <- wPre
      }
    }
  }
  list(beta=beta, w=w, iter=iter, how=sv$how, sigfig=sv$sigfig)
}
