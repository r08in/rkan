#' @export
rkan_grid <- function(x, y, lambda1, lambda2, k, g=p, delta=1e-6, maxIter=50, start.b=NULL, start.w=NULL,...){
  
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
  betaPre <- rep(0, p)
  if(is.null(start.w)){
    wPre <- rep(1, n)
  } else {
    wPre <- start.w
  }
  
  beta <- array(0,dim=c(L1, L2, L3, p))
  w <- array(0,dim=c(L1, L2, L3, n))
  iter <- array(0,dim=c(L1, L2, L3))
  wloss <- array(0,dim=c(L1, L2, L3))
  how <- array(0,dim=c(L1, L2, L3))
  counter <- array(0,dim=c(L1, L2, L3))
  
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
        dvec[(2*p+1):tp] <- lambda1[l1]*floor(k[l3]*g)
        count <- 0
        while(TRUE){
          iter[l1,l2,l3] <- iter[l1,l2,l3]  + 1
          ## update x,y by wPre
          yy <- y * wPre
          xx <- diag(wPre) %*% x # can use apply
          
          ## prepare for quadprog
          Dmat[(1:p), (1:p)] <- 1/n * t(xx) %*% xx
          dvec[1:p] <- -1/n * t(xx) %*% yy
          sv <- ipop(c=dvec,H=Dmat,A=Amat,b=rep(0,2*p),l=l,u=u,r=rep(cmax,2*p),start=start.b,...)
          count <- max(count, sv$counter) 
          #start.b <- sv$sol
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
        wloss[l1, l2, l3] <- t(wPre^2) %*% rseq
        how[l1, l2, l3] <- sv$how
        counter[l1, l2, l3] <- count
      } # end loop for k
      betaPre <- beta[l1, l2, 1, ]
      wPre <- w[l1, l2, 1, ]
    } # end loop for lambda2
    betaPre <- beta[l1, 1, 1, ]
    wPre <- w[l1, 1, 1, ]
  } # end loop for lambda1
  list(beta=beta, w=w, wloss=wloss, iter=iter, how=how, counter=counter)
}

#' @export
rkan_grid2 <- function(x, y, lambda1, lambda2, k, g=p, delta=1e-6, maxIter=50, start.b=rep(0,p), start.w=rep(1,n),...){
  
  ## initialize
  n <- length(y)
  p <- dim(x)[2]
  L1 <- length(lambda1)
  L2 <- length(lambda2)
  L3 <- length(k)
  if(!is.null(g) & sum(g)!=p)
    stop("group info does not match x.")
  beta0 <- rep(1, p)
  w0 <- rep(0, n)
  res <- .Call("rkan_GRID", x, y, lambda1, lambda2, k, g, beta0, w0, delta, maxIter, 
               #ifelse(intercept, 1, 0), 
               startBeta = start.b, startW = start.w)
  res = list(beta = array(res[[1]], dim = c(L1,L2,L3, p)), w = array(res[[2]], dim = c(L1,L2,L3, n)), 
             wloss = array(res[[3]], dim = c(L1,L2,L3)), loss = array(res[[4]], dim = c(L1,L2,L3)), 
             iter = array(res[[5]], dim = c(L1,L2,L3)))
  
}

