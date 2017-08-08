## This file is for tunning parameter selection

## BIC

BIC_grid = function(wloss, beta, w, k, g) {
  
  l1 <- dim(beta)[1]
  l2 <-  dim(beta)[2]
  l3 <- dim(beta)[3]
  p <- dim(beta)[4]
  index <- rep(0, 3)
  BIC.max <- 1e+08
  bicPre <- BIC.max
  bicTemp <- array(0, dim=c(l1, l2, l3))
  bicTemp2<- array(0, dim=c(l1, l2, l3))
  wdf <- array(0, dim=c(l1, l2, l3))
  bdf <- array(0, dim=c(l1, l2, l3))
  n <- dim(w)[4]
  pro <- 0.5
  
  for (i in 1 : l1) {
    for (j in 1 : l2){
      for(q in 1 : l3){
        wdf[i, j, q] <- sum(w[i, j, q,]!=1)
        #bdf[i, j, q] <- getDf4kan(beta[i,j,q,],k[q],g)
        #bdf[i, j, q] <- sum(abs(beta[i,j,q,])>1e-7)
        bdf[i, j, q] <-length(unique(abs(beta[i,j,q,abs(beta[i,j,q,])>1e-7])))
        bicTemp[i, j, q] <- log(wloss[i, j, q]/(n)) + (bdf[i, j, q] + wdf[i, j, q]) * log(n)/(n)
        if(wdf[i, j, q] >= n * pro || bdf[i, j, q] + wdf[i, j, q] >= n){
          bicTemp2[i, j, q] <- BIC.max
        } else{
          bicTemp2[i, j, q]<- bicTemp[i, j, q]
          if (bicTemp[i, j, q] <= bicPre) {
            index <-c(i, j, q)
            bicPre <- bicTemp[i, j, q]
          }
        }
      }
    }
  }
  bicTemp2 <- ifelse(bicTemp2==BIC.max, max(bicTemp), bicTemp2)
  res = list(beta = beta[index[1], index[2], index[3],], 
             w = w[index[1], index[2], index[3],], 
             raw.bic = bicTemp, 
             bic=bicTemp2,
             index=index)
}

getDf4kan <- function(beta, k, g){
  e <- 0
  df <- sum(abs(beta) > 1e-7)
  cn <- 0
  for(i in 1:length(g)){
    kth <- floor(g[i] * k)
    s <- e + 1
    e <- e + g[i]
    sb <- sort(abs(beta[s:e]),decreasing = TRUE)
    if(sb[kth]!=0){
      n <- sum(sb==sb[kth])
      m <- sum(sb[1:kth]==sb[kth])
      cn <- cn + choose(n,m)
    }else{
      cn <- cn + 1
    }
  }
  return(df-cn)
}

rw <- function(x, y, nperb=50, res,...){
  d <- dim(res$beta)
  L1 <- d[1]
  L2 <- d[2]
  L3 <- d[3]
  ka <- array(0,dim=c(L1,L2,L3))
  n <- dim(x)[1]
  for(i in 1: nperb){
    v1 <- rexp(n)
    v2 <- rexp(n)
    res1 <- rkan_grid(x*v1,y*v1,...)
    res2 <- rkan_grid(x*v2,y*v2,...)
    for(l1 in 1:L1){
      for(l2 in 1: L2){
        for(l3 in 1: L3)
        {
          #compute kappa for each pair of pertubtated data
          r1 <- c(ifelse(abs(res1$beta[l1, l2, l3,])<1e-7, 0, 1), ifelse(res1$w[l1, l2, l3,]==1, 0 , 1))
          r2 <- c(ifelse(abs(res2$beta[l1, l2, l3,])<1e-7, 0, 1), ifelse(res2$w[l1, l2, l3,]==1, 0 , 1))
          ka[l1, l2, l3] <- ka[l1, l2, l3] + computeKappa(r1, r2)
        }
      }
    }
  }
  ka <- ka/nperb
  # 
  index <- which(ka==max(ka))
  index <- index[ceiling(length(index)/2)]
  l3 <- ceiling(index/(L1*L2))
  l2 <- ceiling((index-(l3-1)*(L1*L2))/L1)
  l1 <- index-(l3-1)*(L1*L2)-(l2-1)*L1
  #
  list(beta = res$beta[l1, l2, l3,], 
         w = res$w[l1, l2, l3,], index=c(l1,l2,l3), ka=ka) 

}

computeKappa <- function(r1, r2){
  k <- length(r1)
  p0 <- sum(r1==r2)/k
  pe <- sum(r1==1)*sum(r2==1)/k^2+sum(r1==0)*sum(r2==0)/k^2
  if(pe==1){
    return(0)
  } else{
    return((p0-pe)/(1-pe))
  }
   
  
}
