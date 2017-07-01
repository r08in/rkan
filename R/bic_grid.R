## This file is for tunning parameter selection

## BIC

BIC_grid = function(wloss, beta, w, k=NULL, g=p) {
  
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
      for(k in 1 : l3){
        wdf[i, j, k] <- sum(w[i, j, k,]!=1)
        bdf[i, j, k] <- sum(abs(beta[i, j, k,]) > 1e-7)
        bicTemp[i, j, k] <- log(wloss[i, j, k]/(n)) + (bdf[i, j, k] + wdf[i, j, k]) * log(n)/(n)
        if(wdf[i, j, k] >= n * pro || bdf[i, j, k] + wdf[i, j, k] >= n){
          bicTemp2[i, j, k] <- BIC.max
        } else{
          bicTemp2[i, j, k]<- bicTemp[i, j, k]
          if (bicTemp[i, j, k] <= bicPre) {
            index <-c(i, j, k)
            bicPre <- bicTemp[i, j, k]
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

getDf4kan <- function(beta, kg){
  
}
