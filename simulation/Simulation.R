source("Simulation/GenerateData.R")


simulation = function(L, n, beta = NULL, model = c("A", "B", "C", "D"), g=p, pro=0.1, 
                      da=0.5, db=0.5, 
                      useDataFile = FALSE,seed=2017,
                      method = "rkan", standardize = TRUE, intercept = FALSE, 
                      lambda1=NULL, lambda2=NULL,
                      ...) {
    mcount <- length(model)
    
    # define output
    iw = rep(0, L)
    ib = rep(0,L)
    cfr <- rep(0, L)
    ofr <- rep(0, L)
    cfr2 <- rep(0, L)
    pdr <- rep(0, L)
    fdr <- rep(0, L)
    over_size <- 2
    msize <- rep(0, L)
    mses <- rep(0, L)
    times <- rep(0, L)
    nres <- array(list(), mcount)
    nlambda2 <- ifelse(is.null(lambda2),20,length(lambda2))
    nlambda1 <- ifelse(is.null(lambda1),20,length(lambda1))
    crit2 <- matrix(0,nrow=L,ncol=nlambda2)
    crit1 <- matrix(0,nrow=L,ncol=nlambda1)
    wdf <- array(0,dim=c(L,nlambda2,nlambda1))
    bdf <- array(0,dim=c(L,nlambda2,nlambda1))
    bic <- array(0,dim=c(L,nlambda2,nlambda1))
    bic2 <- array(0,dim=c(L,nlambda2,nlambda1))
    lam2 <- crit2
    lam1 <- crit1
    dfw <- rep(0,L)
    dfb <- rep(0,L)
    
    zero <- 1e-7
    
    pb <- txtProgressBar(1, mcount * L, style = 3)
    for (j in 1:mcount) {
        # for each model
        
        # initiate
        if (!is.null(seed)) {
            set.seed(seed)
        }
        p <- length(beta)
        if (useDataFile) {
            f = paste("data\\", model[j], n, "X", p, "_", pro, ".rda", sep = "")
            load(f)
            beta = data[[1]]$beta
            n = length(data[[1]]$y)
        }
        b = array(0, dim = c(L, ifelse(intercept,p+1,p)))
        w = array(0, dim = c(L, n))
        gam = array(0, dim = c(L, n))
        # ROC
        # simulation
        for (i in 1:L) {
            # data
            if (useDataFile) {
                out = data[[i]]
                beta <- out$beta
            } else {
                out = GenerateDataByModel(n = n, beta = beta, model = model[j], pro=pro, g=g, a=da, b=db)
            }
            
            # try different methods
            if (method == "rkan") {
              ptm <- proc.time()
              res = rkan(x=out$x, y=out$y,lambda1=lambda1, lambda2=lambda2, nlambda1=nlambda1, nlambda2=nlambda2, ...)
              times[i] <- (proc.time() - ptm)[1]
              b[i, ] = ifelse(abs(res$beta) >zero, res$beta, 0)
              w[i, ] = res$w
              lam2[i,] =  res$lambda2
              lam1[i,] = res$lambda1
              iter=res$iter
            } 
            
            # record result

            if(intercept) res$beta <- res$beta[-1]
            true_set <- which(beta != 0)
            pnum <- length(true_set)
            active_set <- which(abs(res$beta) > zero)
            msize[i] <- length(active_set)
            common_size <- length(intersect(true_set, active_set))
            cfr[i] <- ifelse(common_size == pnum & msize[i] == pnum, 1, 0)  # correctly fit
            ofr[i] <- ifelse(common_size == pnum & msize[i] > pnum, 1, 0)  # over fit
            cfr2[i] <- ifelse(common_size == pnum & msize[i] <= pnum + over_size, 1, 0)  # over fit by at most 2
            pdr[i] <- common_size/pnum  # positive discover rate
            fdr[i] <- ifelse(msize[i]==0,0,(msize[i] - common_size)/msize[i])  #
            mses[i] <- sum((res$beta - beta)^2)
            setTxtProgressBar(pb, (j - 1) * L + i)
        }
        
        # compute measurement MSE
        MSE <- round(mean(mses), 3)
        # CFR, OFR, PDR, FDR, AN
        CFR <- round(mean(cfr), 3) * 100
        OFR <- round(mean(ofr), 3) * 100
        CFR2 <- round(mean(cfr2), 3) * 100
        PDR <- round(mean(pdr), 3) * 100
        FDR <- round(mean(fdr), 3) * 100
        AN <- round(mean(msize), 3)
        TIME <- sum(times)
        # outlier dectection
        OD <- "not applicable."
        if(method == "rkan"  ){
          if(model[j] == "A" || model[j] == "B"){
            curPro <- 0
          }
          else{
            curPro <- pro
          }
          OD <- OutlierSummary(w, curPro)
          roc <- ComputeROC(w,pro=curPro)
          OD$tpr <- roc$tpr
          OD$fpr <- roc$fpr
        }
       
        # BIC curve
        if(method=="rkan"){
          nres[[j]] <- list(model = model[j], CFR = CFR, CFR2 = CFR2, OFR = OFR, PDR = PDR, FDR = FDR, 
                            AN = AN, MSE = MSE, mses=mses, TIME = TIME,iter=iter,OD=OD,
                            lam2=lam2,lam1=lam1,betas=b,ws=w)
        } else{
          nres[[j]] <- list(model = model[j], CFR = CFR, CFR2 = CFR2, OFR = OFR, PDR = PDR, FDR = FDR, 
                            AN = AN, MSE = MSE, mses=mses, TIME = TIME, iw=iw, ib=ib,OD=OD)
        }
    }
    # return
    nres
}

OutlierSummary = function(w, pro = 0.1) {
  n = dim(w)[1]
  m = dim(w)[2]
  num = round(m * pro)
  if (num == 0) {
    M = 0
    JD = 1
  } else {
    temp = apply(w[, 1:num] == 1, 1, sum)
    M = mean(temp/num)
    JD = mean(temp == 0)
  }
  S = mean(apply(w[, (num + 1):m] != 1, 1, sum)/(m - num))
  
  list(M = M, S = S, JD = JD)
  
}

ComputeROC= function(w, cutoff=seq(0,1.01,by=0.01), pro=0.1)
{
  l <- length(cutoff)
  L <- dim(w)[1]
  n <- dim(w)[2]
  ps <- NULL
  onum <- n*pro
  if(onum!=0) ps <-1:onum 
  tpr <- rep(0,l)
  fpr <- rep(0,l)
  for(j in 1 : L){
    for(i in 1 : l){
      ps_temp <- which(w[j,] < cutoff[i])
      onum_temp <- length(ps_temp)
      m <- length(intersect(ps,ps_temp))
      tpr[i] <- tpr[i] + m/ onum 
      fpr[i] <- fpr[i] + (onum_temp - m) / (n - onum)
    }
  }
  tpr <- tpr/L
  fpr <- fpr/L
  list(tpr=tpr,fpr=fpr)
}

PlotBICs=function(res,model=1,iter=0){
  if(iter==0){ # show grid iteration
    L <- dim(res[[model]]$crit2)[1]
    for(i in 1:L){
      PlotBICs(res,model,i)
      readline(prompt="Press [enter] to continue")
    }
    }else { #show paticular iteration
      res <- res[[model]]
      lam2 <- res$lam2[iter,]
      crit2 <- res$crit2[iter,]
      index2 <- res$iw[iter] 
      lam1 <- res$lam1[iter,]
      crit1 <- res$crit1[iter,]
      index1 <- res$ib[iter]
      par(mfrow=c(1,2))
      x1 <- log(lam2)
      y1 <- crit2
      plot(x1,y1,type="l",main=paste("w iter",iter,":dfw=",res$dfw[iter],"; dfb=", res$dfb[iter]-1,".", sep=""))
      abline(v=log(lam2[index2]), col="grey")
     
      x2 <- log(lam1)
      y2 <- crit1
      plot(x2,y2,type="l",main=paste("beta iter",iter,":dfw=",res$dfw[iter],"; dfb=", res$dfb[iter]-1,".", sep=""))
      abline(v=log(lam1[index1]), col="grey")
    }
}

PlotBIC2D=function(res,model=1,iter=0,tb=4,tw=5){
  if(iter==0){ # show grid iteration
    L <- dim(res[[model]]$bic)[1]
    for(i in 1:L){
      PlotBIC2D(res,model,i,tb=tb,tw=tw)
      readline(prompt="Press [enter] to continue")
    }
  }else { #show paticular iteration
    res <- res[[model]]
    crit <- res$bic[iter,,]
    crit1 <- res$bic2[iter,,]
    lam2 <- res$lam2[iter,]
    lam1 <- res$lam1[iter,]
    wdf <- res$wdf[iter,,]
    bdf <- res$bdf[iter,,]
    ib <- res$ib[iter]
    iw <- res$iw[iter]
    x11()
    par(mfrow = c(2, 3))
    image2D(crit,x=-log(lam2),y=-log(lam1),xlab="-log lambda2", ylab="-log lambda1",main="crit2")
    # x11()
    image2D(crit1,x=-log(lam2),y=-log(lam1),xlab="-log lambda2", ylab="-log lambda1",main="crit1")
    points(-log(lam2[res$iw[iter]]),-log(lam1[res$ib[iter]]),pch=24, col="black")
    # x11()
    image2D(wdf,x=-log(lam2),y=-log(lam1),xlab="-log lambda2", ylab="-log lambda1",
            main=paste("wdf=",wdf[iw,ib],sep = ""))
    points(-log(lam2[res$iw[iter]]),-log(lam1[res$ib[iter]]),pch=24, col="black")
    # x11()
    image2D(bdf,x=-log(lam2),y=-log(lam1),xlab="-log lambda2", ylab="-log lambda1",
            main=paste("bdf=",bdf[iw,ib],sep = ""))
    points(-log(lam2[res$iw[iter]]),-log(lam1[res$ib[iter]]),pch=24, col="black")
    # x11()
    image2D(bdf+wdf,x=-log(lam2),y=-log(lam1),xlab="-log lambda2", ylab="-log lambda1",main="w+bdf")
    points(-log(lam2[res$iw[iter]]),-log(lam1[res$ib[iter]]),pch=24, col="black")
    #last one
    n <- dim(bdf)[1]
    m <- dim(bdf)[2]
    df <- matrix(0, nrow=n, ncol=m)
    for(i in 1:n){
      for(j in 1:m){
        if(bdf[i,j]==tb & wdf[i,j]==tw){
          df[i,j]=1
        }else if( (bdf[i,j]>tb & bdf[i,j]<=tb*1.5) & (wdf[i,j]>tw & wdf[i,j]<= tw*1.5)){
          df[i,j]=0.7
        }else if( (bdf[i,j]>tb*1.5 & bdf[i,j]<=tb*2) & (wdf[i,j]>tw*1.5 & wdf[i,j]<= tw*2)){
          df[i,j]=0.4
        }else if( (bdf[i,j]>tb*2 & bdf[i,j]<=tb*3) & (wdf[i,j]>tw*2 & wdf[i,j]<= tw*3)){
          df[i,j]=0.2
        }else{
          df[i,j]=0
        }
      }
    }
    image2D(df,x=-log(lam2),y=-log(lam1),xlab="-log lambda2", ylab="-log lambda1",main="True df")
    points(-log(lam2[res$iw[iter]]),-log(lam1[res$ib[iter]]),pch=24, col="black")
  }
}

PlotEfficiency <- function(lambda,weight,x, sigma=2,main=NULL)
{
  L <- length(lambda)
  n <- dim(x)[1]
  p <- dim(x)[2]
  varb <- matrix(0,nrow=L,ncol=p)
  for(i in 1 : L){
    varb[i,] <- diag(solve(t(x) %*% diag(weight[i,]^2) %*% x))  * sigma^2
  }
  op <- apply(weight!=1,1,sum)/n
  #plot
  plot(op,varb[,p],type="n",xlim=c(0,0.6), 
       xlab="Percent of Outliers", main=main, ylab=expression(paste("Var(",hat(beta),")",sep = "")))
  qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col <- sample(col_vector, p)
  for(i in 1:p){
    lines(op,varb[,i],lwd=2, lty=1,col=col[i])
  }
}
