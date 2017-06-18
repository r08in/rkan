source("Simulation/CombineData.R")
source("Simulation/SetupMatlab.R")
source("Simulation/GenerateData.R")
source('Simulation/mmnngreg.R')

simulation = function(L, n, beta = NULL, model = c("A", "B", "C", "D"), p = NULL, method = "PAWLS", 
    matlab = NULL, seed = 2014, useDataFile = FALSE, standardize = TRUE, penalty1 = "1-w0", updateInitial = FALSE, 
    criterion = "BIC", initCrit="BIC",intercept = TRUE, initial = "uniform", lambda2 = NULL, 
    lambda1 = NULL,lambda2.min=1e-03, lambda1.min=0.05, search = "cross", type = c("Lasso", 
        "Ridge"), pro=0.1) {
    mcount <- length(model)
    
    # define output
    iter = rep(0, L)
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
    nlambda2 <- ifelse(is.null(lambda2),50,length(lambda2))
    nlambda1 <- ifelse(is.null(lambda1),100,length(lambda1))
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
    pb <- txtProgressBar(1, mcount * L, style = 3)
    for (j in 1:mcount) {
        # for each model
        
        # initiate
        if (!is.null(seed)) {
            set.seed(seed)
        }
        p = length(beta)
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
                out = GenerateDataByModel(n = n, beta = beta, model = model[j], dataType = type, pro=pro)
            }
            
            # try different methods
            if (method == "LAD") {
                init = InitParam(out$x, out$y, method = "LAD")
                b[i, ] = init$beta
            } else if (method == "ROSS") {
              if(intercept) out$x <- AddIntercept(out$x)
              setVariable(matlab, X = out$x)
              setVariable(matlab, y = out$y)
              evaluate(matlab, "[betaRoss time]=RossSimulate(X,y)")
              times[i] <- getVariable(matlab, "time")[[1]]
              betaRoss = getVariable(matlab, "betaRoss")
              res <- list(beta=betaRoss$betaRoss)
              b[i, ] = as.vector(betaRoss$betaRoss)
            } else if (method == "ADL") {
              ptm <- proc.time()
              fit1<- ncvreg(out$x, out$y, penalty='lasso')
              beta.init <- coef(fit1)[-1 ,which.min(BIC(fit1))]
              w <- pmin(1e6, 1/abs(beta.init))
              fit2 <- ncvreg(out$x, out$y,penalty='lasso', penalty.factor=w)
              beta.alasso <- coef(fit2)[, which.min(BIC(fit2))]
              times[i] <- (proc.time() - ptm)[1]
              res <- list(beta=beta.alasso)
            } else if (method == "Lasso") {
              ptm <- proc.time()
              fit1<- ncvreg(out$x, out$y, penalty='lasso')
              
              beta.alasso <- coef(fit1)[, which.min(BIC(fit1))]
              times[i] <- (proc.time() - ptm)[1]
              res <- list(beta=beta.alasso)
            }else if (method == "LTS") {
              require(robustHD)
              ptm <- proc.time()
              res = sparseLTS(out$x, out$y, lambda=seq(from=1, to=ifelse(n>p,0.001,0.01),length.out = ifelse(n>p,100,20)), mode = "fraction",
                              intercept = TRUE)
              times[i] <- (proc.time() - ptm)[1]
              best <- res$crit$best[1]
              b[i, ] = res$coefficients[,best]
              w[i, ] = res$wt[,best]
              res$beta <- res$coefficients[,best]
            } else if (method == "IPOD") {
              ptm <- proc.time()
              H <- out$x %*% solve(t(out$x)%*%out$x)%*%t(out$x)
              res <- IPOD(out$x,out$y,H, method = "soft")
              times[i] <- (proc.time() - ptm)[1]
              w[i,] <- ifelse(res$gamma == 0, 1, 0)
              res$beta <- rep(0, p + intercept)
            }else if (method == "MMNNG") {
              ptm <- proc.time()
              try_res <- try(mmnngreg(out$x,out$y))
              times[i] <- (proc.time() - ptm)[1]
              if(class(try_res)[1]=="try-error"){
                browser()
                res$beta <- rep(0,ifelse(intercept, p+1, p))
              } else{
                res$beta <- try_res$betac
              }
            } else if (method == "MMNNG_DATA") {
                # load data file and result file
                dfile <- paste("data\\", model[j], n, "X", p,"_", pro, ".rda", sep = "")
                rfile <- paste("data\\", model[j], n, "X", p,"_", pro, "_res", ".rda", sep = "")
                lf <- try(load(dfile))
                if (class(lf) == "try-error") {
                  data <- list(out)
                  res_temp <- list(cfr=rep(0,L), cfr2=rep(0,L),ofr=rep(0,L),pdr=rep(0,L),
                                   fdr=rep(0,L), msize=rep(0,L),mses=rep(0,L),times=rep(0,L),ind=1, count=rep(0,L))
                } else {
                  if (length(data) == L) break
                  data = c(data, list(out))
                  load(rfile)
                  res_temp$ind <- res_temp$ind+1
                }
                save(data, file = dfile)
                save(res_temp, file=rfile)
                
                # fit mmnngreg
                ptm <- proc.time()
                res = mmnngreg(out$x, out$y)
                times[i] <- (proc.time() - ptm)[1]
                res$beta <- res$betac
                
            } else if (method == "PAMLS") {
              ptm <- proc.time()
              res = pamls(out$x, out$y, penalty1 = penalty1, nlambda2 = 50, nlambda1 = 100, lambda2 = lambda2,
                            lambda1=lambda1, lambda2.min=lambda2.min, lambda1.min=lambda1.min, delta = 1e-06, 
                            maxIter = 1000, initial = initial, intercept = intercept, standardize = standardize, 
                            updateInitialTimes = 0, criterion = criterion, initCrit=initCrit, search = search)
              times[i] <- (proc.time() - ptm)[1]
              b[i, ] = res$beta
              w[i, ] = 1-res$gam
              iter[i] = res$iter
              iw[i] = res$index2
              ib[i] = res$index1
              crit2[i,] = res$crit2
              crit1[i,] =  res$crit1
              lam2[i,] =  res$lambda2s
              lam1[i,] = res$lambda1s
              dfb[i] = res$bdf
              dfw[i] = res$gdf
              if(search=="grid"){
                bic[i,,] = res$res$bic
                bic2[i,,] = res$res$bic2
                bdf[i,,]= res$res$bdf
                wdf[i,,]= res$res$gdf
              }
              
            }else if (method == "PAWLS") {
              updateInitialTimes <- ifelse(updateInitial, 2, 0)
              ptm <- proc.time()
              res = pawls(out$x, out$y, nlambda2 = 50, nlambda1 = 100, lambda2 = lambda2,
                            lambda1=lambda1, delta = 1e-06, lambda1.min = lambda1.min,lambda2.min = lambda2.min,
                maxIter = 1000, initial = initial, intercept = intercept, standardize = standardize, search = search)
              times[i] <- (proc.time() - ptm)[1]
              b[i, ] = res$beta
              w[i, ] = res$w
              iter[i] = res$iter
              #crit2[i,] = res$crit2
              #crit1[i,] =  res$crit1
              lam2[i,] =  res$lambda2
              lam1[i,] = res$lambda1
              if(search=="grid"){
                bic[i,,] = res$raw.bic
                bic2[i,,] = res$bic
              }
              
            }
            
            # record result
            if(intercept) res$beta <- res$beta[-1]
            true_set <- which(beta != 0)
            pnum <- length(true_set)
            active_set <- which(res$beta != 0)
            msize[i] <- length(active_set)
            common_size <- length(intersect(true_set, active_set))
            cfr[i] <- ifelse(common_size == pnum & msize[i] == pnum, 1, 0)  # correctly fit
            ofr[i] <- ifelse(common_size == pnum & msize[i] > pnum, 1, 0)  # over fit
            cfr2[i] <- ifelse(common_size == pnum & msize[i] <= pnum + over_size, 1, 0)  # over fit by at most 2
            pdr[i] <- common_size/pnum  # positive discover rate
            fdr[i] <- ifelse(msize[i]==0,0,(msize[i] - common_size)/msize[i])  #
            mses[i] <- sum((res$beta - beta)^2)
            setTxtProgressBar(pb, (j - 1) * L + i)
            if(method == "MMNNG_DATA"){
              load(rfile)
              ind <- res_temp$ind
              res_temp$cfr[ind] <- cfr[i]
              res_temp$cfr2[ind] <- cfr2[i]
              res_temp$ofr[ind] <- ofr[i]
              res_temp$pdr[ind] <- pdr[i]
              res_temp$fdr[ind] <- fdr[i]
              res_temp$msize[ind] <- msize[i]
              res_temp$mses[ind] <- mses[i]
              res_temp$times[ind] <- times[i]
              res_temp$count[ind] <- 1
              save(res_temp, file=rfile)
            }
        }
        if(method == "MMNNG_DATA"){
          load(rfile)
          cfr <- res_temp$cfr
          cfr2 <- res_temp$cfr2
          ofr <- res_temp$ofr
          pdr <- res_temp$pdr
          fdr <- res_temp$fdr
          msize <- res_temp$msize
          mses <- res_temp$mses
          times <- res_temp$times
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
        if(method == "PAWLS" || method == "PAMLS" || method == "LTS"|| method=="IPOD" ){
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
        if(method=="PAWLS"||method=="PAMLS"){
          nres[[j]] <- list(model = model[j], CFR = CFR, CFR2 = CFR2, OFR = OFR, PDR = PDR, FDR = FDR, 
                            AN = AN, MSE = MSE, mses=mses, TIME = TIME,iter=iter,OD=OD,
                            crit2=crit2,lam2=lam2,crit1=crit1,lam1=lam1,betas=b,ws=w,
                            bic=bic,bic2=bic2
                            )
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
