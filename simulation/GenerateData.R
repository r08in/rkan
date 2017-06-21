## Generate random data from mvn for linear model

# Generate multi-dataset. each one is in linear model r paramter of variance in generating X
# matrix errorSigma variance of ramdom error offSet indicate the offset of position of non-zero
# beta
GenerateData = function(n, p = NULL, pNum = NULL, dataSetNum = 1, beta = NULL, errorSigma = NULL, 
    errorType = "n", seed = NULL, a=0.5, b=0.5, sigma.opt=c("flat","decay"), g=p) {
    # for test
    if (!is.null(seed)) {
        set.seed(seed)
    }

    require("MASS")
    if (!missing(beta)) {
        if (dataSetNum == 1) {
            p = length(beta)
            pNum = sum(beta != 0 + 0)
        } else {
            if (!is.matrix(beta) || dim(beta)[1] != dataSetNum) 
                stop("dimention of beta should match data set Num!")
            p = dim(beta)[2]
            pNum = sum(beta[1, ] != 0 + 0)
        }
    }
    
    # check data if(n<=0||p<=0||pNum<=0) stop('n or p or pNum cannot smaller than 0.')
    if (p < pNum) 
        stop("p cannot be smaller than pNum.")
    if (dataSetNum <= 0) 
        stop("dataSetNum should be positive integer.")
    if(!(0<= a & a<=1 & 0<= b & b<=1)){
      stop("a or b is not between 0 and 1.")
    }
    if(sum(g)!=p)
      stop("group info does not match x.")
    sigma.opt <- match.arg(sigma.opt)
  
  ## generate sigma
  sigma <- matrix(a*b, nrow = p, ncol = p)
  if(sigma.opt=="flat"){
    i <- 1
    j <- 0
    for(l in 1 : length(g)){
      j <- j + g[l]
      sigma[i:j, i:j] <- a
      i <- j + 1
    }
    diag(sigma)<- 1
  } else if(sigma.opt=="decay"){
    for (i in 1:p) for (j in 1:p) {
      sigma[i, j] = a^abs(i - j)
    }
  }
  
  # generate beta
  tempBeta = array(0, dim = c(dataSetNum, p))
  if (missing(beta)) {
    if (length(offSet) != dataSetNum) {
      offSet = rep(0, dataSetNum)
    }
    for (j in 1:dataSetNum) {
      tempBeta[j, ] = c(rep(0, offSet[j]), rep(1.5, pNum), rep(0, p - offSet[j] - pNum))
    }
  } else {
    tempBeta = matrix(beta, nrow = dataSetNum, ncol = p)
  }
    
    ## generate design matrix x normal parameter for design matrix
    mu = rep(0, p)
    tempx = array(0, dim = c(dataSetNum, n, p))
    for (j in 1:dataSetNum) {
        tempx[j, , ] = mvrnorm(n, mu, sigma)
    }
    
    # generate observation y
    tempy = array(0, dim = c(dataSetNum, n))
    if(is.null(errorSigma)){
      errorSigma <- sqrt(0.01 * t(beta) %*% sigma %*% beta)
    }
    for (j in 1:dataSetNum) {
        if (errorSigma == 0) {
            error = rep(0, n)
        } else if (errorType == "t") {
            error = rt(n, 2)  # df=2 for t distribution
        } else {
            error = rnorm(n, 0, errorSigma)
        }
        tempy[j, ] = tempx[j, , ] %*% tempBeta[j, ] + error
    }
    
    if (dataSetNum == 1) {
        x <- tempx[1, , ]
        beta <- tempBeta[1, ]
        y <- tempy[1, ]
    } else {
        x = tempx
        beta = tempBeta
        y = tempy
    }
    
    list(x = x, y = y, beta = beta, sigma=sigma)
    
}

# Data modification for different model
GenerateDataByModel = function(n, beta, model = c("A", "B", "C", "D", "E"), pro = 0.1, ...) {
    p = length(beta)
    if (model == "A") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, ...) 
    } else if (model == "B") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorType = "t", ...)  #errorSigma is df=2
    } else if (model == "C") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, ...)
        oNum = round(n * pro)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
    } else if (model == "D") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, ...)
        pnum = sum(beta != 0)
        oNum = round(n * pro)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
        out$x[1:oNum, (pnum + 1):(pnum + 5)] = out$x[1:oNum, (pnum + 1):(pnum + 5)] + 10
    } else if (model == "E2") {
      out = GenerateData(n = n, dataSetNum = 1, beta = beta)
      pnum = sum(beta != 0)
      ratio <- 3
      oNum = round(n * pro)
      xx <- out$x
      svd_out <- svd(xx)
      u <- svd_out$u
      v <- svd_out$v
      d <- svd_out$d
     # u[n:(n-oNum+1),] <- u[n:(n-oNum+1),] * ratio
      u[1:oNum,] <- u[1:oNum,] * ratio
      xx <- u %*% diag(d) %*% t(v)
      out$x <- xx
      u1 = runif(oNum, 0, 1)
      u2 = runif(oNum, 0, 1)
      out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
    }else if (model == "E") {
      out = GenerateData(n = n, dataSetNum = 1, beta = beta)
      pnum = sum(beta != 0)
      oNum = round(n * pro)
      xx <- out$x
      svd_out <- svd(xx)
      v <- svd_out$v
      xx[1:oNum,] <- xx[1:oNum,] + matrix(rep(5* v[,dim(v)[2]],oNum),nrow=oNum,byrow = TRUE)
      out$x <- xx
      u1 = runif(oNum, 0, 1)
      u2 = runif(oNum, 0, 1)
      out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
    }
    return(out)
}
