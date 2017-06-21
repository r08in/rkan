## Generate random data from mvn for linear model

# Generate multi-dataset. each one is in linear model r paramter of variance in generating X
# matrix errorSigma variance of ramdom error offSet indicate the offset of position of non-zero
# beta
GenerateData = function(n, p = NULL, pNum = NULL, dataSetNum = 1, beta = NULL, r = 0.9, errorSigma = 1, 
    errorType = "n", seed = NULL) {
    # for test
    if (!is.null(seed)) {
        set.seed(seed)
    }
    # 
    
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
    
    ## generate design matrix x normal parameter for design matrix
    mu = rep(0, p)
    sigma = matrix(0, nrow = p, ncol = p)
    
    ## generate sigma
    for (i in 1:p) for (j in 1:p) {
      sigma[i, j] = r^abs(i - j)
    }
    
    tempx = array(0, dim = c(dataSetNum, n, p))
    for (j in 1:dataSetNum) {
        tempx[j, , ] = mvrnorm(n, mu, sigma)
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
    
    # generate observation y
    tempy = array(0, dim = c(dataSetNum, n))
    for (j in 1:dataSetNum) {
        if (errorSigma == 0) {
            error = rep(0, n)
        } else if (errorType == "t") {
            error = rt(n, 2)  # df=2 for t distribution
        } else {
            error = errorSigma * rnorm(n, 0, 1)
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
    
    list(x = x, y = y, beta = beta)
    
}

# Data modification for different model
GenerateDataByModel = function(n, beta, errorSigma = 2, r = 0.5, model = c("A", "B", "C", "D", "E"), pro = 0.1) {
    p = length(beta)
    if (model == "A") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r)  #errorSigma=2
    } else if (model == "B") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, errorType = "t", 
            r = r)  #errorSigma is df=2
    } else if (model == "C") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r)
        oNum = round(n * pro)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
    } else if (model == "D") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r)
        pnum = sum(beta != 0)
        oNum = round(n * pro)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
        out$x[1:oNum, (pnum + 1):(pnum + 5)] = out$x[1:oNum, (pnum + 1):(pnum + 5)] + 10
    } else if (model == "E2") {
      out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r)
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
      out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r)
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
      #out$y[1:oNum] = xx[1:oNum,]%*% beta + errorSigma * rnorm(oNum, 0, 1) + 
        #ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
    }
    return(out)
}
