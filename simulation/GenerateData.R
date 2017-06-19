## Generate random data from mvn for linear model

# Generate multi-dataset. each one is in linear model r paramter of variance in generating X
# matrix errorSigma variance of ramdom error offSet indicate the offset of position of non-zero
# beta
GenerateData = function(n, p = NULL, pNum = NULL, dataSetNum = 1, beta = NULL, r = 0.9, errorSigma = 1, 
    errorType = "n", offSet = 0, outlier.op = "NONE", outlier.pro = 0.1, outlier.r = 10, dataType = c("Lasso", 
        "Ridge"), seed = NULL) {
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
    if (outlier.op != "NONE" && (outlier.pro > 1 || outlier.pro < 0)) 
        stop("the proportion of outlier is illegal!")
    dataType <- match.arg(dataType)
    
    ## generate design matrix x normal parameter for design matrix
    mu = rep(0, p)
    sigma = matrix(0, nrow = p, ncol = p)
    xx = matrix(0, nrow = n, ncol = p)
    
    if (dataType == "Ridge") {
        for (i in 1:p) for (j in 1:p) {
            if (i == j) {
                sigma[i, j] = 1
            } else {
                sigma[i, j] = r
            }
            
        }
    } else {
        for (i in 1:p) for (j in 1:p) {
            sigma[i, j] = r^abs(i - j)
        }
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
            if (dataType == "Ridge") 
                error = 1/sqrt(3) * rt(n, 3)
        } else {
            error = errorSigma * rnorm(n, 0, 1)
        }
        if (outlier.op != "NONE") {
            # generate outliers
            oNum = round(n * outlier.pro)
            if (outlier.op == "MEANSHIFT") {
                shift = c(rep(outlier.r, oNum), rep(0, n - oNum))
                error = error + shift
            }
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

## Generate X and Y seperately(without beta)
GenerateDataSep = function(n = 150, p = 50, k = 6, r = 0.3, delta = 5, beta = NULL) {
    # check data
    if (n <= 0 || p <= 0 || k <= 0) 
        stop("n or p or k cannot smaller than 0.")
    if (n <= 3 * k) 
        stop("3K cannot be larger or equal to n.")
    
    # generate Y
    mu = rep(0, k)
    sigma = diag(k)
    require("MASS")
    l = mvrnorm(n, mu, sigma)
    errorSigma = sqrt(k)/3
    yerror = rnorm(n, 0, errorSigma)
    if (is.null(beta)) {
        y = apply(l, 1, sum) + yerror
    } else {
        y = l %*% beta[1:k] + yerror
    }
    
    # generate X
    xerror = mvrnorm(n, rep(0, p), diag(p))
    x = matrix(nrow = n, ncol = p, 0)
    range = 1:k
    x[, range] = l + r * xerror[, range]
    range = (k + 1):(2 * k)
    x[, range] = l + delta * xerror[, range]
    range = (2 * k + 1):(3 * k)
    x[, range] = l + delta * xerror[, range]
    range = (3 * k + 1):p
    x[, range] = xerror[, range]
    
    # return
    list(x = x, y = y, beta = c(rep(1, k), rep(0, p - k)))
}


## Generate group data
GenerateGroupData = function(groupSize, groupNum, validGroupNum, dataSize, offSet = 0) {
    n = dataSize
    p = groupNum
    pNum = validGroupNum
    dataSetNum = groupSize
    data = GenerateData(n, p, pNum, dataSetNum, errorSigma = 0.001, offSet = rep(offSet, dataSetNum))
    out = CombineMultiLm(data$x, data$y)
}


## Generate linear model data The design matrix is composed of dummy data
GenerateDummyData = function(n, groupInfo, validGroupNum, errorSigma = 1, offSet = 0) {
    # check data
    p = length(groupInfo)
    pNum = validGroupNum
    if (p < pNum || pNum <= 0) 
        stop("The number of valid group is illegal!")
    if (offSet + pNum > p) 
        stop("The offSet is illegal!")
    # for test
    set.seed(120)
    
    # generate design matrix x
    m = sum(groupInfo)
    x = matrix(0, n, m)
    left = 0
    for (j in 1:p) {
        for (i in 1:n) {
            index = sample(0:groupInfo[j], 1)
            if (index == 0) 
                next
            x[i, (left + index)] = 1
        }
        left = left + groupInfo[j]
    }
    
    # generate beta
    beta = rep(0, m)
    if (offSet == 0) {
        start = 1
        end = sum(groupInfo[1:pNum])
    } else {
        start = sum(groupInfo[1:offSet]) + 1
        end = sum(groupInfo[1:(offSet + pNum)])
    }
    beta[start:end] = 2
    
    # generate y
    
    # x[,1]=0
    if (errorSigma == 0) {
        error = rep(0, n)
    } else {
        error = rnorm(n, 0, errorSigma)
    }
    y = x %*% beta + error
    # return
    list(x = x, y = y, beta = beta)
}

## Generate Dummy model
GenerateDummyModel = function(sizeInfo, groupInfo, validGroupNumInfo, offSet = 0, errorSigma) {
    ## intial
    m = length(sizeInfo)  #dataset Num
    p = sum(groupInfo)
    
    ## generate design matrix
    X = matrix(0, 0, p)
    Y = NULL
    if (length(offSet) != m) {
        offSet = rep(0, m)
    }
    for (i in 1:m) {
        # generate group dummy data
        # GenerateDummyData=function(n,groupInfo,validGroupNum,errorSigma=1,offSet=0)
        out = GenerateDummyData(sizeInfo[i], groupInfo, validGroupNumInfo[i], offSet = offSet[i], 
            errorSigma = errorSigma)
        X = rbind(X, out$x)
        Y = c(Y, out$y)
    }
    
    # combine each dataSet
    CombineDataset(X, sizeInfo, groupInfo, Y)
}

# Data modification for different model
GenerateDataByModel = function(n, beta, errorSigma = 2, r = 0.5, model = c("A", "B", "C", "D", "E"), 
    dataType = c("Lasso", "Ridge"), pro = 0.1) {
    p = length(beta)
    if (model == "A") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)  #errorSigma=2
    } else if (model == "B") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, errorType = "t", 
            r = r, dataType = dataType)  #errorSigma is df=2
    } else if (model == "C") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
        oNum = round(n * pro)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
    } else if (model == "D") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
        pnum = sum(beta != 0)
        oNum = round(n * pro)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
        out$x[1:oNum, (pnum + 1):(pnum + 5)] = out$x[1:oNum, (pnum + 1):(pnum + 5)] + 10
    } else if (model == "E2") {
      out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
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
      out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
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
    }else if (model == "E3") {
      out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
      pnum = sum(beta != 0)
      oNum = round(n * pro)
      xx <- out$x
      svd_out <- svd(xx)
      v <- svd_out$v
      xx[1:oNum,] <- xx[1:oNum,] + matrix(rep(5* v[,p],oNum),nrow=oNum,byrow = TRUE)
      out$x <- xx
      u1 = runif(oNum, 0, 1)
      u2 = runif(oNum, 0, 1)
      #out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
      out$y[1:oNum] = xx[1:oNum,]%*% beta + errorSigma * rnorm(oNum, 0, 1) + 
        ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
    }else if (model == "C2") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
        oNum = round(n * 0.1)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 1, -1, 1) * (20 + 10 * u2)
        ## random change sign method
        u = runif(n, 0, 1)
        index = u < 0.5
        out$y[index] = -out$y[index]
        out$x[index, ] = -out$x[index, ]
        set.seed(round(out$x[1, 1] * 100))
        
        
    } else if (model == "C3") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
        oNum = round(n * 0.1)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 1, -1, 1) * (20 + 10 * u2)
        ## random change sign method
        set.seed(round(out$x[1, 1] * 100))
        
    } else if (model == "LA") {
        # nocontamination
        p = length(beta)
        k = sum(beta != 0)
        out = GenerateDataSep(n = n, p = p, k = k)
    } else if (model == "LB") {
        # verticle outliers
        out = GenerateDataSep(n = n, p = p)
        oNum = round(n * 0.1)
        out$y[1:oNum] = out$y[1:oNum] + 20
    } else if (model == "LC") {
        out = GenerateDataSep(n = n, p = p, beta = beta)
        oNum = round(n * 0.1)
        out$y[1:oNum] = out$y[1:oNum] + 20
        mu = rep(50, p)
        sigma = diag(p)
        require("MASS")
        out$x[1:oNum, ] = mvrnorm(oNum, mu, sigma)
    } else if (model == "LD") {
        out = GenerateDataSep(n = n, p = p)
        oNum = round(n * 0.1)
        # x leverage
        mu = rep(10, p)
        sigma = diag(0.01, p)
        require("MASS")
        out$x[1:oNum, ] = mvrnorm(oNum, mu, sigma)
        # y outlier
        gam = rep(-1/p, p)
        g = 2 * (1:oNum)
        out$y[1:oNum] = out$x[1:oNum, ] %*% gam * g
    } else if (model == "HA") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = 0.5, r = 0.5)
    } else if (model == "HB") {
        # vericle outliers
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r)
        oNum = round(n * 0.1)
        out$y[1:oNum] = out$y[1:oNum] + 20
    } else if (model == "HC") {
        # vericle outliers+leverage
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r)
        oNum = round(n * 0.1)
        out$y[1:oNum] = out$y[1:oNum] + 20
        mu = rep(50, p)
        sigma = diag(p)
        require("MASS")
        out$x[1:oNum, ] = mvrnorm(oNum, mu, sigma)
        
    } else if (model == "HD") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r)
        oNum = round(n * 0.1)
        # x leverage
        mu = rep(10, p)
        sigma = diag(0.01, p)
        require("MASS")
        out$x[1:oNum, ] = mvrnorm(oNum, mu, sigma)
        # y outlier
        gama = rep(-1/p, p)
        g = 2 * (1:oNum)
        out$y[1:oNum] = out$x[1:oNum, ] %*% gamma * g
    } else if (model == "RA") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)  #errorSigma=2
    } else if (model == "RB") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, errorType = "t", 
            r = r, dataType = dataType)  # is df=2
    } else if (model == "RC") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
        oNum = round(n * pro)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
    } else if (model == "RD") {
        out = GenerateData(n = n, dataSetNum = 1, beta = beta, errorSigma = errorSigma, r = r, dataType = dataType)
        oNum = round(n * pro)
        u1 = runif(oNum, 0, 1)
        u2 = runif(oNum, 0, 1)
        out$y[1:oNum] = out$y[1:oNum] + ifelse(u1 < 0.5, -1, 1) * (20 + 10 * u2)
        pIndex = (beta == 0.1)
        # pIndex=(round(p*(1-pro))+1):p
        out$x[1:oNum, pIndex] = out$x[1:oNum, pIndex] + 10
    }
    return(out)
}
