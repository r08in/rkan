# single case
source("simulation/GenerateData.R")
set.seed(2017)
## Low dimension
beta1 <- c(1, 0.5, 0, 0)
beta2 <- c(1, -0.9, -1.3, -0.5)
exam1 <- list(n=100, p=400, g=rep(4,100), beta=c(beta1,beta2,rep(0,392)))
exam2 <- list(n=100, p=16, g=rep(4,4), beta=c(beta1,beta2,rep(0,8)))
#exam1 <- list(n=10,g=rep(2,2), beta=c(1,2,0,0))
exam <- exam2
dout <- GenerateDataByModel(n=exam$n, beta=exam$beta, model="C", g=exam$g)
ptm <- proc.time()
res <- rkan_grid(x=dout$x, y=dout$y, lambda1=0.4, lambda2=1, k=0.5, verb=1)
as.vector(res$beta)
as.vector(res$w)
(proc.time() - ptm)[1]
res$how
res$counter
res$iter

# single case with grid of tunning parameter
source("simulation/GenerateData.R")
set.seed(2017)
beta1 <- c(1, 0.5, 0, 0)
beta2 <- c(1, -0.9, -1.3, -0.5)
exam1 <- list(n=100, p=16, g=rep(4,4), beta=c(beta1,beta2,rep(0,8)))
#exam1 <- list(n=10,g=rep(2,2), beta=c(1,2,0,0))
exam <- exam1
dout <- GenerateDataByModel(n=exam$n, beta=exam$beta, model="C", g=exam$g)
ptm <- proc.time()
res <- rkan(x=dout$x, y=dout$y, nlambda1 = 10, nlambda2 = 10, nk=8)
as.vector(res$beta)
as.vector(res$w)
(proc.time() - ptm)[1]
sum(res$how!="converged")
sum(res$iter)
