# single case
# install.packages("kernlab")
# require(kernlab)

source("simulation/GenerateData.R")
set.seed(2023)
## Low dimension
beta1 <- c(1, 0.5, 0, 0)
beta2 <- c(1, -0.9, -1.3, -0.5)
#beta1 <- c(c(1, 0.5, 0, 0),rep(0,6))
#beta2 <- c(c(1, 0.9, 1.3, 0.5,0.5,0.5,0.5,0.5),rep(0,2))
#exam1 <- list(n=100, p=400, g=rep(4,100), beta=c(beta1,beta2,rep(0,392)))
exam2 <- list(n=100, p=16, g=rep(4,4), beta=c(beta1,beta2,rep(0,8)))
#exam3 <- list(n=100,p=40,g=rep(10,4), beta=c(beta1,beta2,rep(0,20)))
#exam1 <- list(n=10,g=rep(2,2), beta=c(1,2,0,0))
exam <- exam2
dout <- GenerateDataByModel(n=exam$n, beta=exam$beta, model="A", g=exam$g,a=0.8,b=0.5)
ptm <- proc.time()
#k <- 15/16
beta0 <- rep(1,exam$p)
w0 <- rep(0,exam$n)
res <- rkan_grid(x=dout$x, y=dout$y,lambda1=0.1, lambda2=2,k=seq(16/16,1/16,length.out = 16),beta0=beta0,w0=w0,delta=1e-20)
beta <- matrix(res$beta,nrow=16,byrow=F)
View(beta)
w <- matrix(res$w,nrow=16,byrow=F)
View(w)
ifelse(abs(beta)<1e-7,0,beta)


#sum(sort(abs(beta),decreasing = TRUE)[1:(k)])r
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
dout <- GenerateDataByModel(n=exam$n, beta=exam$beta, model="A", g=exam$g)
ptm <- proc.time()
res <- rkan(x=dout$x, y=dout$y, nlambda1 = 20, lambda2 = 100, nk=16, crit = "BIC")
View(matrix(res$raw.bic,nrow=16,byrow=T))
as.vector(res$beta)
as.vector(res$w)
(proc.time() - ptm)[1]
sum(res$how!="converged")
sum(res$iter)
