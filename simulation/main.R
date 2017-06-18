#main
install.packages("quadprog")
require(quadprog)
require(kernlab)
source("simulation/GenerateData.R")
install.packages("Rcplex")
library(Rcplex)

set.seed(2017)
n <- 50
p <- 8
beta <- c(3, 2, 1, 0, 0, 0, 0, 0)
dout <- GenerateDataByModel(n=n, beta=beta, model="C")
res <- rkan(x=dout$x, y=dout$y, lambda1=0.1, lambda2=1, k=2, beta0=rep(0,p),w0=rep(1,n))


m <- 500
set <- sample(1:dim(spam)[1],m)
x <- scale(as.matrix(spam[,-58]))[set,]
y <- as.integer(spam[set,58])
y[y==2] <- -1
##set C parameter and kernel
C <- 5
rbf <- rbfdot(sigma = 0.1)
## create H matrix etc.
H <- kernelPol(rbf,x,,y)
c <- matrix(rep(-1,m))
A <- t(y)
b <- 10
l <- matrix(rep(0,m))
u <- matrix(rep(C,m))
r <- 1e5
sv <- ipop(c,H,A,b,l,u,r)
sv