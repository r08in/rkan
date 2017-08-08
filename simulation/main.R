
# simulation
source("simulation/Simulation.R")
beta1 <- c(1, 0.5, 0, 0)
beta2 <- c(1, -0.9, -1.3, -0.5)
exam1 <- list(n=100, p=16, g=rep(4,4), beta=c(beta1,beta2,rep(0,8)))
exam <- exam1
GenerateDataFile(L=100, n=exam$n,beta=exam$beta, model=c("A","B","C","D","E"),g=exam$g, pro=0.1, da=0.8,db=0.5)

#exam1 <- list(n=100, p=400, g=rep(4,100), beta=c(beta1,beta2,rep(0,392)))
Lres_rkan <- simulation(L=100, n=exam1$n, beta=exam1$beta, g=exam1$g, model=c("A","B","C","D","E"), 
                        useDataFile = T,method="rkan")
Lres_rgkan2 <- simulation(L=20, n=exam1$n, beta=exam1$beta, g=exam1$g, model=c("A","B","C","D","E"), 
                        useDataFile = T,method="rgkan")

install.packages("pawls")
library(pawls)
Lres_pawls <- simulation(L=100, n=exam1$n, beta=exam1$beta, g=exam1$g, model=c("A","B","C","D","E"), 
                        useDataFile = T,method="pawls",
                        nlambda1=20, nlambda2=20,initial="PAWLS",lambda1.min=0.001,lambda2.min=0.05)