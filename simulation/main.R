
# simulation
source("simulation/Simulation.R")
beta1 <- c(1, 0.5, 0, 0)
beta2 <- c(1, -0.9, -1.3, -0.5)
#exam1 <- list(n=100, p=16, g=rep(4,4), beta=c(beta1,beta2,rep(0,8)))
exam1 <- list(n=100, p=400, g=rep(4,100), beta=c(beta1,beta2,rep(0,392)))
Hres_rkan <- simulation(L=100, n=exam1$n, beta=exam1$beta, g=exam1$g, model=c("A"), method="rkan", nk=8)

Lres_rkan10 <- Lres_rkan 
