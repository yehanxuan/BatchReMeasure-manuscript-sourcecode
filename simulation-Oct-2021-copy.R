### Simulation study

rm(list = ls())
library(RConics)
source("./DataSet/Data.R")
source("./Models/ReMeasure_S1.R")
source("./Models/Batch2Only.R")
source("./Models/IgnoreBatch.R")

nc1 = 20
nt2 = 10
nc2 = 2
p = 5
beta = c(1, 1, 1, 1, 1)
a0 = 0
a1 = 0.5
sigma1 = 1
sigma2 = 2
rho = 0.5
reptim = 5000

ztest = a0est = rep(NA, reptim)
rhoest = a1est = s1est = s2est = rep(NA, reptim)
for(i in 1:reptim)
{
S1Data = Generate_S1(nc1, nt2, nc2, p, beta, a0, a1, sigma1, sigma2, rho)

Yc1 = S1Data$Yc1
Yt2 = S1Data$Yt2
Yc2 = S1Data$Yc2

Zc1 = S1Data$Zc1
Zt2 = S1Data$Zt2
Zc2 = S1Data$Zc2

Index = S1Data$Index
out = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-15)
#out = Estimate_OnlyBatch2(Zt2, Zc2, Yt2, Yc2)
#out = Estimate_Ignore_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-15)
ztest[i] = out$a0/sqrt(out$a0Var)
a0est[i] = out$a0
rhoest[i] = out$rho
#a1est[i] = out$a1
s1est[i] = out$sigma1
s2est[i] = out$sigma2
}

a = qnorm(c(0.90,0.95,0.975,0.995))
rej = NULL
for(i in 1:reptim)
{
rej = rbind(rej, abs(ztest[i]) > a)
}
apply(rej,2,mean)

hist(ztest, breaks=30, freq = FALSE)
x = seq(-4, 4, 0.01)
y = dnorm(x)
lines(x, y, lwd=2, col="red")




