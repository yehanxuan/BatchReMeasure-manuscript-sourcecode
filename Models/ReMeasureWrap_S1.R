source("./Models/ReMeasure_S1.R")
## Interface ##
batch.ReMeasure.S1 = function(Yc1, Yt2, X, Zc1, Zt2, ind.r, Y.r) {
  # X: vector for case/control status, 0-control, 1-case 
  # Y: response vector 
  
 # X = c(1: nrow(Zc1))
# X = as.numeric(X %in% ind.r)
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Yc1[ind.r]
  Estimate = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, tol.c = 1e-7)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec))
}

### Example 
p = 5
beta = c(1, 1, 1, 1, 1)


S1Data = Generate_S1(nc1, nt2, nc2, p, beta, a0, a1, sigma1, sigma2, rho)

Yc1 = S1Data$Yc1
Yt2 = S1Data$Yt2
Y.r = S1Data$Yc2
Zc1 = S1Data$Zc1
Zt2 = S1Data$Zt2
Zc2 = S1Data$Zc2
Ind.r = S1Data$Index
X = c(1: nrow(Zc1))
X = as.numeric(X %in% ind.r)

batch.ReMeasure.S1(Yc1, Yt2, X, Zc1, Zt2, ind.r, Y.r)

