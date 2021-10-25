source("./Models/ReMeasure_S1.R")
## Interface ##
batch.ReMeasure.S1 = function(Y, X, Z, ind.r, Y.r) {
  # X: vector for case/control status, 0-control, 1-case 
  # Y: response vector 
  
  # X = c(1: nrow(Zc1))
  # X = as.numeric(X %in% ind.r)
  #Z <- Z[, -1, drop = FALSE]
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, tol.c = 1e-7)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
  Time = Estimate$Time
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, "p.value" = pv,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec, "Time" = Time))
}





