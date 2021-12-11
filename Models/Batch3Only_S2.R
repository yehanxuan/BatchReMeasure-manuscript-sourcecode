##### Batch 3 only S2 ######
Variance_a0_OnlyBatch3_S2 = function(Zc3, Zt3, Yc3, Yt3, sigma3H) {
  nt3 = nrow(Zt3)
  nc3 = nrow(Zc3)
  Zc3_c =  t( t(Zc3) - colMeans(Zc3) )
  Yc3_c = Yc3 - mean(Yc3)
  Zt3_c = t( t(Zt3) - colMeans(Zt3) )
  Yt3_c = Yt3 - mean(Yt3)
  
  Cov =  (t(Zt3_c)%*%Zt3_c + t(Zc3_c)%*%Zc3_c) + 0.001*diag(ncol(Zc3))
  diff = colMeans(Zt3) - colMeans(Zc3)
  c1 = t(rep(1, nt3))/nt3 - t(diff)%*%solve(Cov, t(Zt3_c))
  c2 = -t(rep(1, nc3))/nc3 - t(diff)%*%solve(Cov, t(Zc3_c))
  a0Var = (sigma3H^2)*( c1%*%t(c1) + c2%*%t(c2) )
  return(a0Var)
}


Estimate_OnlyBatch3_S2 = function(Zc3, Zt3, Yc3, Yt3) {
  start = proc.time()[1]
  nt3 = nrow(Zt3)
  nc3 = nrow(Zc3)
  Zc3_c =  t( t(Zc3) - colMeans(Zc3) )
  Yc3_c = Yc3 - mean(Yc3)
  Zt3_c = t( t(Zt3) - colMeans(Zt3) )
  Yt3_c = Yt3 - mean(Yt3) 
  
  Cov =  (t(Zt3_c)%*%Zt3_c + t(Zc3_c)%*%Zc3_c) + 0.001*diag(ncol(Zc3))
  Cor = (t(Zt3_c)%*%Yt3_c + t(Zc3_c)%*%Yc3_c)
  betaH = solve(Cov, Cor)
  a3H = mean(Yc3 - Zc3%*%betaH)
  a0H =  mean(Yt3 - Zt3%*%betaH) - a3H
  
  tmp = sum( (Yc3_c - Zc3_c%*%betaH)^2) + sum((Yt3_c  - Zt3_c%*%betaH)^2)
  sigma3HSq = tmp/(nc3 + nt3)
  sigma3H = sqrt(sigma3HSq)
  
  a0Var = Variance_a0_OnlyBatch3_S2(Zc3, Zt3, Yc3, Yt3, sigma3H)
  Time = proc.time()[1] - start
  
  return(list("a0" = a0H, "a0Var" = a0Var, "a3" = a3H, "beta" = betaH, "sigma3" = sigma3H,
              "Time" = Time))
}


oneReplicate_OnlyBatch3 = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S2.R")
  
  Estimate = batch.Batch3.S2(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  a3H = Estimate$a3
  betaH = Estimate$beta
  sigma1H = NULL
  sigma2H = NULL
  sigma3H = Estimate$sigma3
  rho1H = NULL
  rho2H = NULL
  objVec = Estimate$objVec
  Time = Estimate$Time 
  pv = Estimate$p.value
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H,
              "rho1" = rho1H, "rho2" = rho2H, "beta" = betaH, "objVec" = objVec, 
              "Time" = Time, "p.value" = pv))
}
 
oneReplicateWrap_OnlyBatch3 = function(seedJ) {
  eval = oneReplicate_OnlyBatch3(seedJ)
  return(eval)
} 


batch.Batch3.S2 = function(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  
  Zc3 = Zc1[ind.r1, , drop = F]
  Yc3 = Y.r1
  Zt3 = Zt2[ind.r2, , drop = F]
  Yt3 = Y.r2 
  
  Estimate = Estimate_OnlyBatch3_S2(Zc3, Zt3, Yc3, Yt3)
  a0H = Estimate$a0 
  a0Var = Estimate$a0Var 
  a1H = NULL
  a3H = Estimate$a3 
  betaH = Estimate$beta 
  rho1H = NULL
  rho2H = NULL
  sigma1H = NULL
  sigma2H = NULL
  sigma3H = Estimate$sigma3 
  objVec = NULL
  pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var)))
  Time = Estimate$Time
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H,
              "beta" = betaH, "rho1" = rho1H, "rho2" = rho2H, "p.value" = pv,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H, 
              "objVec" = objVec, "Time" = Time))
}
