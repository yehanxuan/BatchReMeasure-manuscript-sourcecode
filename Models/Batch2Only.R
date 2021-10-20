### Only consider batch 2 ####
Variance_a0_OnlyBatch2 = function(Zt2, Zc2, Yt2, Yc2, sigma2H) {
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  Zc2_c =  t( t(Zc2) - colMeans(Zc2) )
  Yc2_c = Yc2 - mean(Yc2)
  Zt2_c = t( t(Zt2) - colMeans(Zt2) )
  Yt2_c = Yt2 - mean(Yt2)
  
  Cov =  (t(Zt2_c)%*%Zt2_c + t(Zc2_c)%*%Zc2_c)
  
  diff = colMeans(Zt2) - colMeans(Zc2)
  c1 = t(rep(1, nt2))/nt2 - t(diff)%*%solve(Cov, t(Zt2_c))
  c2 = -t(rep(1, nc2))/nc2 - t(diff)%*%solve(Cov, t(Zc2_c))
  a0Var = (sigma2H^2)*( c1%*%t(c1) + c2%*%t(c2) )
  return(a0Var)
}


Estimate_OnlyBatch2 = function(Zt2, Zc2, Yt2, Yc2) {
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  Zc2_c =  t( t(Zc2) - colMeans(Zc2) )
  Yc2_c = Yc2 - mean(Yc2)
  Zt2_c = t( t(Zt2) - colMeans(Zt2) )
  Yt2_c = Yt2 - mean(Yt2)
  
  Cov =  (t(Zt2_c)%*%Zt2_c + t(Zc2_c)%*%Zc2_c)
  Cor = (t(Zt2_c)%*%Yt2_c + t(Zc2_c)%*%Yc2_c)
  betaH = solve(Cov, Cor)
  a1H = mean(Yc2 - Zc2%*%betaH)
  a0H = mean(Yt2 - Zt2%*%betaH) - mean(Yc2 - Zc2%*%betaH)
  
  tmp  = sum( (Yt2_c - Zt2_c%*%betaH)^2) + sum( (Yc2_c - Zc2_c%*%betaH)^2)
  sigma2HSq = tmp/(nc2 + nt2)
  sigma2H = sqrt(sigma2HSq)
  
  a0Var = Variance_a0_OnlyBatch2(Zt2, Zc2, Yt2, Yc2, sigma2H)
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "sigma2" = sigma2H))
}


oneReplicate_OnlyBatch2 = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-S1.R")
  Estimate = Estimate_OnlyBatch2(Zt2, Zc2, Yt2, Yc2)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  sigma1H = NULL
  sigma2H = Estimate$sigma2
  rhoH = NULL
  objVec = NULL
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "rho" =  rhoH, "beta" = betaH, "objVec" = objVec))
}


oneReplicateWrap_OnlyBathc2 = function(seedJ) {
  eval = oneReplicate_OnlyBatch2(seedJ)
  return(eval)
}











