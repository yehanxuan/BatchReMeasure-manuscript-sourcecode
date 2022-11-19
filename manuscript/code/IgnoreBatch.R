### If we ignore the batch effect ####

Ignore_Update_sigma1 = function(Zc1, Yc1, betaH) {
  nc1 = nrow(Zc1)
  sigma1H = sqrt( mean( (Yc1 - Zc1%*%betaH)^2 ) ) 
  return(sigma1H)
}


Ignore_Update_sigma2 = function(Zt2, Yt2, betaH){
  nt2 = nrow(Zt2)
  Zt2c = t( t(Zt2) - colMeans(Zt2) ) 
  Yt2c = Yt2 - mean(Yt2)
  sigma2H = sqrt( mean( (Yt2c - Zt2c%*%betaH)^2 ) )
  return(sigma2H)
}


Ignore_Update_beta = function(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H) {
  nt2 = nrow(Zt2)
  Zt2c = t( t(Zt2) - colMeans(Zt2) ) 
  Yt2c = Yt2 - mean(Yt2)
  Cov = t(Zc1)%*%Zc1/sigma1H^2 + t(Zt2c)%*%Zt2c/sigma2H^2 
  Cor = t(Zc1)%*%Yc1/sigma1H^2 + t(Zt2c)%*% Yt2c/sigma2H^2
  betaH = solve(Cov, Cor)
  return(betaH)            
}


Ignore_Update_a0 = function(Zt2, Yt2, betaH) {
  a0H = mean(Yt2 - Zt2%*%betaH)
  return(a0H)
}


Ignore_Objective_S1 = function(Zc1, Zt2, Yc1, Yt2, a0H, betaH, sigma1H, sigma2H){
  nc1 = length(Yc1)
  nt2 = nrow(Zt2)
  mu1 = Zc1 %*% betaH
  mu2 = a0H + Zt2%*% betaH
  Yc1Scale = (Yc1 - mu1)/sigma1H
  Yt2Scale = (Yt2 - mu2)/sigma2H
  
  obj = as.numeric( nc1*log(sigma1H) + t(Yc1Scale)%*%Yc1Scale/2 +
    nt2*log(sigma2H) + t(Yt2Scale)%*%Yt2Scale/2)
  return(obj)
}


Variance_a0_Ignore = function(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H) {
  nc1 = length(Yc1)
  nt2 = length(Yt2)
  Zt2c = t( t(Zt2) - colMeans(Zt2) )
  Yt2c = Yt2 - mean(Yt2)
  Cov = t(Zc1)%*%Zc1/sigma1H^2 + t(Zt2c)%*%Zt2c/sigma2H^2
  
  c1 = colMeans(Zt2)%*%solve(Cov, t(Zc1)/sigma1H^2)
  c2 = 1/nt2 - colMeans(Zt2)%*%solve(Cov, t(Zt2c)/sigma2H^2)
  a0Var = (sigma1H^2)*c1%*%t(c1) + (sigma2H^2)*c2%*%t(c2)
  return(a0Var)
}


Estimate_Ignore_S1 = function(Zc1, Zt2, Yc1, Yt2, tol.c = 1e-7,
                              a0.Ini = NULL, beta.Ini = NULL, sigma1.Ini = NULL, sigma2.Ini = NULL) {
  nc1 = length(Yc1)
  nt2 = length(Yt2)
  if ( is.null(sigma1.Ini) ) {
         sigma1.Ini = sqrt( mean( (Yc1 - mean(Yc1))^2 ) )
  } 
  
  if ( is.null(sigma2.Ini) ) {
         sigma2.Ini = sqrt( mean( (Yt2 - mean(Yt2))^2 ) )
  }
  
  if ( is.null(beta.Ini) ) {
         beta.Ini = solve( t(Zc1)%*%Zc1, t(Zc1)%*%Yc1)
  }
  
  if ( is.null(a0.Ini) ) {
         a0.Ini = mean(Yt2 - Zt2%*%beta.Ini)
  }
  
  a0H = a0.Ini
  betaH = beta.Ini
  sigma1H = sigma1.Ini
  sigma2H = sigma2.Ini
  
  obj_old = Ignore_Objective_S1(Zc1, Zt2, Yc1, Yt2, a0H, betaH, sigma1H, sigma2H)
  objVec = obj_old 
  i = 0
  gap = 1e7
  
  start = proc.time()[1]
  while ( (i < 100)&&(gap > tol.c) ) {
        i = i + 1
        sigma1H = Ignore_Update_sigma1(Zc1, Yc1, betaH)
        sigma2H = Ignore_Update_sigma2(Zt2, Yt2, betaH)
        betaH =Ignore_Update_beta(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H)
        a0H = Ignore_Update_a0(Zt2, Yt2, betaH)
        obj_new = Ignore_Objective_S1(Zc1, Zt2, Yc1, Yt2, a0H, betaH, sigma1H, sigma2H)
        gap = abs(obj_old - obj_new)
        objVec = c(objVec, obj_new)
        obj_old = obj_new
  }
  Time = proc.time()[1] - start
  a0Var = Variance_a0_Ignore(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H)
  return(list("a0" = a0H, "a0Var" = a0Var, "beta" = betaH, "rho" = NULL, "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec,
                            "Time" = Time))
}

oneReplicate_Ignore = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./Simulation/oneReplicate-New-S1.R")
  Estimate = batch.Ignore.S1(Y, X, Z)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = NULL
  betaH = Estimate$beta
  rhoH = NULL
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "rho" = rhoH, "beta"= betaH, "objVec" = objVec, "Time" = Time))
}

oneReplicateWrap_Ignore = function(seedJ) {
  eval = oneReplicate_Ignore(seedJ) 
  return(eval)
}


batch.Ignore.S1 = function(Y, X, Z) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  Estimate = Estimate_Ignore_S1(Zc1, Zt2, Yc1, Yt2)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  betaH = Estimate$beta
  rhoH = NULL
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
  Time = Estimate$Time
  return(list("a0" = a0H, "a0Var" = a0Var, "beta" = betaH, "rho" = rhoH, "p.value" = pv,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec, "Time" = Time))
}





