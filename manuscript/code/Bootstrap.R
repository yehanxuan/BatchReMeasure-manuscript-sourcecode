
Estimate_ReMeasure_S1_Pair = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7, B, a0.Ini = NULL, a1.Ini=NULL, 
                                      rho.Ini=NULL, beta.Ini = NULL, sigma1.Ini=NULL, sigma2.Ini=NULL) {
  
  start = proc.time()[1]
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  out = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  rho_hat = out$rho
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec
  
  ztest_b = rep(NA, B)
  
  for (j in 1:B){
    ind_nc1_un = sample(1:(nc1-nc2), nc1-nc2, replace=TRUE)
    ind_nt2 = sample(1:nt2, nt2, replace=TRUE)
    ind_nc2 = sample(1:nc2, nc2, replace=TRUE)
    
    
    Zc1_b = Zc1
    Zc1_b[Index, ] = Zc1[Index, ,drop = F][ind_nc2, , drop = F]
    Zc1_b[-Index, ] = Zc1[-Index, , drop = F][ind_nc1_un, , drop = F]
    Yc1_b = as.matrix(Yc1) 
    Yc1_b[Index] = Yc1[Index][ind_nc2]
    Yc1_b[-Index] = Yc1[-Index][ind_nc1_un]
    
    Yt2_b = Yt2[ind_nt2]
    Zt2_b = Zt2[ind_nt2, , drop = F]
    Yc2_b = Yc2[ind_nc2]
    Zc2_b = Zc2[ind_nc2, , drop = F]
    out_b = Estimate_ReMeasure_S1(Zc1_b, Zt2_b, Zc2_b, Yc1_b, Yt2_b, Yc2_b, Index, tol.c)
    ztest_b[j] = (out_b$a0-a0_hat)/sqrt(out_b$a0Var)
    if ( j%%200 == 0) {
      print(j)
    }
  }
  Time = proc.time()[1] - start
  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "beta" = beta_hat,
              "rho" = rho_hat, "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "Time" = Time, 
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


Estimate_ReMeasure_S1_Res = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7, B, a0.Ini = NULL, a1.Ini=NULL, 
                                      rho.Ini=NULL, beta.Ini = NULL, sigma1.Ini=NULL, sigma2.Ini=NULL) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  # B Time of bootstrap 
  start = proc.time()[1]
  out = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  rho_hat = out$rho
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec
  ## Estimation of noise 
  ec1_hat = Yc1 - Zc1%*%beta_hat
  et2_hat = Yt2 - a0_hat - a1_hat - Zt2%*%beta_hat
  ec2_hat = Yc2 - a1_hat - Zc2%*%beta_hat
  ec1_hat_re = ec1_hat[Index, ,drop=F]
  ec1_hat_unre = ec1_hat[-Index, ,drop=F]
  
  ztest_b = rep(NA,B)
  ## Using bootstrap to estimate the standard error
  for (j in 1:B){
    ind_nc1_un = sample(1:(nc1-nc2), nc1-nc2, replace=TRUE)
    ind_nt2 = sample(1:nt2, nt2, replace=TRUE)
    ind_nc2 = sample(1:nc2, nc2, replace=TRUE)
    Yc1_b = as.matrix(Yc1) 
    Yc1_b[-Index,] = Zc1[-Index,]%*%beta_hat + ec1_hat_unre[ind_nc1_un, ,drop=F]
    Yc1_b[Index,] = Zc1[Index,]%*%beta_hat + ec1_hat_re[ind_nc2, ,drop=F] 
    Yt2_b = a0_hat + a1_hat + Zt2%*%beta_hat + et2_hat[ind_nt2, ,drop=F]
    Yc2_b = a1_hat + Zc2%*%beta_hat + ec2_hat[ind_nc2, ,drop=F] 
    out_b = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1_b, Yt2_b, Yc2_b, Index, tol.c)
    ztest_b[j] = (out_b$a0-a0_hat)/sqrt(out_b$a0Var)
    if ( j%%200 == 0) {
      print(j)
    }
  }
  Time = proc.time()[1] - start

  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "beta" = beta_hat,
              "rho" = rho_hat, "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "Time" = Time, 
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


Estimate_ReMeasure_S1_Wild = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7, B, a0.Ini = NULL, a1.Ini=NULL, 
                                         rho.Ini=NULL, beta.Ini = NULL, sigma1.Ini=NULL, sigma2.Ini=NULL) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  start = proc.time()[1]
  out = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  rho_hat = out$rho
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec
  
  ztest_b = rep(NA,B)
  for (j in 1:B) {
    ec1 = rnorm(nc1, mean = 0, sd = sigma1_hat)
    Yc1_b = ec1 + Zc1%*%beta_hat
    et2 = rnorm(nt2, mean = 0, sd = sigma2_hat)
    Yt2_b = a0_hat + a1_hat + et2 + Zt2%*%beta_hat
    ecInd = rnorm(nc2, mean = 0, sd = sigma1_hat)
    ec2 = (rho_hat * ec1[Index] + sqrt(1 - rho_hat^2) * ecInd) * sigma2_hat/sigma1_hat
    Yc2_b = a1_hat + ec2 + Zc2%*%beta_hat
    out_b = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1_b, Yt2_b, Yc2_b, Index, tol.c)
    ztest_b[j] = (out_b$a0-a0_hat)/sqrt(out_b$a0Var)
    if ( j%%200 == 0) {
      print(j)
    }
  }
  
  Time = proc.time()[1] - start
  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "beta" = beta_hat,
              "rho" = rho_hat, "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "Time" = Time, 
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


batch.ReMeasure.S1.pair = function(Y, X, Z, ind.r, Y.r) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = Estimate_ReMeasure_S1_Pair(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, B = 1000)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztest_b = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time, "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}

batch.ReMeasure.S1.res = function(Y, X, Z, ind.r, Y.r) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = Estimate_ReMeasure_S1_Res(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, B = 1000)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztest_b = Estimate$ztest_b
  
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time, "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


batch.ReMeasure.S1.wild = function(Y, X, Z, ind.r, Y.r) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  Estimate = Estimate_ReMeasure_S1_Wild(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r, B = 1000)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztestb = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time, "objVec" = objVec, "ztest" = ztest, "ztestb" = ztestb))
}



oneReplicate_Pair = function(seedJ){
  set.seed(seedJ + repID * 300)
  source("./Simulation/oneReplicate-New-S1.R")
  Estimate = batch.ReMeasure.S1.pair(Y, X, Z, ind.r, Y.r)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var ## variance (Not from bootstrap)
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztestb = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "rho" = rhoH, 
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "ztest" = ztest, "ztestb" = ztestb))
}

oneReplicateWrap_Pair = function(seedJ) {
  eval = oneReplicate_Pair(seedJ)
  return(eval)
}


oneReplicate_Res = function(seedJ){
  set.seed(seedJ + repID * 300)
  source("./Simulation/oneReplicate-New-S1.R")
  Estimate = batch.ReMeasure.S1.res(Y, X, Z, ind.r, Y.r)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var ## variance (Not from bootstrap)
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztestb = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "rho" = rhoH, 
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "ztest" = ztest, "ztestb" = ztestb))
}

oneReplicateWrap_Res = function(seedJ){
  eval = oneReplicate_Res(seedJ)
  return(eval)
}

oneReplicate_Wild = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./Simulation/oneReplicate-New-S1.R")
  Estimate = batch.ReMeasure.S1.wild(Y, X, Z, ind.r, Y.r)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var ## variance (Not from bootstrap)
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztestb = Estimate$ztestb
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "rho" = rhoH, 
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "ztest" = ztest, "ztestb" = ztestb))
}

oneReplicateWrap_Wild = function(seedJ){
  eval = oneReplicate_Wild(seedJ)
  return(eval)
}
  
  
  