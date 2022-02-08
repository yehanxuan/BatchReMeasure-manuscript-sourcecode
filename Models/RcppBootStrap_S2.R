#### RcppBootStrap_S2 ####
RcppReMeasure_Estimate_S2_Pair = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, 
                                          Index_C, Index_T, tol.c = 1e-7, B,
                                          a0.Ini = NULL, a1.Ini = NULL, a3.Ini = NULL, 
                                          rho1.Ini = NULL, rho2.Ini = NULL, beta.Ini = NULL,
                                          sigma1.Ini = NULL, sigma2.Ini = NULL, sigma3.Ini = NULL) {
  start = proc.time()[1]
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  out = RcppReMeasure_Estimate_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, 
                                  Index_C, Index_T, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  a3_hat = out$a3
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  sigma3_hat = out$sigma3
  rho1_hat = out$rho1
  rho2_hat = out$rho2
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec
  
  ztest_b = rep(NA, B)
  for (j in 1:B) {
    ind_nc1_un = sample(1:(nc1-nc3), nc1-nc3, replace=TRUE)
    ind_nc3 = sample(1:nc3, nc3, replace=TRUE)
    Zc1_b = Zc1
    Zc1_b[Index_C, ] = Zc1[Index_C, , drop=F][ind_nc3, , drop = F]
    Zc1_b[-Index_C, ] = Zc1_b[-Index_C, , drop=F][ind_nc1_un, , drop = F]
    Yc1_b = as.matrix(Yc1)
    Yc1_b[Index_C] = Yc1[Index_C][ind_nc3]
    Yc1_b[-Index_C] = Yc1[-Index_C][ind_nc1_un]
    
    Yc3_b = Yc3[ind_nc3]
    Zc3_b = Zc3[ind_nc3, ,drop=F]
    
    ind_nt2_un = sample(1:(nt2 - nt3), nt2-nt3, replace=TRUE)
    ind_nt3 = sample(1:nt3, nt3, replace=TRUE)
    Zt2_b = Zt2
    Zt2_b[Index_T, ] = Zt2[Index_T, , drop=F][ind_nt3, , drop=F]
    Zt2_b[-Index_T, ] = Zt2[-Index_T, , drop = F][ind_nt2_un, , drop=F]
    Yt2_b = Yt2
    Yt2_b[Index_T] = Yt2[Index_T][ind_nt3]
    Yt2_b[-Index_T] = Yt2[-Index_T][ind_nt2_un]
    
    Yt3_b = Yt3[ind_nt3]
    Zt3_b = Zt3[ind_nt3, ,drop=F]
    out_b = RcppReMeasure_Estimate_S2(Zc1_b, Zt2_b, Zc3_b, Zt3_b, 
                                      Yc1_b, Yt2_b, Yc3_b, Yt3_b, Index_C, Index_T, tol.c)
    
    ztest_b[j] = (out_b$a0 - a0_hat)/sqrt(out_b$a0Var)
  }
  Time = proc.time()[1] - start
  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "a3"=a3_hat,
              "beta" = beta_hat, "rho1" = rho1_hat, "rho2" = rho2_hat,
              "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "sigma3"=sigma3_hat, "Time"=Time,
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}

RcppReMeasure_Estimate_S2_Res = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, 
                                     Index_C, Index_T, tol.c = 1e-7, B,
                                     a0.Ini = NULL, a1.Ini = NULL, a3.Ini = NULL, 
                                     rho1.Ini = NULL, rho2.Ini = NULL, beta.Ini = NULL,
                                     sigma1.Ini = NULL, sigma2.Ini = NULL, sigma3.Ini = NULL) {
  
  start = proc.time()[1]
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  out = RcppReMeasure_Estimate_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, Index_C, Index_T, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  a3_hat = out$a3
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  sigma3_hat = out$sigma3
  rho1_hat = out$rho1
  rho2_hat = out$rho2
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec
  
  ec1_hat = Yc1 - Zc1%*%beta_hat
  et2_hat = Yt2 - a0_hat - a1_hat - Zt2%*%beta_hat
  ec3_hat = Yc3 - a3_hat - Zc3%*%beta_hat 
  et3_hat = Yt3 - a0_hat - a3_hat - Zt3%*%beta_hat
  ec1_hat_re = ec1_hat[Index_C, ,drop=F]
  ec1_hat_unre = ec1_hat[-Index_C, , drop = F]
  et2_hat_re = et2_hat[Index_T, ,drop=F]
  et2_hat_unre = et2_hat[-Index_T, , drop = F]
  
  ztest_b = rep(NA, B)
  for (j in 1:B) {
    ind_nc1_un = sample(1:(nc1-nc3), nc1-nc3, replace=TRUE)
    ind_nc3 = sample(1:nc3, nc3, replace=TRUE)
    Yc1_b = Yc1
    Yc1_b[Index_C] = Zc1[Index_C, ,drop=F]%*%beta_hat + ec1_hat_re[ind_nc3]
    Yc1_b[-Index_C] = Zc1[-Index_C, ,drop=F]%*%beta_hat + ec1_hat_unre[ind_nc1_un]
    Yc3_b = a3_hat + Zc3%*%beta_hat + ec3_hat[ind_nc3]
    
    ind_nt2_un = sample(1:(nt2-nt3), nt2-nt3, replace=TRUE)
    ind_nt3 = sample(1:nt3, nt3, replace=TRUE)
    Yt2_b = Yt2 
    Yt2_b[Index_T] = a0_hat + a1_hat+ Zt2[Index_T, ,drop=F]%*%beta_hat + et2_hat_re[ind_nt3]
    Yt2_b[-Index_T] = a0_hat + a1_hat + Zt2[-Index_T, ,drop=F]%*%beta_hat + et2_hat_unre[ind_nt2_un]
    Yt3_b = a0_hat + a3_hat + Zt3%*%beta_hat + et3_hat[ind_nt3]
    out_b = RcppReMeasure_Estimate_S2(Zc1, Zt2, Zc3, Zt3, Yc1_b, Yt2_b, Yc3_b, Yt3_b,
                                      Index_C, Index_T, tol.c)
    ztest_b[j] = (out_b$a0-a0_hat)/sqrt(out_b$a0Var)
  }
  Time = proc.time()[1] - start
  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "a3"=a3_hat,
              "beta" = beta_hat, "rho1" = rho1_hat, "rho2" = rho2_hat,
              "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "sigma3"=sigma3_hat, "Time"=Time,
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


RcppReMeasure_Estimate_S2_Wild = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, 
                                      Index_C, Index_T, tol.c = 1e-7, B,
                                      a0.Ini = NULL, a1.Ini = NULL, a3.Ini = NULL, 
                                      rho1.Ini = NULL, rho2.Ini = NULL, beta.Ini = NULL,
                                      sigma1.Ini = NULL, sigma2.Ini = NULL, sigma3.Ini = NULL) {
  start = proc.time()[1]
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  out = RcppReMeasure_Estimate_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, Index_C, Index_T, tol.c)
  a0_hat = out$a0
  a0Var = out$a0Var
  a1_hat = out$a1
  a3_hat = out$a3
  beta_hat = out$beta
  sigma1_hat = out$sigma1
  sigma2_hat = out$sigma2
  sigma3_hat = out$sigma3
  rho1_hat = out$rho1
  rho2_hat = out$rho2
  ztest = a0_hat/sqrt(out$a0Var)
  objVec = out$objVec
  
  ztest_b = rep(NA, B)
  for (j in 1:B) {
    ec1 = rnorm(nc1, mean = 0, sd = sigma1_hat)
    Yc1_b = ec1 + Zc1%*%beta_hat 
    et2 = rnorm(nt2, mean = 0, sd = sigma2_hat)
    Yt2_b = a0_hat + a1_hat + et2 + Zt2%*%beta_hat
    ecInd = rnorm(nc3, mean = 0, sd = sigma1_hat)
    ec3 = (rho1_hat * ec1[Index_C] + sqrt(1 - rho1_hat^2) * ecInd) * sigma3_hat/sigma1_hat
    Yc3_b = a3_hat + ec3 + Zc3%*%beta_hat
    
    etInd = rnorm(nt3, mean = 0, sd = sigma2_hat) 
    et3 = (rho2_hat * et2[Index_T] + sqrt(1 - rho2_hat^2) * etInd ) * sigma3_hat/sigma2_hat
    Yt3_b = a0_hat + a3_hat + et3 + Zt3%*%beta_hat 
    out_b = RcppReMeasure_Estimate_S2(Zc1, Zt2, Zc3, Zt3, Yc1_b, Yt2_b, Yc3_b, Yt3_b,
                                      Index_C, Index_T, tol.c)
    ztest_b[j] = (out_b$a0-a0_hat)/sqrt(out_b$a0Var)
    # if ( j%%200 == 0) {
    #  print(j)
    # }
  }
  
  Time = proc.time()[1] - start
  return(list("a0" = a0_hat, "a0Var" = a0Var, "a1" = a1_hat, "a3"=a3_hat,
              "beta" = beta_hat, "rho1" = rho1_hat, "rho2" = rho2_hat,
              "sigma1" = sigma1_hat, "sigma2" = sigma2_hat, "sigma3"=sigma3_hat, "Time"=Time,
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b))
}


Rcpp.batch.ReMeasure.S2.Pair = function(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
  ind0 <- X == 0  ## Control Group 
  ind1 <- X == 1
  Yc1 = Y[ind0]
  Yt2 = Y[ind1]
  Zc1 = Z[ind0, , drop=F]
  Zt2 = Z[ind1, , drop=F]
  Zc3 = Zc1[ind.r1, , drop = F]
  Yc3 = Y.r1
  Zt3 = Zt2[ind.r2, , drop = F]
  Yt3 = Y.r2
  
  Estimate = RcppReMeasure_Estimate_S2_Pair(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, 
                                            ind.r1, ind.r2, tol.c = 1e-7, B = 1000)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  a3H = Estimate$a3
  betaH = Estimate$beta
  rho1H = Estimate$rho1
  rho2H = Estimate$rho2 
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  sigma3H = Estimate$sigma3
  objVec = Estimate$objVec
  pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
  Time = Estimate$Time
  ztest = Estimate$ztest 
  ztest_b = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1"=a1H, "a3"=a3H, "beta" = betaH,
              "rho1" = rho1H, "rho2" = rho2H, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H, "Time" = Time, 
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b, "p.value" = pv))
}


Rcpp.batch.ReMeasure.S2.res = function(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
  ind0 <- X == 0  ## Control Group 
  ind1 <- X == 1
  Yc1 = Y[ind0]
  Yt2 = Y[ind1]
  Zc1 = Z[ind0, , drop=F]
  Zt2 = Z[ind1, , drop=F]
  
  Zc3 = Zc1[ind.r1, , drop = F]
  Yc3 = Y.r1
  Zt3 = Zt2[ind.r2, , drop = F]
  Yt3 = Y.r2
  Estimate = RcppReMeasure_Estimate_S2_Res(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, 
                                           ind.r1, ind.r2, tol.c = 1e-7, B = 1000)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  a3H = Estimate$a3
  betaH = Estimate$beta
  rho1H = Estimate$rho1
  rho2H = Estimate$rho2 
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  sigma3H = Estimate$sigma3
  objVec = Estimate$objVec
  pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
  Time = Estimate$Time
  ztest = Estimate$ztest 
  ztest_b = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1"=a1H, "a3"=a3H, "beta" = betaH,
              "rho1" = rho1H, "rho2" = rho2H, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H, "Time" = Time, 
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b, "p.value" = pv))
}


Rcpp.batch.ReMeasure.S2.wild = function(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
  ind0 <- X == 0  ## Control Group 
  ind1 <- X == 1
  Yc1 = Y[ind0]
  Yt2 = Y[ind1]
  Zc1 = Z[ind0, , drop=F]
  Zt2 = Z[ind1, , drop=F]
  
  Zc3 = Zc1[ind.r1, , drop = F]
  Yc3 = Y.r1
  Zt3 = Zt2[ind.r2, , drop = F]
  Yt3 = Y.r2
  Estimate = RcppReMeasure_Estimate_S2_Wild(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3,
                                            ind.r1, ind.r2, tol.c = 1e-7, B = 1000)
  
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  a3H = Estimate$a3
  betaH = Estimate$beta
  rho1H = Estimate$rho1
  rho2H = Estimate$rho2 
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  sigma3H = Estimate$sigma3
  objVec = Estimate$objVec
  pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
  Time = Estimate$Time
  ztest = Estimate$ztest 
  ztest_b = Estimate$ztest_b
  return(list("a0" = a0H, "a0Var" = a0Var, "a1"=a1H, "a3"=a3H, "beta" = betaH,
              "rho1" = rho1H, "rho2" = rho2H, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H, "Time" = Time, 
              "objVec" = objVec, "ztest" = ztest, "ztest_b" = ztest_b, "p.value" = pv))
}


oneReplicate_RcppPair_S2 = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S2.R")
  Estimate = Rcpp.batch.ReMeasure.S2.Pair(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var ## variance (Not from bootstrap)
  a1H = Estimate$a1
  a3H = Estimate$a3
  betaH = Estimate$beta
  rho1H = Estimate$rho1 
  rho2H = Estimate$rho2 
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  sigma3H = Estimate$sigma3
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztestb = Estimate$ztest_b
  pv = Estimate$p.value
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "sigma3" = sigma3H, "rho1" = rho1H, "rho2" = rho2H,
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "ztest" = ztest, "ztestb" = ztestb,
              "p.value" = pv))
}


oneReplicateWrap_RcppPair_S2 = function(seedJ) {
  eval = oneReplicate_RcppPair_S2(seedJ)
  return(eval)
}

oneReplicate_RcppRes_S2 = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S2.R")
  Estimate = Rcpp.batch.ReMeasure.S2.res(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var ## variance (Not from bootstrap)
  a1H = Estimate$a1
  a3H = Estimate$a3
  betaH = Estimate$beta
  rho1H = Estimate$rho1 
  rho2H = Estimate$rho2 
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  sigma3H = Estimate$sigma3
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztestb = Estimate$ztest_b
  pv = Estimate$p.value
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "sigma3" = sigma3H, "rho1" = rho1H, "rho2" = rho2H,
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "ztest" = ztest, "ztestb" = ztestb,
              "p.value" = pv))
}




oneReplicateWrap_RcppRes_S2 = function(seedJ) {
  eval = oneReplicate_RcppRes_S2(seedJ)
  return(eval)
}


oneReplicate_RcppWild_S2 = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S2.R")
  Estimate = Rcpp.batch.ReMeasure.S2.wild(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var ## variance (Not from bootstrap)
  a1H = Estimate$a1
  a3H = Estimate$a3
  betaH = Estimate$beta
  rho1H = Estimate$rho1 
  rho2H = Estimate$rho2 
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  sigma3H = Estimate$sigma3
  objVec = Estimate$objVec
  Time = Estimate$Time
  ztest = Estimate$ztest
  ztestb = Estimate$ztest_b
  pv = Estimate$p.value
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "sigma3" = sigma3H, "rho1" = rho1H, "rho2" = rho2H,
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "ztest" = ztest, "ztestb" = ztestb,
              "p.value" = pv))
}


oneReplicateWrap_RcppWild_S2 = function(seedJ) {
  eval = oneReplicate_RcppWild_S2(seedJ)
  return(eval)
}

