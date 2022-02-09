RcppReMeasure_Estimate_S2 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, Index_C, Index_T, 
                                     tol.c = 1e-7, a0.Ini = NULL, a1.Ini = NULL, a3.Ini = NULL,
                                     rho1.Ini=NULL, rho2.Ini=NULL, beta.Ini = NULL, 
                                     sigma1.Ini = NULL, sigma2.Ini=NULL, sigma3.Ini = NULL) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  if ( is.null(sigma1.Ini) ) {
    sigma1.Ini = sqrt( mean( (Yc1[Index_C] - mean(Yc1[Index_C]) )^2 ) )
  } 
  if (sigma1.Ini <= 0.01 ){
    sigma1.Ini = 0.2
  }
  if ( is.null(sigma2.Ini) ) {
    sigma2.Ini = sqrt( mean( (Yt2[Index_T] - mean(Yt2[Index_T]))^2 ) )
  }
  if (sigma2.Ini <= 0.01) {
    sigma2.Ini = 0.2
  }
  if ( is.null(sigma3.Ini)) {
    sigma3.Ini = sqrt( mean( (Yc3 - mean(Yc3))^2 ) )
  }
  if (sigma3.Ini <= 0.01 ){
    sigma3.Ini = 0.2
  }
  if ( is.null(rho1.Ini) ) {
    rho1.Ini = t(Yc1[Index_C] - mean(Yc1[Index_C]))%*%(Yc3 - mean(Yc3))/(nc3*sigma1.Ini*sigma3.Ini)
  }
  if ( rho1.Ini > 0.99 ) {
    rho1.Ini = 0.95
  } else if (rho1.Ini <= -0.99) {
    rho1.Ini = -0.95
  }
  
  ### Set another estimate for sigma3 in case problem happens in PairBoot
  sigma4.Ini = NULL
  if ( is.null(sigma4.Ini)) {
    sigma4.Ini = sqrt( mean( (Yt3 - mean(Yt3))^2 ) )
  }
  if (sigma4.Ini <= 0.01 ){
    sigma4.Ini = 0.2
  }
  
  if ( is.null(rho2.Ini) ) {
    # rho2.Ini = t( Yt2[Index_T] - mean(Yt2[Index_T]) )%*%(Yt3 - mean(Yt3))/(nt3*sigma2.Ini* (sqrt( mean( (Yt3 - mean(Yt3))^2 ) )) )
    rho2.Ini = t( Yt2[Index_T] - mean(Yt2[Index_T]) )%*%(Yt3 - mean(Yt3))/(nt3*sigma2.Ini* sigma4.Ini)
  }
  if ( rho2.Ini > 0.99 ) {
    rho2.Ini = 0.95
  } else if (rho2.Ini <= -0.99) {
    rho2.Ini = -0.95
  }
  
  if (is.null(beta.Ini)) {
    beta.Ini = solve( t(Zc1)%*%Zc1, t(Zc1)%*%Yc1)
  }
  if (is.null(a3.Ini)) {
    a3.Ini = mean(Yc3 - Zc3%*%beta.Ini)
  }
  if (is.null(a0.Ini)) {
    a0.Ini = mean(Yt3 - Zt3%*%beta.Ini) - a3.Ini
  }
  if (is.null(a1.Ini)) {
    a1.Ini = mean(Yt2 - Zt2%*%beta.Ini) - a0.Ini
  }
  a0H = a0.Ini
  a1H = a1.Ini
  a3H = a3.Ini
  rho1H = rho1.Ini
  rho2H = rho2.Ini
  sigma1H = sigma1.Ini
  sigma2H = sigma2.Ini
  sigma3H = sigma3.Ini
  betaH = beta.Ini
  
  obj_old = Rcpp_Objective_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H, a1H,
                              a3H, rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T)
  objVec = obj_old
  i = 0
  gap = 1e7
  start_S2 = proc.time()[1]
  
  while ( (i < 100)&&(gap > tol.c) ) {
    i = i + 1
    sigma1H = Update_sigma1_S2(Zc1, Zc3, Yc1, Yc3, a3H, betaH, rho1H, sigma3H, Index_C)
    sigma2H = Update_sigma2_S2(Zt2, Zt3, Yt2, Yt3, a0H, a1H, a3H,betaH, rho2H, sigma3H, Index_T)
    sigma3H = Update_sigma3_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, a0H, a1H, a3H, 
                               betaH, rho1H, rho2H, sigma1H, sigma2H, Index_C, Index_T)
    rho1H = Update_rho1_S2(Zc1, Zc3, Yc1, Yc3, a3H, betaH, sigma1H, sigma3H, Index_C)
    rho2H = Update_rho2_S2(Zt2, Zt3, Yt2, Yt3, a0H, a1H, a3H, betaH, sigma2H, sigma3H, Index_T)
    betaH = Update_beta_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, rho1H, rho2H,
                           sigma1H, sigma2H, sigma3H, Index_C, Index_T)
    aVec = Update_a_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, rho1H, rho2H, 
                       sigma1H, sigma2H, sigma3H, Index_C, Index_T)
    a0H = aVec$a0H
    a1H = aVec$a1H
    a3H = aVec$a3H
    obj_new = Rcpp_Objective_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH,
                                a0H, a1H, a3H, rho1H, rho2H, sigma1H, sigma2H, sigma3H,
                                Index_C, Index_T)
    gap = abs(obj_old - obj_new)
    objVec = c(objVec, obj_new)
    obj_old = obj_new 
  }
  Time_S2 = proc.time()[1] - start_S2
  a0Var = S2_Variance_a0(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, sigma1H, sigma2H, 
                         sigma3H, rho1H, rho2H, Index_C, Index_T)
  
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H, "beta" = betaH,
              "rho1" = rho1H, "rho2" = rho2H, "sigma1" = sigma1H, "sigma2" = sigma2H,
              "sigma3" = sigma3H, "Time" = Time_S2, "objVec" = objVec) )
}


Rcpp.batch.ReMeasure.S2 = function(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
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
  
  Estimate = RcppReMeasure_Estimate_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, ind.r1, ind.r2,
                                       tol.c = 1e-7)
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
  return(list("a0" = a0H, "a0Var" = a0Var, "a1"=a1H, "a3"=a3H, "beta" = betaH,
              "rho1" = rho1H, "rho2" = rho2H, "p.value" = pv, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H, "objVec" = objVec, "Time" = Time))
}

oneReplicate_ReMeasure_Rcpp_S2 = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S2.R")
  Estimate =  Rcpp.batch.ReMeasure.S2(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2)
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
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H,
              "rho1" = rho1H, "rho2" = rho2H, 
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "p.value" = pv))
}


oneReplicateWrap_ReMeasure_Rcpp_S2 = function(seedJ) {
  eval = oneReplicate_ReMeasure_Rcpp_S2(seedJ)
  return(eval)
}



