#### IgnoreBatch_S2 ####

S2_Ignore_Update_sigma1 = function(Zc1, Yc1, betaH) {
  nc1 = nrow(Zc1)
  sigma1H = sqrt( mean( (Yc1 - Zc1%*%betaH)^2 ) ) 
  return(sigma1H)
}

S2_Ignore_Update_sigma2 = function(Zt2, Yt2, betaH){
  nt2 = nrow(Zt2)
  Zt2c = t( t(Zt2) - colMeans(Zt2) ) 
  Yt2c = Yt2 - mean(Yt2)
  sigma2H = sqrt(mean( (Yt2c - Zt2c%*%betaH)^2 ))
  return(sigma2H)
}


S2_Ignore_Update_beta = function(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H) {
  nt2 = nrow(Zt2)
  Zt2c = t( t(Zt2) - colMeans(Zt2) ) 
  Yt2c = Yt2 - mean(Yt2)
  Cov = t(Zc1)%*%Zc1/sigma1H^2 + t(Zt2c)%*%Zt2c/sigma2H^2 
  Cor = t(Zc1)%*%Yc1/sigma1H^2 + t(Zt2c)%*% Yt2c/sigma2H^2
  betaH = solve(Cov, Cor)
  return(betaH)            
}



S2_Ignore_Update_a0 = function(Zt2, Yt2, betaH) {
  a0H = mean(Yt2 - Zt2%*%betaH)
  return(a0H)
}


Ignore_Objective_S2 = function(Zc1, Zt2, Yc1, Yt2, a0H, betaH, sigma1H, sigma2H){
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


S2_Variance_a0_Ignore = function(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H) {
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


Estimate_Ignore_S2 = function(Zc1, Zt2, Yc1, Yt2, tol.c = 1e-7,
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
  
  obj_old = Ignore_Objective_S2(Zc1, Zt2, Yc1, Yt2, a0H, betaH, sigma1H, sigma2H)
  objVec = obj_old
  i = 0
  gap = 1e7
  
  start = proc.time()[1]
  while ( (i < 100)&&(gap > tol.c) ) { 
    i = i + 1
    sigma1H = S2_Ignore_Update_sigma1(Zc1, Yc1, betaH)
    sigma2H = S2_Ignore_Update_sigma2(Zt2, Yt2, betaH)
    betaH = S2_Ignore_Update_beta(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H)
    a0H = S2_Ignore_Update_a0(Zt2, Yt2, betaH)
    obj_new = Ignore_Objective_S2(Zc1, Zt2, Yc1, Yt2, a0H, betaH, sigma1H, sigma2H)
    gap = abs(obj_old - obj_new)
    objVec = c(objVec, obj_new)
    obj_old = obj_new
  }
  Time = proc.time()[1] - start
  a0Var = S2_Variance_a0_Ignore(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H)
  return(list("a0" = a0H, "a0Var" = a0Var, "beta" = betaH, "rho" = NULL, "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec,
              "Time" = Time))
}



oneReplicate_Ignore_S2 = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S2.R")
  
  Estimate = batch.Ignore.S2(Y, X, Z)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = NULL
  a3H = NULL
  betaH = Estimate$beta
  rho1H = NULL
  rho2H = NULL
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  sigma3H = Estimate$sigma3
  objVec = Estimate$objVec
  Time = Estimate$Time
  pv <- Estimate$p.value
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H,
              "rho1" = rho1H, "rho2" = rho2H,
              "beta"= betaH, "objVec" = objVec, "Time" = Time, "p.value" = pv))
}



oneReplicateWrap_Ignore_S2 = function(seedJ) {
  eval = oneReplicate_Ignore_S2(seedJ) 
  return(eval)
}


batch.Ignore.S2 = function(Y, X, Z) {
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  Estimate = Estimate_Ignore_S2(Zc1, Zt2, Yc1, Yt2)
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


# S2_Ignore_Update_rho1 = function(Zc1, Zc3, Yc1, Yc3, betaH, sigma1H, sigma3H, Index_C){
#   nc3 = nrow(Zc3)
#   mean1 = Zc1%*%betaH
#   mean3 = Zc3%*%betaH
#   
#   W1p = sum( (Yc1[Index_C] - mean1[Index_C])^2 )
#   W3p = sum( (Yc3 - mean3)^2 )
#   W13p = sum( (Yc1[Index_C] - mean1[Index_C])*(Yc3 - mean3) )
#   # Solve cubic equation
#   c3 = 1
#   c2 = -W13p/(nc3*sigma1H*sigma3H)
#   c1 = W1p/(nc3*sigma1H*sigma1H) + W3p/(nc3*sigma3H*sigma3H) - 1
#   c0 = -W13p/(nc3*sigma1H*sigma3H)
#   
#   rho1H = cubic( c(c3, c2, c1, c0) )[1]
#   rho1H = Re(rho1H) 
#   if (rho1H > 0.999) {
#     rho1H = 0.99
#   } else if (rho1H < -0.999 ) {
#     rho1H = -0.99
#   }
#   return(rho1H)
# }
# 
# S2_Ignore_Update_rho2 = function(Zt2, Zt3, Yt2, Yt3, a0H, betaH, 
#                                  sigma2H, sigma3H, Index_T) {
#   nt3 = nrow(Zt3)
#   mean2 = a0H + Zt2 %*% betaH
#   mean4 = a0H + Zt3 %*% betaH
#   
#   W2p = sum( (Yt2[Index_T] - mean2[Index_T])^2 )
#   W4p = sum( (Yt3 - mean4)^2 )
#   W24p = sum( (Yt2[Index_T] - mean2[Index_T])*(Yt3 - mean4) )
#   
#   c3 = 1
#   c2 = -W24p/(nt3*sigma2H*sigma3H)
#   c1 = W2p/(nt3*sigma2H*sigma2H) + W4p/(nt3*sigma3H*sigma3H) - 1
#   c0 = -W24p/(nt3*sigma2H*sigma3H)
#   
#   rho2H = cubic(c(c3, c2, c1, c0))[1]
#   rho2H = Re(rho2H)
#   
#   if (rho2H > 0.999) {
#     rho2H = 0.99
#   } else if (rho2H < -0.99) {
#     rho2H = -0.99
#   }
#   return(rho2H)
# }
# 
# 
# S2_Ignore_Update_sigma1 = function(Zc1, Zc3, Yc1, Yc3, betaH, rho1H, sigma3H, Index_C){
#   nc1 = nrow(Zc1)
#   mean1 = Zc1%*%betaH
#   mean3 = Zc3%*%betaH
#   
#   W1p = sum( (Yc1[Index_C] - mean1[Index_C])^2)
#   if (length(Index_C) == nc1) {
#     W1s = 0
#   } else {
#     W1s = sum( (Yc1[-Index_C] - mean1[-Index_C])^2 )
#   }
#   
#   W13p = sum( (Yc1[Index_C] - mean1[Index_C])*(Yc3 - mean3) )
#   if (rho1H > 0.999) {
#     rho1H = 0.99
#   } else if (rho1H < -0.999) {
#     rho1H = -0.99
#   }
#   c0 = -W1p/nc1 - (1 - rho1H^2)*W1s/nc1
#   c1 = rho1H*W13p/(nc1*sigma3H)
#   c2 = (1 - rho1H^2)
#   Roots = Re(polyroot(c(c0, c1, c2)))
#   sigma1H = Roots[Roots > 0 ]
#   return(sigma1H)
# }
# 
# S2_Ignore_Update_sigma2 = function(Zt2, Zt3, Yt2, Yt3, a0H, betaH, rho2H, sigma3H, Index_T){
#   nt2 = nrow(Zt2)
#   mean2 = a0H  + Zt2%*%betaH
#   mean4 = a0H  + Zt3%*%betaH
#   W2p = sum( (Yt2[Index_T] - mean2[Index_T])^2 )
#   
#   if (length(Index_T) == nt2) {
#     W2s = 0
#   } else {
#     W2s = sum( (Yt2[-Index_T] - mean2[-Index_T])^2 )  
#   }
#   W24p = sum( (Yt2[Index_T] - mean2[Index_T])*(Yt3 - mean4) )
#   if (rho2H > 0.999) {
#     rho2H = 0.99
#   } else if (rho2H < -0.999) {
#     rho2H = -0.99
#   }
#   
#   c0 = -W2p/nt2 - (1 - rho2H^2)*W2s/nt2
#   c1 = rho2H*W24p/(nt2*sigma3H)
#   c2 = (1 - rho2H^2)
#   Roots = Re(polyroot(c(c0, c1, c2)))
#   sigma2H = Roots[Roots > 0]
#   return(sigma2H)
# }
# 
# 
# S2_Ignore_Update_sigma3 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, a0H,
#                                    betaH, rho1H, rho2H, sigma1H, sigma2H, Index_C, Index_T){
#   nc3 = nrow(Zc3)
#   nt3 = nrow(Zt3)
#   N3 = nc3 + nt3
#   
#   mean1 = Zc1%*%betaH
#   mean2 = a0H + Zt2%*%betaH
#   mean3 = Zc3%*%betaH
#   mean4 = a0H + Zt3%*%betaH
#   
#   W3p = sum( (Yc3 - mean3)^2)
#   W4p = sum( (Yt3 - mean4)^2)
#   W13p = sum( (Yc1[Index_C] - mean1[Index_C])*(Yc3 - mean3) )
#   W24p = sum( (Yt2[Index_T] - mean2[Index_T])*(Yt3 - mean4) )
#   
#   c0 = -(1 - rho2H^2)*W3p/N3 - (1 - rho1H^2)*W4p/N3 
#   c1 = (1 - rho2H^2)*rho1H*W13p/sigma1H  + (1 - rho1H^2)*rho2H*W24p/sigma2H
#   c1 = c1/N3
#   # c2 = 1
#   c2 = (1 - rho1H^2)*(1 - rho2H^2)
#   Roots = Re(polyroot(c(c0, c1, c2)))
#   sigma3H = Roots[Roots > 0]
#   return(sigma3H)
# }
# 
# 
# S2_Ignore_Update_a0 = function( Zt2, Zt3, Yt2, Yt3, betaH, rho1H, rho2H,
#                                sigma1H, sigma2H, sigma3H, Index_C, Index_T) {
#   nt3 = nrow(Zt3)
#   nt2 = nrow(Zt2)
#   R2p = mean( Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
#   R4p = mean(Yt3 - Zt3%*%betaH)
#   if (nt3 == nt2) {
#     R2s = 0
#   } else {
#     R2s = mean(Yt2[-Index_T] - Zt2[-Index_T, , drop = F]%*%betaH)  
#   }
#   w1 = nt3*(1/(sigma2H^2) - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2)
#   w2 = nt3*(1/(sigma3H^2) - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2)
#   w3 = (nt2 - nt3)/(sigma2H^2)
#   w = w1 + w2 + w3 
#   
#   w1 = w1/w 
#   w2 = w2/w
#   w3 = w3/w 
#   
#   a0H = w1*R2p + w2*R4p + w3*R2s  
#   return(a0H)
# }
# 
# S2_Ignore_Update_beta = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3,
#                               rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T) {
#   nc1 = nrow(Zc1)
#   nt2 = nrow(Zt2)
#   nc3 = nrow(Zc3)
#   nt3 = nrow(Zt3)
#   
#   w1 = as.numeric( nt3*(1/(sigma2H^2) - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2) )
#   w2 = as.numeric( nt3*(1/(sigma3H^2) - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2) )
#   w3 = as.numeric( (nt2 - nt3)/(sigma2H^2) )
#   w = w1 + w2 + w3 
#   w1 = w1/w 
#   w2 = w2/w
#   w3 = w3/w
#   
#   Sc1 = as.numeric( (1/sigma1H^2 - rho1H/(sigma1H*sigma3H))/(1 - rho1H^2) )
#   Sc2 = as.numeric( (1/sigma3H^2 - rho1H/(sigma1H*sigma3H))/(1 - rho1H^2) )
#   St1 = as.numeric( (1/sigma2H^2 - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2) )
#   St2 = as.numeric( (1/sigma3H^2 - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2) )
#   
#   wZ = (w1 + w2)*colMeans(Zt3)
#   if (nt3 == nt2) {
#     wZ = wZ
#     Cov4 = 0 
#   } else {
#     Zt2cs = Zt2[-Index_T, , drop = F]
#     wZ = wZ + w3*colMeans(Zt2cs)
#     Cov4 = (t(Zt2cs) - wZ)%*%t( (t(Zt2cs) - wZ) )/sigma2H^2
#   }
#   
#   Cov1 = (Sc1 + Sc2)*t(Zc3)%*%(Zc3)
#   if (nc3 == nc1) {
#     Cov2 = 0
#     Cor2 = 0
#   } else {
#     Zc1cs = Zc1[-Index_C, , drop = F]
#     Cov2 = t(Zc1cs)%*%Zc1cs/(sigma1H^2)
#     Cor2 = t(Zc1cs)%*%Yc1[-Index_C]/(sigma1H^2)
#   }
#   
#   Cov3 = ( t(Zt3) - wZ )%*%t( (t(Zt3) - wZ ) )*(St1 + St2)
#     
#   Cov = Cov1 + Cov2 + Cov3 + Cov4
#   
#   Cor1 = t(Zc3)%*%Yc1[Index_C]*Sc1 +
#     t(Zc3)%*%Yc3*Sc2
# 
#  
#   Cor3 =(t(Zt3) - wZ)%*%(St1*diag(nt3) - (St1 + St2)*w1*rep(1, nt3)%*%t(rep(1, nt3))/nt3) %*%Yt2[Index_T]
#   
#   Cor4 = ( (t(Zt3) - wZ)%*%(St2*diag(nt3) - (St1 + St2)*w2*rep(1, nt3)%*%t(rep(1, nt3))/nt3 ) )%*%Yt3
#   
#   if (nt3 != nt2) {
#     Cor3 = Cor3 - (w1/(sigma2H^2)*(t(Zt2cs) - wZ )%*%rep(1, nt2-nt3)%*%t(rep(1, nt3))/nt3)%*%Yt2[Index_T]
#     Cor4 = Cor4 - (w2/(sigma2H^2)*(t(Zt2cs) - wZ )%*%rep(1, nt2-nt3)%*%t(rep(1, nt3))/nt3)%*%Yt3
#   }
#   
#   if (nt3 == nt2){
#     Cor5 = 0
#   } else {
#     Cor5 = ( ( t(Zt2cs) - wZ )%*%(diag(nt2-nt3) - w3*rep(1, nt2-nt3)%*%t(rep(1, nt2-nt3))/(nt2-nt3))/sigma2H^2-
#       w3*(t(Zt3) - wZ)%*%rep(1, nt3)%*%t(rep(1, nt2-nt3))/(nt2-nt3)*(St1 + St2) )%*%Yt2[-Index_T]
#   }
#   
#   Cor = Cor1 + Cor2 + Cor3 + Cor4 + Cor5
#   
#   betaH = solve(Cov, Cor)
#   return(betaH)
# }
# 
# 
# Ignore_Objective_S2 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H,
#                                rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T){
#   nc3 = nrow(Zc3)
#   nt3 = nrow(Zt3)
#   nc1 = nrow(Zc1)
#   nt2 = nrow(Zt2)
#   mu1 = Zc1 %*% betaH
#   mu2 = a0H + Zt2%*%betaH
#   mu3 = Zc3%*%betaH
#   mu4 = a0H + Zt3%*%betaH
#   
#   Yc1Scale = (Yc1 - mu1)/sigma1H
#   Yt2Scale = (Yt2 - mu2)/sigma2H
#   Yc3Scale = (Yc3 - mu3)/sigma3H
#   Yt3Scale = (Yt3 - mu4)/sigma3H
#   
#   part1 = nc3*log(sigma1H*sigma3H) + nc3*log( 1 - rho1H^2)/2 + 
#     ( t(Yc1Scale[Index_C])%*%Yc1Scale[Index_C] - 2*rho1H*t(Yc1Scale[Index_C])%*%Yc3Scale +
#         t(Yc3Scale)%*%Yc3Scale ) / (2*(1-rho1H^2))
#   
#   if (nc3 == nc1) {
#     part2 = 0
#   } else {
#     part2 = (nc1 - nc3)*log(sigma1H) + t(Yc1Scale[-Index_C])%*%Yc1Scale[-Index_C]/2
#   }
#   
#   
#   part3 = nt3*log(sigma2H*sigma3H) + nt3*log( 1- rho2H^2)/2 + 
#     ( t(Yt2Scale[Index_T])%*%Yt2Scale[Index_T] - 2*rho2H*t(Yt2Scale[Index_T])%*%Yt3Scale +
#         t(Yt3Scale)%*%Yt3Scale )/ (2*(1 - rho2H^2))
#   
#   if (nt3 == nt2) {
#     part4 = 0
#   } else {
#     part4 = (nt2 - nt3)*log(sigma2H) + 
#       t(Yt2Scale[-Index_T]) %*% Yt2Scale[-Index_T]
#   }
#   
#   
#   obj = part1 + part2 + part3 + part4
#   
#   return(obj)
# }
# 
# 
# S2_Variance_a0_Ignore = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, sigma1H, sigma2H, sigma3H, 
#                                  rho1H, rho2H, Index_C, Index_T) {
#   nc1 = nrow(Zc1)
#   nt2 = nrow(Zt2)
#   nc3 = nrow(Zc3)
#   nt3 = nrow(Zt3)
#   
#   w1 = as.numeric( nt3*(1/(sigma2H^2) - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2) )
#   w2 = as.numeric( nt3*(1/(sigma3H^2) - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2) )
#   w3 = as.numeric( (nt2 - nt3)/(sigma2H^2) )
#   w = w1 + w2 + w3 
#   w1 = w1/w 
#   w2 = w2/w
#   w3 = w3/w
#   
#   
#   Sc1 = as.numeric( (1/sigma1H^2 - rho1H/(sigma1H*sigma3H))/(1 - rho1H^2) )
#   Sc2 = as.numeric( (1/sigma3H^2 - rho1H/(sigma1H*sigma3H))/(1 - rho1H^2) )
#   St1 = as.numeric( (1/sigma2H^2 - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2) )
#   St2 = as.numeric( (1/sigma3H^2 - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2) )
#   
#   k = (w1 + w2)*colMeans(Zt3)
#   if (nt3 == nt2) {
#     k = k
#     Cov4 = 0 ### Is it correct
#   } else {
#     Zt2cs = Zt2[-Index_T, , drop = F]
#     k = k + w3*colMeans(Zt2cs)
#     Cov4 = (t(Zt2cs) - k)%*%t( (t(Zt2cs) - k) )/sigma2H^2
#   }
#  
#   Cov1 = (Sc1 + Sc2)*t(Zc3)%*%(Zc3)
#   if (nc3 == nc1) {
#     Cov2 = 0
#     Cor2 = 0
#     B0 = Cor2 
#     B = B0 
#     b = 0
#   } else {
#     Zc1cs = Zc1[-Index_C, , drop = F]
#     Cov2 = t(Zc1cs)%*%Zc1cs/(sigma1H^2)
#     #Cor2 = t(Zc1cs)%*%Yc1[-Index_C]/(sigma1H^2)
#     B0 = t(Zc1cs)/(sigma1H^2)
#   }
#   
#   Cov3 = ( t(Zt3) - k )%*%t( (t(Zt3) - k ) )*
#     (St1 + St2)
#   Cov = Cov1 + Cov2 + Cov3 + Cov4
#   
#   A0 = Sc1*t(Zc3)
#   
#   if (nc3 != nc1) {
#     B = solve(Cov, B0)
#     b = -t(B)%*%k  
#   }
#   
#   
#   C0 =  ( t(Zt3) - k )%*%(St1*diag(nt3) - (St1 + St2)*w1*rep(1, nt3)%*%t(rep(1, nt3))/nt3 ) 
#   E0 = Sc2*t(Zc3)
#   G0 = ( t(Zt3) - k)%*%(St2*diag(nt3) - (St1 + St2)*w2*rep(1, nt3)%*%t(rep(1, nt3))/nt3 )
#   if (nt3 == nt2) {
#     D0 = 0  
#     D = 0
#     d = 0
#   } else {
#     
#     C0 = C0 - w1/(sigma2H^2)*(t(Zt2cs)-k)%*%rep(1, nt2-nt3)%*%t(rep(1, nt3))/nt3
#     D0 =  ( t(Zt2cs) - k )%*%(diag(nt2-nt3) - w3*rep(1, nt2-nt3)%*%t(rep(1, nt2-nt3))/(nt2-nt3))/sigma2H^2
#     D0 = D0 - w3*(t(Zt3) - k)%*%rep(1, nt3)%*%t(rep(1, nt2-nt3))/(nt2-nt3)*(St1 + St2)
#     
#     G0 = G0 - w2/(sigma2H^2)*(t(Zt2cs)-k)%*%rep(1, nt2-nt3)%*%t(rep(1, nt3))/nt3
#     D = solve(Cov, D0)
#     d = (w3/(nt2-nt3) - t(D)%*%k)
#   }
#   
#   
#   
#   A = solve(Cov, A0)
#   C = solve(Cov, C0)
#   
#   E = solve(Cov, E0)
#   G = solve(Cov, G0)
#   
#   a = -t(A)%*%k
#   c = (w1/nt3 - t(C)%*%k)
#   e = -t(E)%*%k
#   f = w2/nt3 - t(G)%*%k
#   
#   Var = (sigma1H^2) * (t(a)%*%a + t(b)%*%b ) + (sigma2H^2)*(t(c)%*%c + t(d)%*%d) +
#     (sigma3H^2)*(t(e)%*%e + t(f)%*%f ) + 2*rho1H*sigma1H*sigma3H*(t(a)%*%e) + 
#     2*rho2H*sigma2H*sigma3H*(t(c)%*%f)
#   return(Var)
# }
# 
# 
# Estimate_Ignore_S2 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, Index_C, Index_T, 
#                               tol.c = 1e-7, a0.Ini = NULL, rho1.Ini=NULL, rho2.Ini=NULL, beta.Ini = NULL, 
#                               sigma1.Ini = NULL, sigma2.Ini = NULL, sigma3.Ini = NULL) {
#   nc1 = nrow(Zc1)
#   nt2 = nrow(Zt2)
#   nc3 = nrow(Zc3)
#   nt3 = nrow(Zt3)
#   if ( is.null(sigma1.Ini) ) {
#     sigma1.Ini = sqrt( mean( (Yc1[Index_C] - mean(Yc1[Index_C]) )^2 ) )
#   } 
#   if (sigma1.Ini <= 0.001 ){
#     sigma1.Ini = 0.5
#   }
#   if ( is.null(sigma2.Ini) ) {
#     sigma2.Ini = sqrt( mean( (Yt2[Index_T] - mean(Yt2[Index_T]))^2 ) )
#   }
#   if (sigma2.Ini <= 0.001) {
#     sigma2.Ini = 0.5
#   }
#   if ( is.null(sigma3.Ini)) {
#     sigma3.Ini = sqrt( mean( (Yc3 - mean(Yc3))^2 ) )
#   }
#   if (sigma3.Ini <= 0.001 ){
#     sigma3.Ini = 0.5
#   }
#   if ( is.null(rho1.Ini) ) {
#     rho1.Ini = t(Yc1[Index_C] - mean(Yc1[Index_C]))%*%(Yc3 - mean(Yc3))/(nc3*sigma1.Ini*sigma3.Ini)
#   }
#   if ( rho1.Ini > 0.99 ) {
#     rho1.Ini = 0.95
#   } else if (rho1.Ini <= -0.99) {
#     rho1.Ini = -0.95
#   }
#   
#   if ( is.null(rho2.Ini) ) {
#     rho2.Ini = t( Yt2[Index_T] - mean(Yt2[Index_T]) )%*%(Yt3 - mean(Yt3))/(nt3*sigma2.Ini* (sqrt( mean( (Yt3 - mean(Yt3))^2 ) )) )
#   }
#   if ( rho2.Ini > 0.99 ) {
#     rho2.Ini = 0.95
#   } else if (rho2.Ini <= -0.99) {
#     rho2.Ini = -0.95
#   }
#   
#   if (is.null(beta.Ini)) {
#     beta.Ini = solve( t(Zc1)%*%Zc1, t(Zc1)%*%Yc1)
#   }
#   
#   rho1H = rho1.Ini
#   rho2H = rho2.Ini
#   sigma1H = sigma1.Ini
#   sigma2H = sigma2.Ini
#   sigma3H = sigma3.Ini
#   betaH = beta.Ini
#   
#   w1 = nt3*(1/(sigma2H^2) - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2)
#   w2 = nt3*(1/(sigma3H^2) - rho2H/(sigma2H*sigma3H))/(1 - rho2H^2)
#   w3 = (nt2 - nt3)/(sigma2H^2)
#   w = w1 + w2 + w3 
#   w1 = w1/w 
#   w2 = w2/w
#   w3 = w3/w
#   
#   R2p = mean( Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
#   R4p = mean(Yt3 - Zt3%*%betaH)
#   if (nt3 == nt2) {
#     R2s = 0
#   } else {
#     R2s = mean(Yt2[-Index_T] - Zt2[-Index_T, , drop = F]%*%betaH)  
#   }
#   
#   if (is.null(a0.Ini)) {
#     a0.Ini = w1*R2p + w2*R4p + w3*R2s
#   }
#   
#   a0H = as.numeric(a0.Ini)
#   
#   
#   obj_old = Ignore_Objective_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH,
#                                 a0H, rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T)
#   objVec = obj_old
#   i = 0
#   gap = 1e7 
#   start_S2 = proc.time()[1]
#   while ( (i < 100)&&(gap > tol.c)  ) {
#     i = i + 1 
#     sigma1H = S2_Ignore_Update_sigma1(Zc1, Zc3, Yc1, Yc3, betaH, rho1H, sigma3H, Index_C)
#     sigma2H = S2_Ignore_Update_sigma2(Zt2, Zt3, Yt2, Yt3, a0H, betaH, rho2H, sigma3H, Index_T)
#     sigma3H = S2_Ignore_Update_sigma3(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, a0H, betaH, rho1H,
#                                       rho2H,sigma1H, sigma2H, Index_C, Index_T)
#     
#     rho1H = S2_Ignore_Update_rho1(Zc1, Zc3, Yc1, Yc3, betaH, sigma1H, sigma3H, Index_C)
#     rho2H = S2_Ignore_Update_rho2(Zt2, Zt3, Yt2, Yt3, a0H, betaH, sigma2H, sigma3H,
#                                   Index_T)
#     betaH = S2_Ignore_Update_beta(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, rho1H, rho2H,
#                                sigma1H, sigma2H, sigma3H, Index_C, Index_T)
#     a0H = S2_Ignore_Update_a0(Zt2, Zt3, Yt2, Yt3, betaH, rho1H, rho2H, sigma1H,
#                               sigma2H, sigma3H, Index_C, Index_T)
#     obj_new = Ignore_Objective_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H,
#                                   rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T)
#     gap = abs(obj_old - obj_new)
#     objVec = c(objVec, obj_new)
#     obj_old = obj_new
#   }
#   Time_S2 = proc.time()[1] - start_S2
#   
#   a0Var = S2_Variance_a0_Ignore(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, sigma1H,
#                                 sigma2H, sigma3H, rho1H, rho2H,
#                                 Index_C, Index_T)
#   return(list("a0" = a0H, "a0Var" = a0Var, "beta" = betaH,
#               "rho1" = rho1H, "rho2" = rho2H, "sigma1" = sigma1H, "sigma2" = sigma2H,
#               "sigma3" = sigma3H, "Time" = Time_S2, "objVec" = objVec) )
# }
# 
# 
# oneReplicate_Ignore_S2 = function(seedJ) {
#   set.seed(seedJ + repID * 300)
#   source("./oneReplicate/oneReplicate-New-S2.R")
#   Estimate = batch.Ignore.S2(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2)
#   a0H = Estimate$a0
#   a0Var = Estimate$a0Var
#   a1H = NULL
#   a3H = NULL
#   betaH = Estimate$beta
#   rho1H = Estimate$rho1
#   rho2H = Estimate$rho2
#   sigma1H = Estimate$sigma1
#   sigma2H = Estimate$sigma2
#   sigma3H = Estimate$sigma3
#   objVec = Estimate$objVec
#   pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
#   Time = Estimate$Time
#   
#   return(list("a0" = a0H, "a0Var" = a0Var, "a1"=a1H, "a3"=a3H, 
#               "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H,
#               "rho1" = rho1H, "rho2" = rho2H, 
#               "beta" = betaH, "objVec" = objVec, 
#               "Time" = Time, "p.value" = pv))
# }
# 
# 
# oneReplicateWrap_Ignore_S2 = function(seedJ){
#   eval = oneReplicate_Ignore_S2(seedJ)
#   return(eval)
# }
# 
# 
# batch.Ignore.S2 = function(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
#   ind0 <- X == 0
#   ind1 <- X == 1
#   Yc1 <- Y[ind0]
#   Yt2 <- Y[ind1]
#   Zc1 <- Z[ind0, , drop = F]
#   Zt2 <- Z[ind1, , drop = F]
#   
#   Zc3 = Zc1[ind.r1, , drop = F]
#   Yc3 = Y.r1
#   Zt3 = Zt2[ind.r2, , drop = F]
#   Yt3 = Y.r2
#    
#   
#   Estimate = Estimate_Ignore_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, ind.r1, ind.r2, 
#                                 tol.c = 1e-7)
#   a0H = Estimate$a0
#   a0Var = Estimate$a0Var
#   a1H = NULL
#   a3H = NULL
#   betaH = Estimate$beta
#   rho1H = Estimate$rho1
#   rho2H = Estimate$rho2 
#   sigma1H = Estimate$sigma1
#   sigma2H = Estimate$sigma2
#   sigma3H = Estimate$sigma3
#   objVec = Estimate$objVec
#   pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
#   Time = Estimate$Time
#   return(list("a0" = a0H, "a0Var" = a0Var, "a1"=a1H, "a3"=a3H, "beta" = betaH,
#               "rho1" = rho1H, "rho2" = rho2H, "p.value" = pv, 
#               "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H, "objVec" = objVec, "Time" = Time))
# }
# 

