### If we ignore the batch effect ####

# Ignore_Update_rho = function(Zc1, Zc2, Yc1, Yc2, a0H, betaH, sigma1H, sigma2H, Index) {
#   nc2 = nrow(Zc2)
#   
#   mean1 = Zc1%*%betaH
#   meanCon = Zc2 %*% betaH
#   W1p = sum( (Yc1[Index] - mean1[Index])^2 )
#   W3p = sum( (Yc2 - meanCon)^2 )
#   W13p = sum( (Yc1[Index] - mean1[Index])*(Yc2 - meanCon) )
#   c3 = 1
#   c2 = -W13p/(nc2*sigma1H*sigma2H)
#   c1 = (W1p/(nc2*sigma1H*sigma1H) + W3p/(nc2*sigma2H*sigma2H) - 1)
#   c0 = -W13p/(nc2*sigma1H*sigma2H)
#   
#   rhoH = cubic(c(c3, c2, c1, c0))[1]
#   rhoH = Re(rhoH)
#   if (rhoH > 0.999) {
#     rhoH = 0.99
#   } else if (rhoH < - 0.999) {
#     rhoH = -0.99
#   }
#   return(rhoH)
# }

Ignore_Update_sigma1 = function(Zc1, Yc1, betaH) {
  nc1 = nrow(Zc1)
  sigma1H = sqrt( mean( (Yc1 - Zc1%*%betaH)^2 ) ) 
  return(sigma1H)
}


# Ignore_Update_sigma1 = function(Zc1, Zc2, Yc1, Yc2, a0H, betaH, rhoH, sigma2H, Index) {
#   nc1 = nrow(Zc1)
#   mean1 = Zc1%*%betaH
#   mean3 = Zc2%*%betaH
#   
#   W1p = sum( (Yc1[Index] - mean1[Index])^2)
#   if (length(Index) == nc1) {
#     W1s = 0
#   } else {
#     W1s = sum( (Yc1[-Index] - mean1[-Index])^2 )
#   }
#   W13p = sum( (Yc1[Index] - mean1[Index])*(Yc2 - mean3) )
#   ## In case rho is larger than 1
#   if (rhoH > 0.999) {
#     rhoH = 0.99
#   } else if (rhoH < -0.999) {
#     rhoH = -0.99
#   }
#   c0 = -W1p/nc1 - (1 - rhoH^2)*W1s/nc1
#   c1 = rhoH*W13p/(nc1*sigma2H)
#   c2 = (1 - rhoH^2)
#   Roots = Re(polyroot(c(c0, c1, c2)))
#   sigma1H = Roots[Roots > 0]
#   return(sigma1H)
# }

Ignore_Update_sigma2 = function(Zt2, Yt2, betaH){
  nt2 = nrow(Zt2)
  Zt2c = t( t(Zt2) - colMeans(Zt2) ) 
  Yt2c = Yt2 - mean(Yt2)
  sigma2H = sqrt( mean( (Yt2c - Zt2c%*%betaH)^2 ) )
  return(sigma2H)
}

# Ignore_Update_sigma2 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, 
#                                 a0H, betaH, rhoH, sigma1H, Index) {
#   nc2 = nrow(Zc2)
#   nt2 = nrow(Zt2)
#   N2 = nc2 + nt2 
#   mean1 = Zc1%*%betaH
#   mean2 = a0H + Zt2%*%betaH
#   mean3 = Zc2%*%betaH 
#   
#   W2 = sum( (Yt2 - mean2)^2 )
#   W3p = sum( (Yc2 - mean3)^2 )
#   W13p = sum( (Yc1[Index] - mean1[Index])*(Yc2 - mean3) )
#   if (rhoH > 0.999) {
#     rhoH = 0.99
#   } else if (rhoH < -0.999) {
#     rhoH = -0.99
#   }
#   c0 = -W3p/N2 - (1-rhoH^2)*W2/N2
#   c1 = rhoH*W13p/(N2*sigma1H)
#   c2 = (1 - rhoH^2)
#   Roots = Re( polyroot(c(c0, c1, c2)) )
#   sigma2H = Roots[Roots > 0]
#   return(sigma2H)
# }

Ignore_Update_beta = function(Zc1, Zt2, Yc1, Yt2, sigma1H, sigma2H) {
  nt2 = nrow(Zt2)
  Zt2c = t( t(Zt2) - colMeans(Zt2) ) 
  Yt2c = Yt2 - mean(Yt2)
  Cov = t(Zc1)%*%Zc1/sigma1H^2 + t(Zt2c)%*%Zt2c/sigma2H^2 
  Cor = t(Zc1)%*%Yc1/sigma1H^2 + t(Zt2c)%*% Yt2c/sigma2H^2
  betaH = solve(Cov, Cor)
  return(betaH)            
}

# Ignore_Update_beta = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, rhoH, sigma1H, sigma2H, Index) {
#   nc1 = nrow(Zc1)
#   Cor1 = ( 1/(sigma1H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Yc1[Index]
#   Cor2 = ( 1/(sigma2H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Yc2
#   if (length(Index) == nc1) {
#     Cor3 = 0
#   } else {
#     Cor3 = 1/(sigma1H^2)*t(Zc1[-Index, , drop = F])%*%Yc1[-Index]
#   }
#   
#   Zt2_c = t( t(Zt2) - colMeans(Zt2))
#   Yt2_c = Yt2 - mean(Yt2)
#   Cor4 = t(Zt2_c)%*%Yt2_c/(sigma2H^2)
#   Cor = Cor1 + Cor2 + Cor3 + Cor4 
#   
#   
#   if (length(Index) == nc1) {
#     Cov = ( 1/(sigma1H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Zc2 +
#       ( 1/(sigma2H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Zc2 + 
#       t(Zt2_c)%*%Zt2_c/(sigma2H^2)  
#   } else {
#     Cov = ( 1/(sigma1H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Zc2 +
#       ( 1/(sigma2H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Zc2 + 
#       t(Zc1[-Index, , drop = F])%*%Zc1[-Index, , drop = F]/(sigma1H^2) + 
#       t(Zt2_c)%*%Zt2_c/(sigma2H^2)  
#   }
#   
#   
#   betaH = solve(Cov, Cor)
#   return(betaH)
# }


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

# Ignore_Objective_S1 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, betaH, rhoH, sigma1H, sigma2H, Index) {
#   nc1 = length(Yc1)
#   nc2 = length(Yc2) 
#   nt2 = nrow(Zt2)
#   mu1 = Zc1 %*% betaH
#   mu2 = a0H + Zt2%*% betaH
#   mu3 = Zc2%*%betaH 
#   
#   Yc1Scale = (Yc1 - mu1)/sigma1H
#   Yt2Scale = (Yt2 - mu2)/sigma2H
#   Yc2Scale = (Yc2 - mu3)/sigma2H
#   
#   part1 = nc2*log(sigma1H) + nc2*log(sigma2H) + (nc2/2)*log(1 - rhoH^2) + 
#     1/(2*(1-rhoH^2)) * ( t(Yc1Scale[Index])%*%Yc1Scale[Index] - 
#                            2*rhoH*t(Yc1Scale[Index])%*%Yc2Scale + 
#                            t(Yc2Scale)%*%Yc2Scale )
#   if (length(Index) == nc1) {
#     part2 = 0
#   } else {
#     part2 = (nc1 - nc2)*log(sigma1H) + t(Yc1Scale[-Index])%*%Yc1Scale[-Index]/2 
#   }
#   
#   part3 = nt2*log(sigma2H) + t(Yt2Scale)%*%Yt2Scale/2
#   obj = part1 + part2 + part3 
#   
#   return(obj)
# }

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

# Variance_a0_Ignore = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, rhoH, sigma1H, sigma2H, Index) {
#   nc1 = nrow(Zc1)
#   nt2 = nrow(Zt2)
#   Zt2_c = t( t(Zt2) - colMeans(Zt2))
#   if (length(Index) == nc1) {
#     Cov = ( 1/(sigma1H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Zc2 +
#       ( 1/(sigma2H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Zc2 + 
#       t(Zt2_c)%*%Zt2_c/(sigma2H^2)  
#   } else {
#     Cov = ( 1/(sigma1H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Zc2 +
#       ( 1/(sigma2H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)%*%Zc2 + 
#       t(Zc1[-Index, , drop = F])%*%Zc1[-Index, , drop = F]/(sigma1H^2) + 
#       t(Zt2_c)%*%Zt2_c/(sigma2H^2)  
#   }
#   
#   A0 = ( 1/(sigma1H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)
#   B0 = ( 1/(sigma2H^2*(1-rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*t(Zc2)
#   D0 = t(Zt2_c)/(sigma2H^2)
#   
#   if (length(Index) == nc1 ){
#     C0 = 0
#     C = 0 
#   } else {
#     C0 = 1/(sigma1H^2)*t(Zc1[-Index, , drop = F])
#     C = solve(Cov, C0)
#   }
#   
#   A = solve(Cov, A0)
#   B = solve(Cov, B0)
#   D = solve(Cov, D0)
#   
#   Zt2Mean = colMeans(Zt2)
#   a = -t(Zt2Mean)%*%A
#   b = -t(Zt2Mean)%*%B
#   if (length(Index) == nc1 ) {
#     c = 0
#   } else {
#     c = -t(Zt2Mean)%*%C
#   }
#   d = t(rep(1, nt2))/nt2 - t(Zt2Mean)%*%D
#   
#   a0Var = (sigma1H^2)*( a%*%t(a) + c%*%t(c) ) + 
#     (sigma2H^2)*( b%*%t(b) + d%*%t(d) ) + 
#     (2*rhoH*sigma1H*sigma2H)*(a%*%t(b))
#   
#   return(a0Var)
# }

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

# Estimate_Ignore_S1 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7,
#                               a0.Ini = NULL, beta.Ini = NULL, rho.Ini = NULL, sigma1.Ini = NULL, sigma2.Ini = NULL){
#   nc2 = nrow(Zc2)
#   if ( is.null(sigma1.Ini) ) {
#     sigma1.Ini = sqrt( mean( (Yc1[Index] - mean(Yc1[Index]) )^2 ) )
#   } 
#   if ( is.null(sigma2.Ini) ) {
#     sigma2.Ini = sqrt( mean( (Yc2 - mean(Yc2))^2 ) )
#   }
#   if ( is.null(rho.Ini) ) {
#     # All use Control in case rho larger than 1 
#     rho.Ini = as.numeric(t( Yc1[Index] - mean(Yc1[Index]) )%*%(Yc2 - mean(Yc2))/(nc2*sigma1.Ini*sigma2.Ini))
#   }
#   if (rho.Ini > 0.999) {
#     rho.Ini = 0.99
#   } else if (rho.Ini < -0.999) {
#       rho.Ini = -0.99}  
#   if ( is.null(beta.Ini) ) {
#     beta.Ini = solve( t(Zc1)%*%Zc1, t(Zc1)%*%Yc1)
#   }
#   if ( is.null(a0.Ini) ) {
#     a0.Ini = mean(Yt2 - Zt2%*%beta.Ini)
#   }
#   
#   a0H = a0.Ini
#   betaH = beta.Ini
#   sigma1H = sigma1.Ini
#   sigma2H = sigma2.Ini
#   rhoH = rho.Ini
#   
#   obj_old = Ignore_Objective_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, 
#                                 a0H, betaH, rhoH, sigma1H, sigma2H, Index)
#   objVec = obj_old 
#   i = 0
#   gap = 1e7
#   
#   start = proc.time()[1]
#   while ( (i < 100)&&(gap > tol.c) ) {
#     i = i + 1
#     sigma1H = Ignore_Update_sigma1(Zc1, Zc2, Yc1, Yc2, a0H, betaH, rhoH, sigma2H, Index)
#     sigma2H = Ignore_Update_sigma2(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, betaH, rhoH, sigma1H, Index)
#     rhoH = Ignore_Update_rho(Zc1, Zc2, Yc1, Yc2, a0H, betaH, sigma1H, sigma2H, Index)
#     betaH = Ignore_Update_beta(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, rhoH, sigma1H, sigma2H, Index)
#     a0H = Ignore_Update_a0(Zt2, Yt2, betaH)
#     obj_new = Ignore_Objective_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, betaH, rhoH, sigma1H, sigma2H, Index)
#     gap = abs(obj_old - obj_new)
#     objVec = c(objVec, obj_new)
#     obj_old = obj_new
#   }
#   
#   Time = proc.time()[1] - start
#   a0Var = Variance_a0_Ignore(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, rhoH, sigma1H, sigma2H, Index)
#   return(list("a0" = a0H, "a0Var" = a0Var, "beta" = betaH, "rho" = rhoH, "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec,
#               "Time" = Time))
# }



oneReplicate_Ignore = function(seedJ) {
  set.seed(seedJ + repID * 300)
  #source("./oneReplicate/oneReplicate-S1.R")
  #Estimate = Estimate_Ignore_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index)
  source("./oneReplicate/oneReplicate-New-S1.R")
  # Estimate = batch.Ignore.S1(Y, X, Z, ind.r, Y.r)
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

# batch.Ignore.S1 = function(Y, X, Z, ind.r, Y.r) {
#   ind0 <- X == 0
#   ind1 <- X == 1
#   Yc1 <- Y[ind0]
#   Yt2 <- Y[ind1]
#   Zc1 <- Z[ind0, , drop = F]
#   Zt2 <- Z[ind1, , drop = F]
#   
#   Zc2 = Zc1[ind.r, , drop = F]
#   Yc2 = Y.r
#   Estimate = Estimate_Ignore_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, ind.r)
#   a0H = Estimate$a0
#   a0Var = Estimate$a0Var
#   betaH = Estimate$beta
#   rhoH = Estimate$rho
#   sigma1H = Estimate$sigma1
#   sigma2H = Estimate$sigma2
#   objVec = Estimate$objVec
#   pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
#   Time = Estimate$Time
#   return(list("a0" = a0H, "a0Var" = a0Var, "beta" = betaH, "rho" = rhoH, "p.value" = pv,
#               "sigma1" = sigma1H, "sigma2" = sigma2H, "objVec" = objVec, "Time" = Time))
# }




