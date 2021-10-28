## MLE estimator for S1.
library(RConics) # Solve cubic equation

Update_rho = function(Zc1, Zc2, Yc1, Yc2, a0H, a1H, betaH, sigma1H, sigma2H, Index) {
   nc2 = nrow(Zc2)
   # Using estimated mean 
   mean1 = Zc1%*%betaH 
   mean3 = a1H + Zc2%*%betaH
   
   W1p = sum( (Yc1[Index] - mean1[Index])^2)
   W3p = sum( (Yc2 - mean3)^2 )
   W13p = sum( (Yc1[Index] - mean1[Index])*(Yc2 - mean3) )
   # Solve the cubic equation c0 + c1\rho + c2 \rho^2 + c3 \rho^3
   c3 = 1
   c2 = -W13p/(nc2*sigma1H*sigma2H)
   c1 = (W1p/(nc2*sigma1H*sigma1H) +  W3p/(nc2*sigma2H*sigma2H) - 1)
   c0 = -W13p/(nc2*sigma1H*sigma2H)
   
   rhoH = cubic(c(c3, c2, c1, c0))[1]
   rhoH = Re(rhoH)
   if (rhoH > 0.999) {
     rhoH = 0.99
   } else if (rhoH < -0.999) {
     rhoH = -0.99
   }
   # or we can try polyroot J, T algorithm
   #polyroot(c(c0, c1, c2, c3))
   
   return(rhoH)
}

Update_sigma1 = function(Zc1, Zc2, Yc1, Yc2, a0H, a1H, betaH, rhoH, sigma2H, Index) {
  #nc2 = nrow(Zc2) # True?
  nc1 = nrow(Zc1)
  mean1 = Zc1%*%betaH 
  mean3 = a1H + Zc2%*%betaH
  
  W1p = sum( (Yc1[Index] - mean1[Index])^2)
  if (length(Index) == nc1) {
    W1s = 0
  } else {
    W1s = sum( (Yc1[-Index] - mean1[-Index])^2 )
  }
  W13p = sum( (Yc1[Index] - mean1[Index])*(Yc2 - mean3) )
  ## In case rho is larger than 1
  if (rhoH > 0.999) {
    rhoH = 0.99
  } else if (rhoH < -0.999) {
    rhoH = -0.99
  }
  c0 = -W1p/nc1 - (1 - rhoH^2)*W1s/nc1
  c1 = rhoH*W13p/(nc1*sigma2H)
  c2 = (1 - rhoH^2)
  Roots = Re(polyroot(c(c0, c1, c2)))
  sigma1H = Roots[Roots > 0] 
  return(sigma1H)
}

Update_sigma2 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH, sigma1H, Index) {
  nc2 = nrow(Zc2)
  nt2 = nrow(Zt2)
  N2 = nc2 + nt2
  mean1 = Zc1%*%betaH
  mean2 = a0H + a1H + Zt2%*%betaH
  mean3 = a1H + Zc2%*%betaH
  
  W2 = sum( (Yt2 - mean2)^2 )
  W3p = sum( (Yc2 - mean3)^2 )
  W13p = sum( (Yc1[Index] - mean1[Index])*(Yc2 - mean3) )
  # In case rho is larger than 1
  if (rhoH > 0.999) {
    rhoH = 0.99
  } else if (rhoH < -0.999) {
    rhoH = -0.99
  }
  c0 = -W3p/N2 - (1-rhoH^2)*W2/N2
  c1 = rhoH*W13p/(N2*sigma1H)
  c2 = (1 - rhoH^2)
  Roots = Re(polyroot(c(c0, c1, c2)))
  sigma2H = Roots[Roots > 0]
  return(sigma2H)
}

Update_beta = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, rhoH, sigma1H, sigma2H, Index) {
  Cor1 = t(Zc1)%*%Yc1/(sigma1H^2)
  Cor2 = t(Zt2)%*%( Yt2 - mean(Yt2) )/(sigma2H^2)
  ratio = sigma2H/sigma1H
  Ztilde = Zc2 - rhoH*ratio*Zc1[Index, , drop = F]
  Ytilde = Yc2 - rhoH*ratio*Yc1[Index]
  
  Cor3 = 1/( (sigma2H^2)*(1 - rhoH^2) ) * t(Ztilde)%*%(Ytilde - mean(Ytilde))
  t = Cor1 + Cor2 + Cor3
  
  Cov1 = t(Zc1)%*%Zc1/(sigma1H^2) 
  Cov2 = t(Zt2)%*%t( t(Zt2) - colMeans(Zt2) )/(sigma2H^2)
  Cov3 = t(Ztilde)%*%t( t(Ztilde) -colMeans(Ztilde) )/( (sigma2H^2)*(1-rhoH^2) ) # Very small
  S = Cov1 + Cov2 + Cov3
  betaH = solve(S, t)
  return(betaH)
}

Update_a = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH, sigma1H, sigma2H, Index) {
  ratio = sigma2H/sigma1H
  Ytilde = Yc2 - rhoH*ratio*Yc1[Index]
  Ztilde = Zc2 - rhoH*ratio*Zc1[Index, , drop = F]
  
  a1H = mean(Ytilde - Ztilde%*%betaH)
  a0H = mean(Yt2 - Zt2%*%betaH) - a1H
  return(list("a1" = a1H, "a0" = a0H))
}

## Negative log likelihood for S1 # where is a0H, a1H ??
ObjectiveValue_S1 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH, sigma1H, sigma2H, Index) {
  nc1 = nrow(Zc1)
  nc2 = length(Yc2)
  nt2 = nrow(Zt2)
  mu1 = Zc1 %*% betaH
  mu2 = a0H + a1H + Zt2%*%betaH
  mu3 = a1H + Zc2%*%betaH
  
  Yc1Scale = (Yc1 - mu1)/sigma1H
  Yt2Scale = (Yt2 - mu2)/sigma2H
  Yc2Scale = (Yc2 - mu3)/sigma2H
  
  part1 = nc2*log(sigma1H) + nc2*log(sigma2H) + nc2*log(1 - rhoH^2)/2 + 
    1/(2*(1-rhoH^2))*( t(Yc1Scale[Index])%*%Yc1Scale[Index] - 
                        2*rhoH*t(Yc1Scale[Index])%*%Yc2Scale + 
                        t(Yc2Scale)%*%Yc2Scale )
  if (length(Index) == nc1) {
    part2 = 0 
  } else {
    part2 = (nc1 - nc2)*log(sigma1H) + t(Yc1Scale[-Index])%*%Yc1Scale[-Index]/2 
  }
  
  part3 = nt2*log(sigma2H) + t(Yt2Scale)%*%Yt2Scale/2
  
  obj = part1 + part2 + part3
  return(obj)
}


Variance_a0 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sigma1H, sigma2H, rhoH, Index) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  Ztilde = (1 - rhoH*(sigma2H/sigma1H) )*Zc2
  Cov1 = t(Zc1)%*%Zc1/(sigma1H^2) 
  Cov2 = t(Zt2)%*%t( t(Zt2) - colMeans(Zt2) )/(sigma2H^2)
  Cov3 = t(Ztilde)%*%t( t(Ztilde) -colMeans(Ztilde) )/( (sigma2H^2)*(1-rhoH^2) ) # Very small
  S = Cov1 + Cov2 + Cov3
  
  # Identity matrix, do not make mistake
 A0 = t(Zc2)%*%( diag(nc2)/((sigma1H^2)*(1 - rhoH^2)) - 
   rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) * (diag(nc2) - rep(1, nc2)%*%t(rep(1, nc2))/nc2 ) -
   (rhoH^2)/( (sigma1H^2)*( 1 - rhoH^2))*(rep(1, nc2)%*%t(rep(1, nc2))/nc2 ) )
  
  
  B0 = (1/(sigma2H^2) )*t(Zt2)%*%(diag(nt2) - rep(1, nt2)%*%t(rep(1, nt2))/nt2)
  
  C0 = ( 1/((sigma2H^2)*(1 - rhoH^2)) - rhoH/(sigma1H*sigma2H*(1 - rhoH^2)) )*
    t(Zc2)%*%(diag(nc2) - rep(1, nc2)%*%t(rep(1, nc2))/nc2 )
  
  Zt2Mean = colMeans(Zt2)
  ZtildeMean = colMeans(Ztilde)
  
  if (length(Index) == nc1){
    D0 = 0
    d = 0 
  } else {
    D0 = t(Zc1[-Index, , drop = F])/(sigma1H^2)
    D = solve(S, D0)
    d = t(D)%*%(ZtildeMean - Zt2Mean)
  }
  
  A = solve(S, A0)
  B = solve(S, B0)
  C = solve(S, C0)
  
 # (ZtildeMean - Zt2Mean)%*%A + (rhoH*sigma2H/sigma1H)*t(rep(1, nc2)/nc2)
              
  a = t(A)%*%(ZtildeMean - Zt2Mean) + (rhoH*sigma2H/sigma1H)*rep(1, nc2)/nc2
  b = t(B)%*%(ZtildeMean - Zt2Mean) + rep(1, nt2)/nt2
 # c = t(C)%*%(ZtildeMean - Zt2Mean) + rep(1, nc2)/nc2 ## modified here 10/12/2021
  c = t(C)%*%(ZtildeMean - Zt2Mean) - rep(1, nc2)/nc2
  
  Var = (sigma1H^2)*( t(a)%*%a + t(d)%*%d ) + 
    (sigma2H^2)*(t(b)%*%b + t(c)%*%c) + 2*rhoH*sigma1H*sigma2H*t(a)%*%c
  
  ### Hindsight variance 
  #Var = (sigma1^2)*( t(a)%*%a + t(d)%*%d ) + 
  #  (sigma2^2)*(t(b)%*%b + t(c)%*%c) + 2*rho*sigma1*sigma2*t(a)%*%c
  
  Var = as.numeric(Var)
  return(Var)
}


Estimate_ReMeasure_S1 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7, a0.Ini = NULL, a1.Ini=NULL, 
                                 rho.Ini=NULL, beta.Ini = NULL, sigma1.Ini=NULL, sigma2.Ini=NULL){
  nc2 = nrow(Zc2)
  if ( is.null(sigma1.Ini) ) {
    sigma1.Ini = sqrt( mean( (Yc1[Index] - mean(Yc1[Index]) )^2 ) )
  } 
  if ( is.null(sigma2.Ini) ) {
    sigma2.Ini = sqrt( mean( (Yc2 - mean(Yc2))^2 ) )
  }
  if ( is.null(rho.Ini) ) {
    # All use Control in case rho larger than 1 
    rho.Ini = as.numeric( t(  (Yc1[Index] - mean(Yc1[Index])) ) %*%(Yc2 - mean(Yc2))/(nc2*sigma1.Ini*sigma2.Ini))
  }
  
  if ( rho.Ini > 0.99 ) {
    rho.Ini = 0.95
  } else if (rho.Ini <= -0.99) {
    rho.Ini = -0.95
  }
  if ( is.null(beta.Ini) ) {
    beta.Ini = solve( t(Zc1)%*%Zc1, t(Zc1)%*%Yc1)
  }
  if ( is.null(a1.Ini) ) {
    a1.Ini = mean(Yc2 - Zc2%*%beta.Ini) - rho.Ini*sigma2.Ini/sigma1.Ini*mean(Yc1[Index] - Zc2%*%beta.Ini)
  } 
  if (is.null(a0.Ini)) {
    a0.Ini = mean(Yt2 - Zt2%*%beta.Ini) - a1.Ini
  }
  
  a0H = a0.Ini
  a1H = a1.Ini
  betaH = beta.Ini
  sigma1H = sigma1.Ini
  sigma2H = sigma2.Ini
  rhoH = rho.Ini
  
  obj_old = ObjectiveValue_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, 
                              betaH, rhoH, sigma1H, sigma2H, Index)
  objVec = obj_old 
  i = 0
  gap = 1e7
  
  start = proc.time()[1]
  while ( (i < 100)&&(gap > tol.c) ) {
    i = i + 1
    sigma1H = Update_sigma1(Zc1, Zc2, Yc1, Yc2, a0H, a1H, betaH, rhoH, sigma2H, Index)
    sigma2H = Update_sigma2(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH,
                            sigma1H, Index)
    rhoH = Update_rho(Zc1, Zc2, Yc1, Yc2, a0H, a1H, betaH, sigma1H, sigma2H, Index)
    
    betaH = Update_beta(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, rhoH, sigma1H, sigma2H, Index)
    aVec = Update_a(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH, sigma1H, sigma2H, Index)
    a0H = aVec$a0
    a1H = aVec$a1
    obj_new = ObjectiveValue_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, a0H, a1H, betaH, rhoH,
                                sigma1H, sigma2H, Index)
    gap = abs(obj_old - obj_new)
    objVec = c(objVec, obj_new)
    obj_old = obj_new 
  }
  
  Time = proc.time()[1] - start

  a0Var = Variance_a0(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sigma1H, sigma2H, rhoH, Index)
  
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time, "objVec" = objVec))
}



oneReplicate_ReMeasure = function(seedJ) {
  set.seed(seedJ + repID * 300)
  #source("./oneReplicate/oneReplicate-S1.R")
 #Estimate = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index)
  source("./oneReplicate/oneReplicate-New-S1.R")
  Estimate = batch.ReMeasure.S1(Y, X, Z, ind.r, Y.r)
  a0H = Estimate$a0
  a0Var = Estimate$a0Var
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time 
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "rho" = rhoH, 
              "beta" = betaH,"objVec" = objVec, "Time" = Time))
}


oneReplicateWrap_ReMeasure = function(seedJ){
  eval = oneReplicate_ReMeasure(seedJ)
  return(eval)
}









