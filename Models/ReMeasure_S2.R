library(RConics) 

S2_Update_rho1 = function(Zc1, Zc3, Yc1, Yc3, a3H, betaH, 
                          sigma1H, sigma3H, Index_C){
  nc3 = nrow(Zc3)
  mean1 = Zc1%*%betaH
  mean3 = a3H + Zc3%*%betaH
  
  W1p = sum( (Yc1[Index_C] - mean1[Index_C])^2 )
  W3p = sum( (Yc3 - mean3)^2 )
  W13p = sum( (Yc1[Index_C] - mean1[Index_C])*(Yc3 - mean3) )
  # Solve quadratic equation
  c3 = 1
  c2 = -W13p/(nc3*sigma1H*sigma3H)
  c1 = W1p/(nc3*sigma1H*sigma1H) + W3p/(nc3*sigma3H*sigma3H) - 1
  c0 = -W13p/(nc3*sigma1H*sigma3H)
  
  rho1H = cubic( c(c3, c2, c1, c0) )[1]
  rho1H = Re(rho1H) 
  if (rho1H > 0.999) {
    rho1H = 0.99
  } else if (rho1H < -0.999 ) {
    rho1H = -0.99
  }
  return(rho1H)
}
 
S2_Update_rho2 = function(Zt2, Zt3, Yt2, Yt3, a0H, a1H, a3H, betaH, sigma2H, sigma3H, Index_T) {
  nt3 = nrow(Zt3)
  mean2 = a0H + a1H + Zt2 %*% betaH
  mean4 = a0H + a3H + Zt3 %*% betaH
  
  W2p = sum( (Yt2[Index_T] - mean2[Index_T])^2 )
  W4p = sum( (Yt3 - mean4)^2 )
  W24p = sum( (Yt2[Index_T] - mean2[Index_T])*(Yt3 - mean4) )
  
  c3 = 1
  c2 = -W24p/(nt3*sigma2H*sigma3H)
  c1 = W2p/(nt3*sigma2H*sigma2H) + W4p/(nt3*sigma3H*sigma3H) - 1
  c0 = -W24p/(nt3*sigma2H*sigma3H)
  
  rho2H = cubic(c(c3, c2, c1, c0))[1]
  rho2H = Re(rho2H)
  
  if (rho2H > 0.999) {
    rho2H = 0.99
  } else if (rho2H < -0.99) {
    rho2H = -0.99
  }
  return(rho2H)
}


S2_Update_sigma1 = function(Zc1, Zc3, Yc1, Yc3, a3H, betaH, rho1H, sigma3H, Index_C) {
  nc1 = nrow(Zc1)
  mean1 = Zc1%*%betaH
  mean3 = a3H + Zc3%*%betaH
  
  W1p = sum( (Yc1[Index_C] - mean1[Index_C])^2)
  if (length(Index) == nc1) {
    W1s = 0
  } else {
    W1s = sum( (Yc1[-Index_C] - mean1[-Index_C])^2 )
  }
 
  W13p = sum( (Yc1[Index_C] - mean1[Index_C])*(Yc3 - mean3) )
  if (rho1H > 0.999) {
    rho1H = 0.99
  } else if (rho1H < -0.999) {
    rho1H = -0.99
  }
  c0 = -W1p/nc1 - (1 - rho1H^2)*W1s/nc1
  c1 = rho1H*W13p/(nc1*sigma3H)
  c2 = (1 - rho1H^2)
  Roots = Re(polyroot(c(c0, c1, c2)))
  sigma1H = Roots[Roots > 0 ]
  return(sigma1H)
}


S2_Update_sigma2 = function(Zt2, Zt3, Yt2, Yt3, a0H, a1H, a3H, betaH, rho2H, sigma3H, Index_T) {
  nt2 = nrow(Zt2)
  mean2 = a0H + a1H + Zt2%*%betaH
  mean4 = a0H + a3H + Zt3%*%betaH
  
  W2p = sum( (Yt2[Index_T] - mean2[Index_T])^2 )
  if (length(Index_T) == nt2) {
    W2s = 0
  } else {
    W2s = sum( (Yt2[-Index_T] - mean2[-Index_T])^2 )  
  }
  W24p = sum( (Yt2[Index_T] - mean2[Index_T])*(Yt3 - mean4) )
  if (rho2H > 0.999) {
    rho2H = 0.99
  } else if (rho2H < -0.999) {
    rho2H = -0.99
  }
  c0 = -W2p/nt2 - (1 - rho2H^2)*W2s/nt2
  c1 = rho2H*W24p/(nt2*sigma3H)
  c2 = (1 - rho2H^2)
  Roots = Re(polyroot(c(c0, c1, c2)))
  sigma2H = Roots[Roots > 0]
  return(sigma2H)
}


S2_Update_sigma3 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, a0H, a1H, a3H,
                            betaH, rho1H, rho2H, sigma1H, sigma2H, Index_C, Index_T) {
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  N3 = nc3 + nt3
  
  mean1 = Zc1%*%betaH
  mean2 = a0H + a1H + Zt2%*%betaH
  mean3 = a3H + Zc3%*%betaH
  mean4 = a0H + a3H + Zt3%*%betaH
  
  W3p = sum( (Yc3 - mean3)^2)
  W4p = sum( (Yt3 - mean4)^2)
  W13p = sum( (Yc1[Index_C] - mean1[Index_C])*(Yc3 - mean3) )
  W24p = sum( (Yt2[Index_T] - mean2[Index_T])*(Yt3 - mean4) ) 
  
  c0 = -W3p/( N3*(1 - rho1H^2) ) - W4p/(N3*(1 - rho2H^2)) 
  c1 = rho1H*W13p/( sigma1H*( 1 - rho1H^2) ) + rho2H*W24p/(sigma2H*(1 - rho2H^2))
  c1 = c1/N3
  c2 = 1
  Roots = Re(polyroot(c(c0, c1, c2)))
  sigma3H = Roots[Roots > 0]
  return(sigma3H)
}


S2_Update_a0 = function(Zt2, Zt3, Yt2, Yt3, a1H, a3H, rho2H, sigma2H, sigma3H, Index_T) {
  nt3 = nrow(Zt3)
  nt2 = nrow(Zt2)
  w1 = nt3/( (sigma2H^2)*( 1 - rho2H^2) )
  w2 = nt3*rho2H/(sigma2H*sigma3H*(1 - rho2H^2)) 
  w3 = nt3/( (sigma3H^2)*(1 - rho2H^2) )
  w4 = nt3*rho2H/(sigma2H*sigma3H*( 1 - rho2H^2) )
  w5 = (nt2 - nt3)/(sigma2H^2)
  
  R2p = mean( Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
  R4p = mean(Yt3 - Zt3%*%betaH)
  if (nt3 == nt2) {
    R2s = 0
  } else {
    R2s = mean(Yt2[-Index_T] - Zt2[-Index_T, , drop = F]%*%betaH)  
  }
  
  tmp = w1*(R2p - a1H) - w2*(R4p - a3H) + w3*(R4p - a3H) - w4*(R2p - a1H) + w5*(R2s - a1H)
  a0H = tmp/(w1 - w2 + w3 - w4 + w5)
  return(a0H)
}


S2_Update_a1 = function(Zt2, Zt3, Yt2, Yt3, betaH, a0H, a3H, rho2H, sigma2H, sigma3H, Index_T) {
  nt3 = nrow(Zt3)
  nt2 = nrow(Zt2)
  w1 = nt3/( (sigma2H^2)*(1 - rho2H^2) )
  w2 = nt3*rho2H/( sigma2H*sigma3H*(1 - rho2H^2) )
  w3 = (nt2 - nt3)/(sigma2H^2)
  R2p = mean(Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
  R4p = mean(Yt3 - Zt3%*%betaH)
  if (nt3 == nt2) {
    R2s = 0
  } else {
    R2s = mean(Yt2[-Index_T] - Zt2[-Index_T, , drop = F]%*%betaH)  
  }
  
  a1H = ( w1*(R2p - a0H) - w2*(R4p - a0H - a3H) + w3*(R2s - a0H) )/(w1 + w3)
  return(a1H)
}

S2_Update_a3 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H, a1H, 
                        rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T) {
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  
  w1 = nc3/((sigma3H^2)*(1 - rho1H^2))
  w2 = nc3*rho1H/(sigma1H*sigma3H*(1 - rho1H^2))
  w3 = nt3/(sigma3H^2*(1 - rho2H^2))
  w4 = nt3*rho2H/(sigma2H*sigma3H*(1 - rho2H^2))
  
  R1p = mean(Yc1[Index_C] - Zc1[Index_C, , drop = F]%*%betaH)
  R3p = mean(Yc3 - Zc3%*%betaH)
  R2p = mean(Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
  R4p = mean(Yt3 - Zt3%*%betaH)
  
  a3H = ( w1*R3p - w2*R1p + w3*(R4p - a0H) - w4*(R2p - a0H -a1H) )/(w1 + w3)
  return(a3H)
}


S2_Update_b = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, a0H, a1H, a3H,
                       rho1H, rho2H, sigma1H, sigma2H, sigma3H) {
  
  Cov1 = t(Zc3)%*% Zc3 * ( 1/( (sigma1H^2)*(1 - rho1H^2)) + 1/( (sigma3H^2)*(1 - rho1H^2)) -
    2*rho1H/(sigma1H*sigma3H*(1 - rho1H^2)) )
  
  Cov2 = t(Zc1[-Index_C, , drop = F])%*%Zc1[-Index_C, , drop = F]/(sigma1H^2)
  Cov3 = t(Zt3)%*%Zt3*( 1/( (sigma2H^2)*(1 - rho2H^2 )) + 1/( (sigma3H^2)*(1 - rho2H^2 )) - 
    2*rho2H/(sigma2H*sigma3H*(1-rho2H^2)) )
  Cov4 = t(Zt2[-Index_T, , drop = F])%*%Zt2[-Index_T, , drop = F]/(sigma2H^2)
  
  Cov = Cov1 + Cov2 + Cov3 + Cov4
  
  Cor1 = t(Zc3)%*%Yc1[Index_C] * 1/( (sigma1H^2)*(1 - rho1H^2)) + 
    t(Zc3)%*%(Yc3 - a3H) * 1/( (sigma3H^2)*(1 - rho1H^2)) - t(Zc3)%*%(Yc3 - a3H) * rho1H/(sigma1H*sigma3H*(1 - rho1H^2))-
    t(Zc3)%*%Yc1[Index_C] * rho1H/(sigma1H*sigma3H*(1 - rho1H^2)) + 
    t(Zc1[-Index_C, , drop = F]) %*% Yc1[-Index_C]/(sigma1H^2)
  
  
  Cor2 = t(Zt3)%*%(Yt2[Index_T] - a0H - a1H) * 1/(sigma2H^2*(1 - rho2H^2 )) + 
    t(Zt3)%*%(Yt3 - a0H - a3H) * 1/(sigma3H^2*(1 - rho2H^2 )) - 
    t(Zt3)%*%(Yt3 - a0H - a3H) * rho2H/(sigma2H*sigma3H*(1-rho2H^2)) -
    t(Zt3)%*%(Yt2[Index_T] - a0H - a1H) * rho2H/(sigma2H*sigma3H*(1-rho2H^2)) +
    t(Zt2[-Index_T, , drop = F])%*%(Yt2[-Index_T] - a0H - a1H)/(sigma2H^2)
    
  Cor = Cor1 + Cor2
  betaH = solve(Cov, Cor)

  return(betaH)
}

#Negative log likelihood
ObjectiveValue_S2 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H, a1H, a3H,
                             rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T) {
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  mu1 = Zc1 %*% betaH
  mu2 = a0H + a1H + Zt2%*%betaH
  mu3 = a3H + Zc3%*%betaH
  mu4 = a0H + a3H + Zt3%*%betaH
  
  Yc1Scale = (Yc1 - mu1)/sigma1H
  Yt2Scale = (Yt2 - mu2)/sigma2H
  Yc3Scale = (Yc3 - mu3)/sigma3H
  Yt3Scale = (Yt3 - mu4)/sigma3H
  
  part1 = nc3*log(sigma1H*sigma3H) + nc3*log( 1 - rho1H^2)/2 + 
    ( t(Yc1Scale[Index_C])%*%Yc1Scale[Index_C] - 2*rho1H*t(Yc1Scale[Index_C])%*%Yc3Scale +
    t(Yc3Scale)%*%Yc3Scale ) / (2*(1-rho1H^2))
  
  part2 = (nc1 - nc3)*log(sigma1H) + t(Yc1Scale[-Index_C])%*%Yc1Scale[-Index_C]/2
  
  part3 = nt3*log(sigma2H*sigma3H) + nt3*log( 1- rho2H^2)/2 + 
    ( t(Yt2Scale[Index_T])%*%Yt2Scale[Index_T] - 2*rho2H*t(Yt2Scale[Index_T])%*%Yt3Scale +
    t(Yt3Scale)%*%Yt3Scale )/ (2*(1 - rho2H^2))
  
  part4 = (nt2 - nt3)*log(sigma2H) + 
    t(Yt2Scale[-Index_T]) %*% Yt2Scale[-Index_T]
  
  obj = part1 + part2 + part3 + part4
  
  return(obj)
}

S2_Variance_a0 = function() {
  
}


