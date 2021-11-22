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
  if (length(Index_C) == nc1) {
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
  
  c0 = -(1 - rho2H^2)*W3p/N3 - (1 - rho1H^2)*W4p/N3 
  c1 = (1 - rho2H^2)*rho1H*W13p/sigma1H  + (1 - rho1H^2)*rho2H*W24p/sigma2H
  c1 = c1/N3
  # c2 = 1
  c2 = (1 - rho1H^2)*(1 - rho2H^2)
  Roots = Re(polyroot(c(c0, c1, c2)))
  sigma3H = Roots[Roots > 0]
  return(sigma3H)
}

S2_Update_a0 = function(Zc1, Zt2, Zt3, Yc1, Yt2, Yt3, betaH, rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T) {
  nt3 = nrow(Zt3)
  nt2 = nrow(Zt2)
  R1p = mean(Yc1[Index_C] - Zc1[Index_C, , drop = F]%*%betaH)
  R3p = mean(Yc3 - Zc3%*%betaH)
  R2p = mean( Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
  R4p = mean(Yt3 - Zt3%*%betaH)
  if (nt3 == nt2) {
    R2s = 0
  } else {
    R2s = mean(Yt2[-Index_T] - Zt2[-Index_T, , drop = F]%*%betaH)  
  }
  
  a0H = R4p - R3p + rho1H*(sigma3H/sigma1H)*R1p - 
    rho2H*(sigma3H/sigma2H)*(1 - nt3/nt2)*(R2p - R2s)
  
  return(a0H)
}

# S2_Update_a0 = function(Zt2, Zt3, Yt2, Yt3, a1H, a3H, rho2H, sigma2H, sigma3H, Index_T) {
#   nt3 = nrow(Zt3)
#   nt2 = nrow(Zt2)
#   w1 = nt3/( (sigma2H^2)*( 1 - rho2H^2) )
#   w2 = nt3*rho2H/(sigma2H*sigma3H*(1 - rho2H^2)) 
#   w3 = nt3/( (sigma3H^2)*(1 - rho2H^2) )
#   w4 = nt3*rho2H/(sigma2H*sigma3H*( 1 - rho2H^2) )
#   w5 = (nt2 - nt3)/(sigma2H^2)
#   
#   R2p = mean( Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
#   R4p = mean(Yt3 - Zt3%*%betaH)
#   if (nt3 == nt2) {
#     R2s = 0
#   } else {
#     R2s = mean(Yt2[-Index_T] - Zt2[-Index_T, , drop = F]%*%betaH)  
#   }
#   
#   tmp = w1*(R2p - a1H) - w2*(R4p - a3H) + w3*(R4p - a3H) - w4*(R2p - a1H) + w5*(R2s - a1H)
#   a0H = tmp/(w1 - w2 + w3 - w4 + w5)
#   return(a0H)
# }

S2_Update_a1 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, rho1H, rho2H, 
                        sigma1H, sigma2H, sigma3H,
                        Index_C, Index_T) {
  nt2 = nrow(Zt2)
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  R1p = mean(Yc1[Index_C] - Zc1[Index_C, , drop = F]%*%betaH)
  R3p = mean(Yc3 - Zc3%*%betaH)
  R2p = mean(Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
  R4p = mean(Yt3 - Zt3%*%betaH)
  if (nt3 == nt2) {
    R2s = 0
  } else {
    R2s = mean(Yt2[-Index_T] - Zt2[-Index_T, , drop = F]%*%betaH)  
  }
  
  a1H = R3p - rho1H*sigma3H*R1p/sigma1H + (nt3/nt2 + (1 - nt3/nt2)*rho2H*sigma3H/sigma2H )*R2p - R4p +
    (1 - nt3/nt2)*(1 - rho2H*sigma3H/sigma2H)*R2s
  return(a1H)
}

# S2_Update_a1 = function(Zt2, Zt3, Yt2, Yt3, betaH, a0H, a3H, rho2H, sigma2H, sigma3H, Index_T) {
#   nt3 = nrow(Zt3)
#   nt2 = nrow(Zt2)
#   w1 = nt3/( (sigma2H^2)*(1 - rho2H^2) )
#   w2 = nt3*rho2H/( sigma2H*sigma3H*(1 - rho2H^2) )
#   w3 = (nt2 - nt3)/(sigma2H^2)
#   R2p = mean(Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
#   R4p = mean(Yt3 - Zt3%*%betaH)
#   if (nt3 == nt2) {
#     R2s = 0
#   } else {
#     R2s = mean(Yt2[-Index_T] - Zt2[-Index_T, , drop = F]%*%betaH)  
#   }
#   
#   a1H = ( w1*(R2p - a0H) - w2*(R4p - a0H - a3H) + w3*(R2s - a0H) )/(w1 + w3)
#   return(a1H)
# }

S2_Update_a3 = function(Zc1, Yc1, Zc3, Yc3, betaH, rho1H, sigma1H, sigma3H, Index_C){
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  R1p = mean(Yc1[Index_C] - Zc1[Index_C, , drop = F]%*%betaH)
  R3p = mean(Yc3 - Zc3%*%betaH)
  a3H = R3p - rho1H*sigma3H*R1p/sigma1H
  return(a3H)
}

# S2_Update_a3 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H, a1H, 
#                         rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T) {
#   nc3 = nrow(Zc3)
#   nt3 = nrow(Zt3)
#   
#   w1 = nc3/((sigma3H^2)*(1 - rho1H^2))
#   w2 = nc3*rho1H/(sigma1H*sigma3H*(1 - rho1H^2))
#   w3 = nt3/(sigma3H^2*(1 - rho2H^2))
#   w4 = nt3*rho2H/(sigma2H*sigma3H*(1 - rho2H^2))
#   
#   R1p = mean(Yc1[Index_C] - Zc1[Index_C, , drop = F]%*%betaH)
#   R3p = mean(Yc3 - Zc3%*%betaH)
#   R2p = mean(Yt2[Index_T] - Zt2[Index_T, , drop = F]%*%betaH)
#   R4p = mean(Yt3 - Zt3%*%betaH)
#   
#   a3H = ( w1*R3p - w2*R1p + w3*(R4p - a0H) - w4*(R2p - a0H -a1H) )/(w1 + w3)
#   return(a3H)
# }

# S2_Update_b = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, rho1H, rho2H,  )


S2_Update_b = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, rho1H, rho2H, 
                       sigma1H, sigma2H, sigma3H) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  
  Ztilde1 = (1 - rho1H*sigma3H/sigma1H)*Zc3
  Cov1 = t(Zc1)%*%Zc1/(sigma1H^2)
  Cov2 = t(Ztilde1)%*%t(t(Ztilde1) - colMeans(Ztilde1) )/( (sigma3H^2)*(1 - rho1H^2) )
  
  ratioC = nc3/nc1
  ratioT = nt3/nt2
  R1 = rho1H*sigma3H/sigma1H
  R2 = rho2H*sigma3H/sigma2H
  
  Zt3c =  rep(1, nt3)%*% t(rep(1, nt3))%*%Zt3/nt3
  Cov3 =  ( 1/(sigma2H^2*(1 - rho2H^2)) - rho2H/(sigma2H*sigma3H*(1 - rho2H^2))  ) *
    t(Zt3) %*% (Zt3 - nt3*Zt3c/nt2 ) 
  
  Cov3 = Cov3 + (1 - R2)*t(Zt3)%*%( Zt3 - Zt3c + R2*(1 - ratioT)*Zt3c )/(sigma3H^2*(1 - rho2H^2))    
      
  if (nt3 == nt2) {
    Cov4 = 0 
    Cov5 = 0
  } else {
    Zt2cs =  Zt2[-Index_T, , drop = F]
    Cov4 = t(Zt2cs)%*%Zt2cs/(sigma3H^2)  # sigma2H??
    Cov5 = -(1 - ratioT)*t(Zt3)%*%rep(1, nt3)%*%colMeans(Zt2cs)/(sigma2H^2)
  }
  
  Cov = Cov1 + Cov2 + Cov3 + Cov4 + Cov5 
  
  if (nc3 == nc1) {
    Cor1 = t(Zc3)%*%Yc1[Index_C]/(sigma1H^2) + 
    (1 - R1) * t(Zc3)%*% ( Yc3 - mean(Yc3) - R1*(Yc1[Index_C] - mean(Yc1[Index_C]) ) )/(sigma3H^2*(1 - rho1H^2))
  } else {
    Cor1 = t(Zc3)%*%Yc1[Index_C]/(sigma1H^2) + t(Zc1[-Index_C, , drop = F]) %*% Yc1[-Index_C]/(sigma1H^2) + 
    (1 - R1) * t(Zc3)%*% ( Yc3 - mean(Yc3) - R1*(Yc1[Index_C] - mean(Yc1[Index_C]) ) )/(sigma3H^2*(1 - rho1H^2))
  }
  
  
  Yt2s_mean = mean(Yt2[Index_T]) 
  if (nt3 == nt2) {
    Cor2 =  t(Zt3)%*%(Yt2[Index_T] - ratioT*Yt2s_mean ) *
      (1/(sigma2H^2*(1 - rho2H^2)) - rho2H/(sigma2H*sigma3H*(1 - rho2H^2))) 
    Cor3 = (1 - R2) * t(Zt3) %*% (Yt3 - mean(Yt3) + R2*(1 - ratioT)*Yt2s_mean )/(sigma3H^2*(1 - rho2H^2))
    Cor4 = 0
  } else {
    Yt2un_mean = mean(Yt2[-Index_T])
    Cor2 = ( t(Zt3)%*%(Yt2[Index_T] - ratioT*Yt2s_mean - (1-ratioT)*Yt2un_mean)) *
      (1/(sigma2H^2*(1 - rho2H^2)) - rho2H/(sigma2H*sigma3H*(1 - rho2H^2))) 
    Cor3 = (1 - R2) * t(Zt3) %*% (Yt3 - mean(Yt3) + R2*(1 - ratioT)*(Yt2s_mean - Yt2un_mean) )/(sigma3H^2*(1 - rho2H^2))
    Cor4 = t(Zt2[-Index_T, , drop = F])%*%Yt2[-Index_T]/(sigma3H^2)
  }
  
  
  Cor = Cor1 + Cor2 + Cor3 + Cor4
  betaH = solve(Cov, Cor)
  return(betaH)
}


# S2_Update_b = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, a0H, a1H, a3H,
#                        rho1H, rho2H, sigma1H, sigma2H, sigma3H) {
#   
#   Cov1 = t(Zc3)%*% Zc3 * ( 1/( (sigma1H^2)*(1 - rho1H^2)) + 1/( (sigma3H^2)*(1 - rho1H^2)) -
#     2*rho1H/(sigma1H*sigma3H*(1 - rho1H^2)) )
#   
#   Cov2 = t(Zc1[-Index_C, , drop = F])%*%Zc1[-Index_C, , drop = F]/(sigma1H^2)
#   Cov3 = t(Zt3)%*%Zt3*( 1/( (sigma2H^2)*(1 - rho2H^2 )) + 1/( (sigma3H^2)*(1 - rho2H^2 )) - 
#     2*rho2H/(sigma2H*sigma3H*(1-rho2H^2)) )
#   Cov4 = t(Zt2[-Index_T, , drop = F])%*%Zt2[-Index_T, , drop = F]/(sigma2H^2)
#   
#   Cov = Cov1 + Cov2 + Cov3 + Cov4
#   
#   Cor1 = t(Zc3)%*%Yc1[Index_C] * 1/( (sigma1H^2)*(1 - rho1H^2)) + 
#     t(Zc3)%*%(Yc3 - a3H) * 1/( (sigma3H^2)*(1 - rho1H^2)) - t(Zc3)%*%(Yc3 - a3H) * rho1H/(sigma1H*sigma3H*(1 - rho1H^2))-
#     t(Zc3)%*%Yc1[Index_C] * rho1H/(sigma1H*sigma3H*(1 - rho1H^2)) + 
#     t(Zc1[-Index_C, , drop = F]) %*% Yc1[-Index_C]/(sigma1H^2)
#   
#   
#   Cor2 = t(Zt3)%*%(Yt2[Index_T] - a0H - a1H) * 1/(sigma2H^2*(1 - rho2H^2 )) + 
#     t(Zt3)%*%(Yt3 - a0H - a3H) * 1/(sigma3H^2*(1 - rho2H^2 )) - 
#     t(Zt3)%*%(Yt3 - a0H - a3H) * rho2H/(sigma2H*sigma3H*(1-rho2H^2)) -
#     t(Zt3)%*%(Yt2[Index_T] - a0H - a1H) * rho2H/(sigma2H*sigma3H*(1-rho2H^2)) +
#     t(Zt2[-Index_T, , drop = F])%*%(Yt2[-Index_T] - a0H - a1H)/(sigma2H^2)
#     
#   Cor = Cor1 + Cor2
#   betaH = solve(Cov, Cor)
# 
#   return(betaH)
# }

#Negative log likelihood
ObjectiveValue_S2 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H, a1H, a3H,
                             rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T) {
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
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
  
  if (nc3 == nc1) {
    part2 = 0
  } else {
    part2 = (nc1 - nc3)*log(sigma1H) + t(Yc1Scale[-Index_C])%*%Yc1Scale[-Index_C]/2
  }
  
  
  part3 = nt3*log(sigma2H*sigma3H) + nt3*log( 1- rho2H^2)/2 + 
    ( t(Yt2Scale[Index_T])%*%Yt2Scale[Index_T] - 2*rho2H*t(Yt2Scale[Index_T])%*%Yt3Scale +
    t(Yt3Scale)%*%Yt3Scale )/ (2*(1 - rho2H^2))
  
  if (nt3 == nt2) {
    part4 = 0
  } else {
    part4 = (nt2 - nt3)*log(sigma2H) + 
      t(Yt2Scale[-Index_T]) %*% Yt2Scale[-Index_T]
  }
  
  
  obj = part1 + part2 + part3 + part4
  
  return(obj)
}

S2_Variance_a0 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, sigma1H, sigma2H, sigma3H, 
                          rho1H, rho2H, Index_C, Index_T) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc3 = nrow(Zc3)
  nt3 = nrow(Zt3)
  
  
  
  ratioC = nc3/nc1
  ratioT = nt3/nt2
  R1 = rho1H*sigma3H/sigma1H
  R2 = rho2H*sigma3H/sigma2H
  
  Ztilde1 = (1 - R1)*Zc3
  Cov1 = t(Zc1)%*%Zc1/(sigma1H^2)
  Cov2 = t(Ztilde1)%*%t(t(Ztilde1) - colMeans(Ztilde1) )/( (sigma3H^2)*(1 - rho1H^2) )
  
  # Cov3 = (1 - ratioT)/(sigma2H^2)*t(Zt3)%*%Zt3
  # Cov3 = Cov3 + ( (1 - rho2H*sigma3H/sigma2H)^2*(1 - ratioC)/(sigma3H^2*(1 - rho2H^2)) + 
  #   ratioT/(sigma2H^2*(1 - rho2H^2)) + ratioT/(sigma3H^2*(1 - rho2H^2)) - 
  #   2*rho2H*ratioT/(sigma2H*sigma3H*(1 - rho2H^2)) )*t(Zt3)%*% t( t(Zt3) - colMeans(Zt3))
  Zt3c =  rep(1, nt3)%*% t(rep(1, nt3))%*%Zt3/nt3
  Cov3 =  (1/(sigma2H^2*(1 - rho2H^2)) - rho2H/(sigma2H*sigma3H*(1 - rho2H^2))  ) *
    t(Zt3) %*% (Zt3 - nt3*Zt3c/nt2 ) 
  
  Cov3 = Cov3 + (1 - R2)*t(Zt3)%*%( Zt3 - Zt3c + R2*(1 - ratioT)*Zt3c )/(sigma3H^2*(1 - rho2H^2))    
  
  if (nt3 == nt2) {
    Cov4 = 0 
    Cov5 = 0
  } else {
    Zt2cs =  Zt2[-Index_T, , drop = F]
    Cov4 = t(Zt2cs)%*%Zt2cs/(sigma3H^2)  # sigma2H??
    Cov5 = -(1 - ratioT)*t(Zt3)%*%rep(1, nt3)%*%colMeans(Zt2cs)/(sigma2H^2)
  }
  
  S = Cov1 + Cov2 + Cov3 + Cov4 + Cov5 
  
  
  
  if (nt3 == nt2) {
    k = colMeans(Zt3) - colMeans(Zc3) + R1*colMeans(Zc3)
  } else {
    k = colMeans(Zt3) - colMeans(Zc3) + R1*colMeans(Zc3) -
      R2*(1 - ratioT)*(colMeans(Zt3) - colMeans(Zt2[-Index_T, , drop = F]))
  }
  
  A0 = t(Zc3)/(sigma2H^2) - 
    R1*(1 - R1)*t(Zc3)%*%(diag(nc3) - rep(1,nc3)%*%t(rep(1,nc3))/nc3 )/(sigma3H^2*(1 - rho1H^2))
  if (nc3 == nc1) {
    B0= 0
    b = 0
  } else {
    B0 = t(Zc1[-Index_C, , drop = F])/(sigma1H^2)
    B = solve(S, B0)
    b = -t(B)%*%k
  }
  
  
  C0 = t(Zt3)%*%(diag(nt3) - ratioT*rep(1,nt3)%*%t(rep(1,nt3))/nt3)*
    (1/(sigma2H^2*(1-rho2H^2)) - rho2H/(sigma2H*sigma3H*(1 - rho2H^2)) )
  C0 = C0 + R2*(1 - R2)*(1 - ratioT)*t(Zt3)%*%(rep(1,nt3)%*%t(rep(1,nt3))/nt3)/
    (sigma3H^2*(1 - rho2H^2))
  
  if (nt3 == nt2) {
    D0 = 0
    d = 0
  } else {
    D0 = (1 - ratioT)*t(Zt3)%*%( rep(1, nt3)%*%t(rep(1, nt2 - nt3))/(nt2 - nt3) )*
      (-1/(sigma2H^2*(1-rho2H^2)) + rho2H/(sigma2H*sigma3H*(1 - rho2H^2)))
    D0 = D0 - (1 - R2)*R2*(1 - ratioT)*t(Zt3)%*%( rep(1, nt3)%*%t(rep(1, nt2 - nt3))/(nt2 - nt3) )/
      ( sigma3H^2*(1 - rho2H^2) )
    D = solve(S, D0)
    d = R2*(1 - ratioT)*rep(1, nt2 - nt3)/(nt2 - nt3) - t(D)%*%k
  }
  
  
  E0 = (1 - R1)*t(Zc3)%*%(diag(nc3) - rep(1, nc3)%*%t(rep(1, nc3))/nc3 )/(sigma3H^2*(1 - rho1H^2))
  G0 = (1 - R2)*t(Zt3)%*%(diag(nt3) - rep(1, nt3)%*%t(rep(1, nt3))/nt3 )/(sigma3H^2*(1 - rho2H^2))
  
  A = solve(S, A0)
  C = solve(S, C0)
  E = solve(S, E0)
  G = solve(S, G0)
  
  
  
  a = (R1*rep(1, nc3)/nc3 - t(A)%*%k)
  c = -R2*(1 - ratioT)*rep(1, nt3)/nt3 - t(C)%*%k
  e = -rep(1, nc3)/nc3 - t(E)%*%k
  f = rep(1, nt3)/nt3 - t(G)%*%k
  
  Var = (sigma1H^2) * (t(a)%*%a + t(b)%*%b) + 
    (sigma2H^2) * (t(c)%*%c + t(d)%*%d ) + (sigma3H^2) * (t(e)%*%e + t(f)%*%f)+
    2*rho1H*sigma1H*sigma3H*(t(a)%*%e) + 2*rho2H*sigma2H*sigma3H*(t(c)%*%f)
  
#  Var = (sigma1^2) * (t(a)%*%a + t(b)%*%b) + 
#    (sigma2^2) * (t(c)%*%c + t(d)%*%d ) + (sigma3^2) * (t(e)%*%e + t(f)%*%f)+
#    2*rho1*sigma1*sigma3*(t(a)%*%e) + 2*rho2*sigma2*sigma3*(t(c)%*%f)
  
  Var = as.numeric(Var)
  return(Var)
}

Estimate_ReMeasure_S2 = function(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, Index_C, Index_T, 
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
  if (sigma1.Ini <= 0.001 ){
    sigma1.Ini = 0.5
  }
  if ( is.null(sigma2.Ini) ) {
    sigma2.Ini = sqrt( mean( (Yt2[Index_T] - mean(Yt2[Index_T]))^2 ) )
  }
  if (sigma2.Ini <= 0.001) {
    sigma2.Ini = 0.5
  }
  if ( is.null(sigma3.Ini)) {
    sigma3.Ini = sqrt( mean( (Yc3 - mean(Yc3))^2 ) )
  }
  if (sigma3.Ini <= 0.001 ){
    sigma3.Ini = 0.5
  }
  if ( is.null(rho1.Ini) ) {
    rho1.Ini = t(Yc1[Index_C] - mean(Yc1[Index_C]))%*%(Yc3 - mean(Yc3))/(nc3*sigma1.Ini*sigma3.Ini)
  }
  if ( rho1.Ini > 0.99 ) {
    rho1.Ini = 0.95
  } else if (rho1.Ini <= -0.99) {
    rho1.Ini = -0.95
  }
  
  if ( is.null(rho2.Ini) ) {
    rho2.Ini = t( Yt2[Index_T] - mean(Yt2[Index_T]) )%*%(Yt3 - mean(Yt3))/(nt3*sigma2.Ini* (sqrt( mean( (Yt3 - mean(Yt3))^2 ) )) )
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
  
  obj_old = ObjectiveValue_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H, a1H, a3H, rho1H, rho2H, sigma1H, sigma2H, sigma3H,
                              Index_C, Index_T)
  objVec = obj_old 
  i = 0
  gap = 1e7
  
  start_S2 = proc.time()[1]
  while ( (i < 100)&&(gap > tol.c) ) {
    i = i + 1
    sigma1H = S2_Update_sigma1(Zc1, Zc3, Yc1, Yc3, a3H, betaH, rho1H, sigma3H, Index_C)
    sigma2H = S2_Update_sigma2(Zt2, Zt3, Yt2, Yt3, a0H, a1H, a3H, betaH, rho2H, sigma3H, Index_T)
    sigma3H = S2_Update_sigma3(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, a0H, a1H, a3H, 
                               betaH, rho1H, rho2H, sigma1H, sigma2H, Index_C, Index_T)
    rho1H = S2_Update_rho1(Zc1, Zc3, Yc1, Yc3, a3H, betaH, sigma1H, sigma3H, Index_C)
    rho2H = S2_Update_rho2(Zt2, Zt3, Yt2, Yt3, a0H, a1H, a3H, betaH, sigma2H, sigma3H, Index_T)
    # betaH = S2_Update_b(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, a0H, a1H, a3H, 
    #                     rho1H, rho2H, sigma1H, sigma2H, sigma3H)
    betaH = S2_Update_b(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, rho1H, rho2H,
                        sigma1H, sigma2H, sigma3H)
    a0H = S2_Update_a0(Zc1, Zt2, Zt3, Yc1, Yt2, Yt3, betaH, rho1H, rho2H, sigma1H, sigma2H,
                       sigma3H, Index_C, Index_T)
    a1H = S2_Update_a1(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, rho1H, rho2H, sigma1H, sigma2H,
                       sigma3H, Index_C, Index_T)
    a3H = S2_Update_a3(Zc1, Yc1, Zc3, Yc3, betaH, rho1H, sigma1H, sigma3H, Index_C)
    
    obj_new = ObjectiveValue_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, betaH, a0H, a1H, a3H,
                                rho1H, rho2H, sigma1H, sigma2H, sigma3H, Index_C, Index_T)
    
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


batch.ReMeasure.S2 = function(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
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
  
  Estimate = Estimate_ReMeasure_S2(Zc1, Zt2, Zc3, Zt3, Yc1, Yt2, Yc3, Yt3, ind.r1, ind.r2,
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





