#### Approximate Power Calculation ###
### This is the function for compute the power
Calculate_Power_S1 = function(nc1, nt2, r1, r2, a0, a1, nc2, seedJ = 2 ){
  set.seed(seedJ)
  n = nc1 + nt2 
  X =  as.numeric(gl(2, n / 2)) - 1
  Z <- cbind(rep(1, n), rnorm(n))
  b <- c(0, -0.5)
  
  v1 = 2/(1+r1)
  v2 = r1 * v1
  Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
  Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
  
  Z.r.a <- Z[1 : (n / 2), ]
  Et.r.a <- Et[1 : (n / 2)]
  Y.r.a <- a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) +
    rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )
  
  ind.r <- 1:nc2 
  Y.r = Y.r.a[ind.r]
  
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  
  a0Var = Oracle_Variance_a0(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sqrt(v1), sqrt(v2), r2, ind.r)
  
  C = a0/sqrt(a0Var)
  
  C1 = qnorm(alpha/2, lower.tail = TRUE) - C
  C2 = qnorm(alpha/2, lower.tail = TRUE) + C
  
  Power = pnorm(C1) + pnorm(C2)
  #qnorm(1 - alpha/2) - C
  return(Power)
}

oneReplicate_theory = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S1.R")
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  
  a0Var = Oracle_Variance_a0(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sqrt(v1), sqrt(v2), r2, ind.r)
  # C = sqrt(n)*a0/sqrt(a0Var)
  C = a0/sqrt(a0Var)
  C1 = qnorm(alpha/2, lower.tail = TRUE) - C
  C2 = qnorm(alpha/2, lower.tail = TRUE) + C
  Power = pnorm(C1) + pnorm(C2)
  return("Power" = Power)
}


oneReplicateWrap_theory = function(seedJ) {
  eval = oneReplicate_theory(seedJ)
  return(eval)
}






