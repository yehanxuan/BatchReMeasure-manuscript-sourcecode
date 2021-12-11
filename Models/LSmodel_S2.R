# Matching batch1 and batch3, batch2 and batch 3 
batch.correct.r2.naive2 = function(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
  start_LS = proc.time()[1]
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
  
  Y.o1 = Yc1[ind.r1]
  Y.o2 = Yt2[ind.r2]
  
  I1 = mean(Y.r1 - Y.o1)  ## control
  I2 = mean(Y.r2 - Y.o2)
  s1 <- sd(Y.r1)/sd(Y.o1) ## 
  s2 <- sd(Y.r2)/sd(Y.o2)
  
  Y[X == 0] <- (Y[X==0] + I1)
  Y[X == 1] <- (Y[X==1] + I2)
  m1 <- mean(Y[X==0])
  m2 <- mean(Y[X==1])
  Y[X == 0] <- ( Y[X == 0] - m1 )*s1 + m1
  Y[X == 1] <- ( Y[X == 1] - m2)*s2 + m2
  
  lm.obj <- lm(Y ~ X + Z - 1)
  pv <- coef(summary(lm.obj))[1, 4]
  
  a0H = coef(summary(lm.obj))[1, 1]
  a0Var =  coef(summary(lm.obj))[1, 2]^2
  a1H = NULL
  a3H =  coef(summary(lm.obj))[2, 1]
  betaH = c(0,  coef(summary(lm.obj))[3,1])
  sigma1H = sd(Y.o1)
  sigma2H = sd(Y.o2)
  nc3 = length(ind.r1)
  nt3 = length(ind.r2)
  sigma3H = ( nc3*sd(Y.r1) + nt3*sd(Y.r2) )/(nc3 + nt3)
  rho1H = NULL
  rho2H = NULL
  objVec = NULL
  Time_LS = proc.time()[1] - start_LS
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H,
              "beta" = betaH, "rho1" = rho1H, "rho2" = rho2H, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H,
              "Time" = Time_LS, "objVec" = objVec, 
              "p.value" = pv))
}


oneReplicate_LS_S2 = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S2.R")
  Estimate = batch.correct.r2.naive2(Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2)
  a0H = Estimate$a0  ## Fill it up later 
  a0Var = Estimate$a0var 
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
  pv = Estimate$p.value
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "a3" = a3H,
              "sigma1" = sigma1H, "sigma2" = sigma2H, "sigma3" = sigma3H,
              "rho1" = rho1H, "rho2" = rho2H, 
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "p.value" = pv))
}

oneReplicateWrap_LS_S2 = function(seedJ){
  eval = oneReplicate_LS_S2(seedJ)
  return(eval)
}