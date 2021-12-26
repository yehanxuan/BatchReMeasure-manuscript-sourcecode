# The goal of L/S is making alignment of control sample in batch1 and batch2?
batch.correct.r1.naive1 <- function(Y, X, Z, ind.r, Y.r) {
  # location/scale model 
  start_LS = proc.time()[1]
  Y.o <- Y[ind.r]
  
  I <- mean(Y.r - Y.o)
  s <- sd(Y.r)/sd(Y.o)

  Y[X == 0] <- (Y[X==0] + I)
  m <- mean(Y[X==0])
  Y[X == 0] <- ( Y[X == 0] - m )*s + m 
  
  lm.obj <- lm(Y ~ X + Z - 1)  ### omit the intercept
  # after adjustion, the coeeficient the Z1 becomes a1
  pv <- coef(summary(lm.obj))[1, 4]
  
  a0H = coef(summary(lm.obj))[1, 1]
  a0Var =  coef(summary(lm.obj))[1, 2]^2
  a1H = coef(summary(lm.obj))[2, 1]
  betaH = c(0, coef(summary(lm.obj))[3, 1]) # seems non-identifiable
  rhoH = NULL
  sigma1H = sd(Y.o)
  sigma2H = sd(Y.r)
  objVec = NULL
  Time_LS = proc.time()[1] - start_LS
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "beta" = betaH, "rho" = rhoH, 
              "sigma1" = sigma1H, "sigma2" = sigma2H, "Time" = Time_LS, "objVec" = objVec, 
              "p.value" = pv))
}

oneReplicate_LS = function(seedJ) {
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S1.R")
  Estimate = batch.correct.r1.naive1(Y, X, Z, ind.r, Y.r)
  a0H = Estimate$a0  ## Fill it up later 
  a0Var = Estimate$a0var 
  a1H = Estimate$a1
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time
  pv = Estimate$p.value
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "rho" = rhoH, 
              "beta" = betaH,"objVec" = objVec, "Time" = Time, "p.value" = pv))
}

oneReplicateWrap_LS = function(seedJ) {
  eval = oneReplicate_LS(seedJ)
  return(eval)
}

