### Theoretical curve ####
rm(list = ls())
library(ggplot2)
library(reshape)
require(cowplot)
library(latex2exp)
library(doSNOW)
library(foreach)
source("./code/ReMeasure_S1.R")
method = "ReMeasure"
ns = c(100, 200, 400)
a1s=a1= 0.5
repID = 2
a0s = c( 0.5, 0.8)
RecoverSeq = c(0.8, 0.9, 0.95)
r2s = seq(0, 0.9, len = 10)
reptim = 10
r1s = c(0.5, 1, 2)
alpha = 0.05

ratio.a = array(NA, c(length(a0s), length(a1s), length(r1s), length(r2s), length(ns), length(RecoverSeq) ),
                  dimnames = list(TrueEffect = paste(a0s), LocationEffect = paste(a1s), ScaleEffect = paste(r1s),
                                  SNR = paste(r2s), samples = paste(ns/2), Ratio = paste(RecoverSeq)) )
  
Path = paste0("./Oracle_S1/Ratio-", method, "-repID-", repID, "-a1-", a1, "-a0-", paste(a0s, collapse = "-"), ".RData")
for (n in ns) {  
  revn1s = rev(floor( seq(5, n/2, length = 30)) )
  for (a0 in a0s) {
    for (r1 in r1s ) {
      PowerMatrix = matrix(0, length(r2s), length(revn1s) )
      ratioMatrix = matrix(0, length(r2s), length(revn1s) )
      ratioMatrix[, 1] = 1
      for (j in 1:length(r2s)) {
        r2 = r2s[j]
        for (k in 1:length(revn1s)) {
          n1 = revn1s[k]
          result1 = c()
          result2 = c()
           for (i in 1:reptim) {
             seedJ = i
             set.seed(seedJ + repID * 300)
             source("./Simulation/oneReplicate-New-S1.R")
             ind0 <- X == 0
             ind1 <- X == 1
             Yc1 <- Y[ind0]
             Yt2 <- Y[ind1]
             Zc1 <- Z[ind0, , drop = F]
             Zt2 <- Z[ind1, , drop = F]
             Zc2 = Zc1[ind.r, , drop = F]
             Yc2 = Y.r
             
             a0Var_ora = Oracle_Variance_a0(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sqrt(v1), sqrt(v2), r2, ind.r)
             C = a0/sqrt(a0Var_ora)
             C1 = qnorm(alpha/2, lower.tail = TRUE) - C
             C2 = qnorm(alpha/2, lower.tail = TRUE) + C
             Power = pnorm(C1) + pnorm(C2)
             result1 = c(result1, a0Var_ora)
             result2 = c(result2, Power)
           }
      
          PowerMatrix[j, k] = mean(result2)
          ratioMatrix[j, k] = mean(result2)/PowerMatrix[j, 1]
        }
        for (ratio in RecoverSeq) {
          ratio.a[as.character(a0), as.character(a1), as.character(r1), as.character(r2), as.character(n/2), as.character(ratio)]  = revn1s[length( which( ratioMatrix[j, ] >= ratio) )]/(n/2)
        }
      }
    }  
  }
}  

save(ratio.a, file = Path)





