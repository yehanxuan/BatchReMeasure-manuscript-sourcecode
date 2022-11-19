#### Calculate the Oracle Variance #####
rm(list = ls())
library(ggplot2)
library(reshape)
source("./code/ReMeasure_S1.R")
source("./code/Batch2Only.R")

method = "ReMeasure"
n = 100
a1 = 0.5
repID = 2
n1s <- seq(5, 50, len = 10)
a0s = c(0, 0.5, 0.8)
r1s <- c(0.5, 1, 2)
r2s <- c(0.3, 0.6, 0.9)
reptim = 1000

for ( n1 in n1s ) {
  for (a0 in a0s ) {
    for (r1 in r1s ) {
      for (r2 in r2s ){
        Path = paste0("./Oracle_S1/Oracle-S1-method-", method,"-repID-",repID,"-nc1-", n/2 ,"-nt2-", n/2,"-nc2-", n1 ,"-r1-",
                      r1,"-r2-", r2,"-a0-", a0, "-a1-", a1, ".RData")
        result1 = list()
        result2 = list()
        
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
        result_ora = list("a0Var_ora" = result1, "Power" = result2)
        save(result_ora, file = Path)
      }
    }
  }
}


