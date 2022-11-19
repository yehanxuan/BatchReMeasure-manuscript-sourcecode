### real data Null hypothesis ###
rm(list = ls())
source("./code/OvarianCancer.R")
source("./code/ReMeasure_S1.R")
source("./code/Batch2Only.R")
source("./code/IgnoreBatch.R")
source("./code/LSmodel.R")

allMethods = c("ReMeasure", "Batch2", "Ignore", "LS")
Datalist = c("C1C2C4C5_RNAseq", "C1C2C4C5_Agilent", "C1C2C4C5_Agilent")

C1_Agilent_re = which(colnames(C1MES_Agilent) %in% ReMeasure_C1MES_name[1:10] == TRUE)
C1_Agilent_nonre = which(!(colnames(C1MES_Agilent) %in% ReMeasure_C1MES_name) == TRUE)[1:59]
C1_RNAseq_re = which(colnames(C1MES_RNAseq) %in% ReMeasure_C1MES_name[1:10] == TRUE)
C1_RNAseq_nonre =  which(!(colnames(C1MES_RNAseq) %in% ReMeasure_C1MES_name) == TRUE)[1:7]
sort( union(C1_Agilent_re, C1_Agilent_nonre) )
sort(union( C1_RNAseq_re, C1_RNAseq_nonre ))
C1_Agilent = C1MES_Agilent[, sort( union(C1_Agilent_re, C1_Agilent_nonre) )]
C1_RNAseq = C1MES_RNAseq[ , sort(union( C1_RNAseq_re, C1_RNAseq_nonre ))]

C2_Agilent_re = which(colnames(C2IMM_Agilent) %in% ReMeasure_C2IMM_name[1:10] == TRUE)
C2_Agilent_nonre = which( !(colnames(C2IMM_Agilent) %in% ReMeasure_C2IMM_name[1:10]) == TRUE)[1:59]
C2_RNAseq_re = which(colnames(C2IMM_RNAseq) %in% ReMeasure_C2IMM_name[1:10] == TRUE)
C2_RNAseq_nonre  = which(!(colnames(C2IMM_RNAseq) %in% ReMeasure_C2IMM_name) == TRUE)[1:7]
sort( union(C2_Agilent_re, C2_Agilent_nonre))
sort( union(C2_RNAseq_re, C2_RNAseq_nonre) )
C2_Agilent = C2IMM_Agilent[ , sort( union(C2_Agilent_re, C2_Agilent_nonre)) ]
C2_RNAseq = C2IMM_RNAseq[ , sort( union(C2_RNAseq_re, C2_RNAseq_nonre) ) ]

C4_Agilent_re =  which(colnames(C4DIF_Agilent) %in% ReMeasure_C4DIF_name[1:10] == TRUE)
C4_Agilent_nonre = which( !(colnames(C4DIF_Agilent) %in% ReMeasure_C4DIF_name[1:10]) == TRUE)[1:59]
C4_RNAseq_re = which(colnames(C4DIF_RNAseq) %in% ReMeasure_C4DIF_name[1:10] == TRUE)
C4_RNAseq_nonre = which(!(colnames(C4DIF_RNAseq) %in% ReMeasure_C4DIF_name) == TRUE)[1:7]
sort( union(C4_Agilent_re, C4_Agilent_nonre) )
sort( union(C4_RNAseq_re, C4_RNAseq_nonre))
C4_Agilent = C4DIF_Agilent[ , sort( union(C4_Agilent_re, C4_Agilent_nonre) )]
C4_RNAseq = C4DIF_RNAseq[ , sort( union(C4_RNAseq_re, C4_RNAseq_nonre))]

C5_Agilent_re = which(colnames(C5PRO_Agilent) %in% ReMeasure_C5PRO_name[1:10] == TRUE)
C5_Agilent_nonre = which( !(colnames(C5PRO_Agilent) %in% ReMeasure_C5PRO_name[1:10]) == TRUE)[1:59]
C5_RNAseq_re = which(colnames(C5PRO_RNAseq) %in% ReMeasure_C5PRO_name[1:10] == TRUE)
C5_RNAseq_nonre = which(!(colnames(C5PRO_RNAseq) %in% ReMeasure_C5PRO_name) == TRUE)[1:7]
sort( union(C5_Agilent_re, C5_Agilent_nonre) )
sort( union(C5_RNAseq_re, C5_RNAseq_nonre))
C5_Agilent = C5PRO_Agilent[ , sort( union(C5_Agilent_re, C5_Agilent_nonre) )]
C5_RNAseq = C5PRO_RNAseq[ , sort( union(C5_RNAseq_re, C5_RNAseq_nonre)) ]

Yc1mat = cbind(C1_RNAseq, C2_RNAseq, C4_RNAseq, C5_RNAseq)
Yt2mat = cbind(C1_Agilent, C2_Agilent, C4_Agilent, C5_Agilent)

Ind_all = c( rbind( ReMeasure_C1MES_name[1:10], ReMeasure_C2IMM_name[1:10],
       ReMeasure_C4DIF_name[1:10], ReMeasure_C5PRO_name[1:10]) )
n1s <- seq(5, length(Ind_all), by = 5)

library(doSNOW)
nCPUS = 10 
cl <- makeCluster(nCPUS, type = "SOCK", outfile = 'OvarianCancer_log.txt')
registerDoSNOW(cl)
library(foreach)

for ( m in 1:length(allMethods) ) {
  method = allMethods[m]
  res <- foreach(n1=iter(n1s), .combine = 'cbind') %do% {
    Ind = Ind_all[1:n1]
    Yc2mat = Agilent.dat[ ,Ind]
    
    Zc1 = matrix(1, nrow = ncol(Yc1mat), 1)
    Zt2 = matrix(1, nrow = ncol(Yt2mat), 1)
    Zc2 = matrix(1, nrow = ncol(Yc2mat), 1)
    output <- foreach (i = 1:nrow(Yc1mat), .combine = 'c') %dopar% {
      library(RConics) 
      Yc1 = Yc1mat[i, ]
      Yt2 = Yt2mat[i, ]
      Yc2 = Yc2mat[i, ]
      nc2 = length(Yc2)
      Index = which(names(Yc1) %in% Ind)
      if (method == "ReMeasure") {
        Estimate = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7)
        rhoH = Estimate$rho
        sigma1H = Estimate$sigma1
        objVec = Estimate$objVec
      }  else if (method == "Batch2") {
        Estimate = Estimate_OnlyBatch2(Zt2, Zc2, Yt2, Yc2)
        rhoH = NULL
        sigma1H = NULL
        objVec = NULL
      }  else if (method == "Ignore") {
        Estimate = Estimate_Ignore_S1(Zc1, Zt2, Yc1, Yt2)
        rhoH = NULL
        sigma1H = Estimate$sigma1
        objVec = Estimate$objVec
      } else if (method == "LS") {
        Z = rbind(Zc1, Zt2)
        Y = c(Yc1,  Yt2)
        X = c( rep(0, length(Yc1)), rep(1, length(Yt2)))
        try({ 
          Estimate = batch.correct.r1.naive1(Y, X, Z, Index, Yc2)
        })
      }
      
      a0H = Estimate$a0
      a0Var = Estimate$a0Var
      pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
      return(pv)
    }
    return(output)
  }
  colnames(res) <- n1s
  # specify the path where data is stored
  write.csv(res, file = paste0( './mixorder-', paste(Datalist, collapse = '-'), "-method-", method,
                               ".csv"), row.names = FALSE)
}

stopCluster(cl)





