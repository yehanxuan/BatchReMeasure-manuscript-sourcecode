## Main NULL ##
rm(list = ls())
library(sva)
library(pamr)
library(limma)
library(stringi)
# This file preprocess the ovarian cancer data for us 
source("./code/OvarianCancer.R")
source("./code/ReMeasure_S1.R")

combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}
folder = "Data_real"
if (file.exists(folder)) {
  cat("The folder already exists") 
} else {
  dir.create(folder)
}

Datalist = c("C1C2C4C5_RNAseq", "C1C2C4C5_Agilent", "C1C2C4C5_Agilent")

C1_Agilent_re = which(colnames(C1MES_Agilent) %in% ReMeasure_C1MES_name[1:10] == TRUE)
C1_Agilent_nonre = which(!(colnames(C1MES_Agilent) %in% ReMeasure_C1MES_name) == TRUE)[1:59]
C1_RNAseq_re = which(colnames(C1MES_RNAseq) %in% ReMeasure_C1MES_name[1:10] == TRUE)
C1_RNAseq_nonre =  which(!(colnames(C1MES_RNAseq) %in% ReMeasure_C1MES_name) == TRUE)[1:7]
C1_Agilent = C1MES_Agilent[, sort( union(C1_Agilent_re, C1_Agilent_nonre) )]
C1_RNAseq = C1MES_RNAseq[ , sort(union( C1_RNAseq_re, C1_RNAseq_nonre ))]

C2_Agilent_re = which(colnames(C2IMM_Agilent) %in% ReMeasure_C2IMM_name[1:10] == TRUE)
C2_Agilent_nonre = which( !(colnames(C2IMM_Agilent) %in% ReMeasure_C2IMM_name[1:10]) == TRUE)[1:59]
C2_RNAseq_re = which(colnames(C2IMM_RNAseq) %in% ReMeasure_C2IMM_name[1:10] == TRUE)
C2_RNAseq_nonre  = which(!(colnames(C2IMM_RNAseq) %in% ReMeasure_C2IMM_name) == TRUE)[1:7]
C2_Agilent = C2IMM_Agilent[ , sort( union(C2_Agilent_re, C2_Agilent_nonre)) ]
C2_RNAseq = C2IMM_RNAseq[ , sort( union(C2_RNAseq_re, C2_RNAseq_nonre) ) ]

C4_Agilent_re =  which(colnames(C4DIF_Agilent) %in% ReMeasure_C4DIF_name[1:10] == TRUE)
C4_Agilent_nonre = which( !(colnames(C4DIF_Agilent) %in% ReMeasure_C4DIF_name[1:10]) == TRUE)[1:59]
C4_RNAseq_re = which(colnames(C4DIF_RNAseq) %in% ReMeasure_C4DIF_name[1:10] == TRUE)
C4_RNAseq_nonre = which(!(colnames(C4DIF_RNAseq) %in% ReMeasure_C4DIF_name) == TRUE)[1:7]
C4_Agilent = C4DIF_Agilent[ , sort( union(C4_Agilent_re, C4_Agilent_nonre) )]
C4_RNAseq = C4DIF_RNAseq[ , sort( union(C4_RNAseq_re, C4_RNAseq_nonre))]

C5_Agilent_re = which(colnames(C5PRO_Agilent) %in% ReMeasure_C5PRO_name[1:10] == TRUE)
C5_Agilent_nonre = which( !(colnames(C5PRO_Agilent) %in% ReMeasure_C5PRO_name[1:10]) == TRUE)[1:59]
C5_RNAseq_re = which(colnames(C5PRO_RNAseq) %in% ReMeasure_C5PRO_name[1:10] == TRUE)
C5_RNAseq_nonre = which(!(colnames(C5PRO_RNAseq) %in% ReMeasure_C5PRO_name) == TRUE)[1:7]
C5_Agilent = C5PRO_Agilent[ , sort( union(C5_Agilent_re, C5_Agilent_nonre) )]
C5_RNAseq = C5PRO_RNAseq[ , sort( union(C5_RNAseq_re, C5_RNAseq_nonre)) ]

ds1 = cbind(C1_RNAseq, C2_RNAseq, C4_RNAseq, C5_RNAseq)
ds2 = cbind(C1_Agilent, C2_Agilent, C4_Agilent, C5_Agilent)
Ind_all = c( ReMeasure_C1MES_name[1:10], ReMeasure_C2IMM_name[1:10], 
          ReMeasure_C4DIF_name[1:10], ReMeasure_C5PRO_name[1:10] )

n1s = seq(10, length(Ind_all), by = 5)
library(doSNOW)
library(foreach)
nCPUS = 20
cl <- makeCluster(nCPUS, type = "SOCK", outfile = 'sva_ComBat_OvarianNULL.txt')
registerDoSNOW(cl)
 res <- foreach( n1=iter(n1s), .combine = combine) %do% {
   Ind = Ind_all[1:n1]
   ds3 = Agilent.dat[ ,Ind]
   dat = as.matrix( t( rbind(t(ds1), t(ds2), t(ds3)) ) )
   batch = c(rep(1, ncol(ds1)), rep(2, ncol(ds2)), rep(2, ncol(ds3 )) )
   batch = as.factor(batch)
   case_control = c(rep(0, ncol(ds1)), rep(1, ncol(ds2)), rep(0, ncol(ds3)) )
   treatment_frame = data.frame(case = case_control)
   mod = model.matrix(~as.factor(case), data = treatment_frame)
   mod0 = model.matrix(~1, data = treatment_frame)
   # We create model that only preserve remeasured data in batch 2
   ds1_ind = ds1[ , which( !(colnames(ds1) %in% Ind)), drop = F]
   dat_only = as.matrix( t( rbind( t(ds1_ind), t(ds2), t(ds3)) ) )
   batch_only = as.factor( c(rep(1, ncol(ds1_ind)), rep(2, ncol(ds2)), rep(2, ncol(ds3 ))) )
   cc_only = c( rep(0, ncol(ds1_ind)), rep(1, ncol(ds2)), rep(0, ncol(ds3)) )
   frame_only = data.frame(case = cc_only)
   mod_only = model.matrix(~as.factor(case), data = frame_only)
   mod0_only = model.matrix(~1, data = frame_only)
   # ComBat
   dat.par = ComBat(dat = dat_only, batch = batch_only, mod = mod_only, par.prior = TRUE, prior.plots = FALSE)
   pValuesComBat_par = f.pvalue(dat.par, mod_only, mod0_only)
   pValuesComBat_par = as.vector(pValuesComBat_par)
   # sva
   n.sv = num.sv(dat_only, mod_only, method = 'be')
   svobj_ovariance = sva(dat_only, mod_only, n.sv = n.sv)
  # combine with SV
   modSv = cbind(mod_only, svobj_ovariance$sv)
   mod0Sv = cbind(mod0_only, svobj_ovariance$sv)
   pValuesSv = f.pvalue(dat_only, modSv, mod0Sv)
   pValuesSv = as.vector(pValuesSv)
   # Our remeasure method
   Yc1mat = ds1
   Yt2mat = ds2
   Yc2mat = ds3
   Zc1 = matrix(1, nrow = ncol(Yc1mat), 1)
   Zt2 = matrix(1, nrow = ncol(Yt2mat), 1)
   Zc2 = matrix(1, nrow = ncol(Yc2mat), 1)
   output <- foreach (i=1:nrow(Yc1mat), .combine = 'c') %dopar% {
     library(RConics)
     Yc1 = Yc1mat[i, ]
     Yt2 = Yt2mat[i, ]
     Yc2 = Yc2mat[i, ]
     nc2 = length(Yc2)
     Index = which(names(Yc1) %in% Ind)
     Estimate = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7)
     rhoH = Estimate$rho
     sigma1H = Estimate$sigma1
     objVec = Estimate$objVec
     a0H = Estimate$a0
     a0Var = Estimate$a0Var
     a1H = Estimate$a1
     betaH = Estimate$beta
     sigma2H = Estimate$sigma2
     pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var)))
     return(pv)
   }
   return(list("pv" = output, "pv_par" = pValuesComBat_par,
                "pv_sva" = pValuesSv
              ))
 }

 stopCluster(cl)
 savepath = paste0("./Data_real/IndependSvaComBatNULL-", paste(Datalist, collapse = '-'), ".RData")
 save(res, file = savepath)






