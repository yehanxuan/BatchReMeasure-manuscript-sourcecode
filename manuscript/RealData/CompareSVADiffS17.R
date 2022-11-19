rm(list = ls())
library(sva)
library(pamr)
library(limma)
library(stringi)
source("./code/OvarianCancer.R")
source("./code/ReMeasure_S1.R")
DList = list(
  c("C1C2C5_RNAseq", "C4DIF_Agilent", "C1C2C5_Agilent"),
  c("C1C4C5_RNAseq", "C2IMM_Agilent", "C1C4C5_Agilent"),
  c("C1C2C4_RNAseq", "C5PRO_Agilent", "C1C2C4_Agilent"),
  c("C2C4C5_RNAseq", "C1MES_Agilent", "C2C4C5_Agilent")
)

combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}

folder = "Data_real"
if (file.exists(folder)) {
  cat("The folder already exists") 
} else {
  dir.create(folder)
}

for (l in 1:length(DList)) {
  Datalist = DList[[l]]
  ds1 = get(Datalist[1])
  ds2 = get(Datalist[2])
  if ( isTRUE(grepl("C1MES", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C1MES_name
  } else if ( isTRUE(grepl("C2IMM", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C2IMM_name
  } else if ( isTRUE(grepl("C4DIF", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C4DIF_name
  } else if ( isTRUE(grepl("C5PRO", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C5PRO_name
  } else if ( isTRUE(grepl("C1C2C4", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C1C2C4_name
  } else if ( isTRUE(grepl("C1C2C5", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C1C2C5_name
  } else if ( isTRUE(grepl("C1C4C5", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C1C4C5_name
  } else if (isTRUE(grepl("C2C4C5", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C2C4C5_name
  }
  n1s = seq(10, length(Ind_all), by = 5)
  library(doSNOW)
  library(foreach)
  nCPUS = 10
  cl <- makeCluster(nCPUS)
  registerDoSNOW(cl)
  
  res <- foreach( n1=iter(n1s), .combine = combine) %do% {
    Ind = Ind_all[1:n1]
    ds3 = get(Datalist[3])[ ,Ind]  
    dat <- rbind(t(ds1), t(ds2),t((ds3)))
    dat <- as.matrix(t(dat))
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
    output <-foreach (i=1:nrow(Yc1mat), .combine = 'c') %dopar% {
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
                 "pv_sva" = pValuesSv))
  }
  stopCluster(cl)
  savepath = paste0("./Data_real/IndependSvaComBat-", paste(Datalist, collapse = '-'), ".RData")
  save(res, file = savepath)
}


