rm(list = ls())
source("./code/OvarianCancer.R")
source("./code/ReMeasure_S1.R")
source("./code/Batch2Only.R")
source("./code/IgnoreBatch.R")
source("./code/LSmodel.R")

DList = list(
  c("C1C2C5_RNAseq", "C4DIF_Agilent", "C1C2C5_Agilent"),
  c("C1C4C5_RNAseq", "C2IMM_Agilent", "C1C4C5_Agilent"),
  c("C1C2C4_RNAseq", "C5PRO_Agilent", "C1C2C4_Agilent"),
  c("C2C4C5_RNAseq", "C1MES_Agilent", "C2C4C5_Agilent")
)

allMethods = c("ReMeasure", "Batch2", "Ignore", "LS")
# Two sample t-test
T_pvVec = NULL
for (i in 1:nrow(C1MES_Agilent)) {
  x = C1MES_Agilent[i, ]
  y = C2C4C5_Agilent[i, ]
  if ( (0.95 < var(x)/var(y) ) && ( var(x)/var(y) < 1.05 ) ) {
    Equal = TRUE
  } else {
    Equal = FALSE
  }
  Ttest = t.test(x, y, mu = 0, alternative = "two.sided", var.equal = Equal, paired = F)
  T_pvVec = c(T_pvVec, Ttest$p.value)
}

alpha_t = 0.2
Rej_Bon_t_Ind = (p.adjust(T_pvVec, method = "bonferroni") < alpha_t)
sum(Rej_Bon_t_Ind)
Sig_Ind = which(Rej_Bon_t_Ind == TRUE)

folder = "Data_real"
if (file.exists(folder)) {
  cat("The folder already exists") 
} else {
  dir.create(folder)
}

strList = c("Strong")
#strList = c("Weak", "Strong")
for (strength in strList) {
  for (l in 1:length(DList)) {
    
    Datalist = DList[[l]]
    Yc1mat = get(Datalist[1])
    Yt2mat = get(Datalist[2])
    Gene_Index = get(paste0("Gene_", strength))
    Yc1mat = Yc1mat[Gene_Index, ]
    Yt2mat = Yt2mat[Gene_Index, ]
  
    if (isTRUE( grepl("recur", Datalist[1], fix=TRUE) ) ) {
      Ind_all = ReMeasure_recur_name
    } else if ( isTRUE(grepl("non",Datalist[1], fix = TRUE ) ) ) {
      Ind_all = ReMeasure_non_name
    } else if ( isTRUE(grepl("C1MES", Datalist[1], fix = TRUE)) ) {
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
    
    n1s = seq(5, length(Ind_all), by = 5)
    library(doSNOW)
    nCPUS = 10 
    cl <- makeCluster(nCPUS, type = "SOCK", outfile = 'OvarianCancer_log.txt')
    registerDoSNOW(cl)
    library(foreach)

    for (m in 1:length(allMethods) ) {
      method = allMethods[m]
      res <- foreach (n1=iter(n1s), .combine = 'cbind') %do% {
        Ind = Ind_all[1:n1]
        Yc2mat = get(Datalist[3])[ ,Ind]
        # with strength 
        Yc2mat = Yc2mat[Gene_Index, ]
        Zc1 = matrix(1, nrow = ncol(Yc1mat), 1)
        Zt2 = matrix(1, nrow = ncol(Yt2mat), 1)
        Zc2 = matrix(1, nrow = ncol(Yc2mat), 1)
        output <- foreach (i = 1:nrow(Yc1mat), .combine = 'c') %dopar% {
          library(RConics)  
          Yc1 = Yc1mat[i, ]
          Yt2 = Yt2mat[i, ]
          Yc2 = Yc2mat[i, ]
          length(Yc1); length(Yt2);length(Yc2)
          nc2 = length(Yc2)
          length(Ind)
          Index = which(names(Yc1) %in% Ind)
          if (method == "ReMeasure") {
            Estimate = Estimate_ReMeasure_S1(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, Index, tol.c = 1e-7)
            rhoH = Estimate$rho
            sigma1H = Estimate$sigma1
            objVec = Estimate$objVec
          } else if (method == "Batch2") {
            Estimate = Estimate_OnlyBatch2(Zt2, Zc2, Yt2, Yc2)
            rhoH = NULL
            sigma1H = NULL
            objVec = NULL
          } else if (method == "Ignore") {
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
          a1H = Estimate$a1
          betaH = Estimate$beta
          sigma2H = Estimate$sigma2
          pv <- 2 * stats::pnorm(-abs(a0H / sqrt(a0Var))) 
          Time = Estimate$Time
          return(pv)
        }
        return(output)
      }
      colnames(res) <- n1s
      if (is.null(strength)) {
        write.csv(res, file = paste0( './Data_real/', paste(Datalist, collapse = '-'), "-method-", method,
                                      ".csv"), row.names = FALSE)
      } else {
        write.csv(res, file = paste0( './Data_real/strength-', strength, "-", paste(Datalist, collapse = '-'), "-method-", method,
                                      ".csv"), row.names = FALSE)
      }
    }
    stopCluster(cl)
  }
}
  


