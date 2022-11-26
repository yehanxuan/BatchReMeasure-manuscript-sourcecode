rm(list = ls())
library(snowfall)
library(RConics)
library(RcppArmadillo)
library(Rcpp)
source("./code/ReMeasure_S1.R")
source("./code/Batch2Only.R")
source("./code/IgnoreBatch.R")
source("./code/Bootstrap.R")
source("./code/LSmodel.R")
source("./code/CorFunc.R")
Rcpp::sourceCpp("code/RcppReMeasure_S1.cpp")
source("./code/RcppReMeasure_S1.R")
  
args <- commandArgs(trailingOnly = TRUE)
method = args[1]
n = as.numeric(args[2])
n1 = as.numeric(args[3])
r1 = as.numeric(args[4])
r2 = as.numeric(args[5])
a0 = as.numeric(args[6])
a1 = as.numeric(args[7])
repID = as.numeric(args[8])

 method = "ReMeasure"
 n = 100
 n1 = 40
 r1 = 1
 r2 = 0.6
 a0 = 0.25
 a1 = 0.5
 repID = 2

source(paste0("./Simulation/oneReplicate-New-S1.R"))

# specify the save path
if (method != "RcppReMeasure") {
  savePath = paste0("./S1Data_Revise/S1-method-", method,"-repID-",repID,"-nc1-", n/2 ,"-nt2-", n/2,"-nc2-", n1 ,"-r1-",
                    r1,"-r2-", r2,"-a0-", a0, "-a1-", a1, ".RData")
} else {
  savePath = paste0("./S1Data_Rcpp_Revise/S1-method-", method,"-repID-",repID,"-nc1-", n/2 ,"-nt2-", n/2,"-nc2-", n1 ,"-r1-",
                    r1,"-r2-", r2,"-a0-", a0, "-a1-", a1, ".RData")
}
 



nCPUS = 10
maxIter = 10

result1 = list()
result2 = list()
result3 = list()
result4 = list()
result5 = list()
result6 = list()
result7 = list()
result8 = list()
result9 = list()
result10 = list()
result11 = list()
result12 = list()
a0Var_ora = list()
Power_th = list()
if ((method == "Ignore") | (method == "Gen") ) {
  for (i in 1:ceiling(maxIter/nCPUS)){
    print(i)
    sfInit(parallel = TRUE, cpus = nCPUS)
    sfExportAll()
    sfLibrary(RConics)
    sfLibrary(MASS)
    sfLibrary(mvtnorm)
    sfLibrary(lmvar)
    sBegin = (i-1)*nCPUS +1
    sEnd = min(i*nCPUS, maxIter)
    seedSeq = seq(sBegin, sEnd, by = 1)
    if (method == "Ignore") {
      tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_Ignore)
    } else if (method == "Gen") {
      tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_Gen)
    }
    sfStop()
    tmp1 = lapply(tmp, function(x) x[[1]])
    tmp2 = lapply(tmp, function(x) x[[2]])
    tmp3 = lapply(tmp, function(x) x[[3]])
    tmp4 = lapply(tmp, function(x) x[[4]])
    tmp5 = lapply(tmp, function(x) x[[5]])
    tmp6 = lapply(tmp, function(x) x[[6]])
    tmp7 = lapply(tmp, function(x) x[[7]])
    tmp8 = lapply(tmp, function(x) x[[8]])
    tmp9 = lapply(tmp, function(x) x[[9]])
    result1 = c(result1, tmp1)
    result2 = c(result2, tmp2)
    result3 = c(result3, tmp3)
    result4 = c(result4, tmp4)
    result5 = c(result5, tmp5)
    result6 = c(result6, tmp6)
    result7 = c(result7, tmp7)
    result8 = c(result8, tmp8)
    result9 = c(result9, tmp9)
  }
} else if (method == "Batch2") {
  for (i in 1:ceiling(maxIter/nCPUS)){
    print(i)
    sfInit(parallel = TRUE, cpus = nCPUS)
    sfExportAll()
    sfLibrary(RConics)
    sfLibrary(MASS)
    sfLibrary(mvtnorm)
    sBegin = (i-1)*nCPUS +1
    sEnd = min(i*nCPUS, maxIter)
    seedSeq = seq(sBegin, sEnd, by = 1)
    tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_OnlyBatch2)
    sfStop()
    tmp = Filter(function(x) length(x) >= 3, tmp)
    tmp1 = lapply(tmp, function(x) x[[1]])
    tmp2 = lapply(tmp, function(x) x[[2]])
    tmp3 = lapply(tmp, function(x) x[[3]])
    tmp4 = lapply(tmp, function(x) x[[4]])
    tmp5 = lapply(tmp, function(x) x[[5]])
    tmp6 = lapply(tmp, function(x) x[[6]])
    tmp7 = lapply(tmp, function(x) x[[7]])
    tmp8 = lapply(tmp, function(x) x[[8]])
    tmp9 = lapply(tmp, function(x) x[[9]])
    result1 = c(result1, tmp1)
    result2 = c(result2, tmp2)
    result3 = c(result3, tmp3)
    result4 = c(result4, tmp4)
    result5 = c(result5, tmp5)
    result6 = c(result6, tmp6)
    result7 = c(result7, tmp7)
    result8 = c(result8, tmp8)
    result9 = c(result9, tmp9)
  }
} else if ( (method == "ReMeasure") | (method == "RcppReMeasure")) {
  for (i in 1:ceiling(maxIter/nCPUS)){
    print(i)
    sfInit(parallel = TRUE, cpus = nCPUS)
    sfExportAll()
    sfLibrary(RConics)
    sfLibrary(MASS)
    sfLibrary(mvtnorm)
    sfLibrary(RcppArmadillo)
    sfLibrary(Rcpp)
    sfClusterEval(sourceCpp("code/RcppReMeasure_S1.cpp"))
    sBegin = (i-1)*nCPUS +1
    sEnd = min(i*nCPUS, maxIter)
    seedSeq = seq(sBegin, sEnd, by = 1)
    if (method == "ReMeasure") {
      tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_ReMeasure)
    } else if (method == "RcppReMeasure") {
      tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_ReMeasure_Rcpp)
    }
    sfStop()
    tmp = Filter(function(x) length(x) >= 3, tmp)
    tmp1 = lapply(tmp, function(x) x[[1]])
    tmp2 = lapply(tmp, function(x) x[[2]])
    tmp3 = lapply(tmp, function(x) x[[3]])
    tmp4 = lapply(tmp, function(x) x[[4]])
    tmp5 = lapply(tmp, function(x) x[[5]])
    tmp6 = lapply(tmp, function(x) x[[6]])
    tmp7 = lapply(tmp, function(x) x[[7]])
    tmp8 = lapply(tmp, function(x) x[[8]])
    tmp9 = lapply(tmp, function(x) x[[9]])
    tmp10 = lapply(tmp, function(x) x[[10]])
    tmp11 = lapply(tmp, function(x) x[[11]])
    result1 = c(result1, tmp1)
    result2 = c(result2, tmp2)
    result3 = c(result3, tmp3)
    result4 = c(result4, tmp4)
    result5 = c(result5, tmp5)
    result6 = c(result6, tmp6)
    result7 = c(result7, tmp7)
    result8 = c(result8, tmp8)
    result9 = c(result9, tmp9)
    a0Var_ora = c(a0Var_ora, tmp10)
    Power_th = c(Power_th, tmp11)
  }
} else if (method == "ResBoot") {
  for (i in 1:ceiling(maxIter/nCPUS)){
  print(i)
  sfInit(parallel = TRUE, cpus = nCPUS)
  sfExportAll()
  sfLibrary(RConics)
  sfLibrary(MASS)
  sfLibrary(mvtnorm)
  sBegin = (i-1)*nCPUS +1
  sEnd = min(i*nCPUS, maxIter)
  seedSeq = seq(sBegin, sEnd, by = 1)
  tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_Res)
  sfStop()
  tmp = Filter(function(x) length(x) >= 3, tmp)
  tmp1 = lapply(tmp, function(x) x[[1]])
  tmp2 = lapply(tmp, function(x) x[[2]])
  tmp3 = lapply(tmp, function(x) x[[3]])
  tmp4 = lapply(tmp, function(x) x[[4]])
  tmp5 = lapply(tmp, function(x) x[[5]])
  tmp6 = lapply(tmp, function(x) x[[6]])
  tmp7 = lapply(tmp, function(x) x[[7]])
  tmp8 = lapply(tmp, function(x) x[[8]])
  tmp9 = lapply(tmp, function(x) x[[9]])
  tmp10 = lapply(tmp, function(x) x[[10]])
  tmp11 = lapply(tmp, function(x) x[[11]])
  result1 = c(result1, tmp1)
  result2 = c(result2, tmp2)
  result3 = c(result3, tmp3)
  result4 = c(result4, tmp4)
  result5 = c(result5, tmp5)
  result6 = c(result6, tmp6)
  result7 = c(result7, tmp7)
  result8 = c(result8, tmp8)
  result9 = c(result9, tmp9)
  result10 = c(result10, tmp10)
  result11 = c(result11, tmp11)
  }
} else if (method == "WildBoot") {
  for (i in 1:ceiling(maxIter/nCPUS)){
  print(i)
  sfInit(parallel = TRUE, cpus = nCPUS)
  sfExportAll()
  sfLibrary(RConics)
  sfLibrary(MASS)
  sfLibrary(mvtnorm)
  sBegin = (i-1)*nCPUS +1
  sEnd = min(i*nCPUS, maxIter)
  seedSeq = seq(sBegin, sEnd, by = 1)
  tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_Wild)
  sfStop()
  tmp = Filter(function(x) length(x) >= 3, tmp)
  tmp1 = lapply(tmp, function(x) x[[1]])
  tmp2 = lapply(tmp, function(x) x[[2]])
  tmp3 = lapply(tmp, function(x) x[[3]])
  tmp4 = lapply(tmp, function(x) x[[4]])
  tmp5 = lapply(tmp, function(x) x[[5]])
  tmp6 = lapply(tmp, function(x) x[[6]])
  tmp7 = lapply(tmp, function(x) x[[7]])
  tmp8 = lapply(tmp, function(x) x[[8]])
  tmp9 = lapply(tmp, function(x) x[[9]])
  tmp10 = lapply(tmp, function(x) x[[10]])
  tmp11 = lapply(tmp, function(x) x[[11]])
  result1 = c(result1, tmp1)
  result2 = c(result2, tmp2)
  result3 = c(result3, tmp3)
  result4 = c(result4, tmp4)
  result5 = c(result5, tmp5)
  result6 = c(result6, tmp6)
  result7 = c(result7, tmp7)
  result8 = c(result8, tmp8)
  result9 = c(result9, tmp9)
  result10 = c(result10, tmp10)
  result11 = c(result11, tmp11)
  }
} else if (method == "PairBoot") {
  for (i in 1:ceiling(maxIter/nCPUS)){
    print(i)
    sfInit(parallel = TRUE, cpus = nCPUS)
    sfExportAll()
    sfLibrary(RConics)
    sfLibrary(MASS)
    sfLibrary(mvtnorm)
    sBegin = (i-1)*nCPUS +1
    sEnd = min(i*nCPUS, maxIter)
    seedSeq = seq(sBegin, sEnd, by = 1)
    tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_Pair)
    sfStop()
    tmp = Filter(function(x) length(x) >= 3, tmp)
    tmp1 = lapply(tmp, function(x) x[[1]])
    tmp2 = lapply(tmp, function(x) x[[2]])
    tmp3 = lapply(tmp, function(x) x[[3]])
    tmp4 = lapply(tmp, function(x) x[[4]])
    tmp5 = lapply(tmp, function(x) x[[5]])
    tmp6 = lapply(tmp, function(x) x[[6]])
    tmp7 = lapply(tmp, function(x) x[[7]])
    tmp8 = lapply(tmp, function(x) x[[8]])
    tmp9 = lapply(tmp, function(x) x[[9]])
    tmp10 = lapply(tmp, function(x) x[[10]])
    tmp11 = lapply(tmp, function(x) x[[11]])
    result1 = c(result1, tmp1)
    result2 = c(result2, tmp2)
    result3 = c(result3, tmp3)
    result4 = c(result4, tmp4)
    result5 = c(result5, tmp5)
    result6 = c(result6, tmp6)
    result7 = c(result7, tmp7)
    result8 = c(result8, tmp8)
    result9 = c(result9, tmp9)
    result10 = c(result10, tmp10)
    result11 = c(result11, tmp11)
  }
} else if ( (method == "LS") ) {
  for (i in 1:ceiling(maxIter/nCPUS)){
    print(i)
    sfInit(parallel = TRUE, cpus = nCPUS)
    sfExportAll()
    sfLibrary(RConics)
    sfLibrary(MASS)
    sfLibrary(mvtnorm)
    sBegin = (i-1)*nCPUS +1
    sEnd = min(i*nCPUS, maxIter)
    seedSeq = seq(sBegin, sEnd, by = 1)
    tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_LS)
    sfStop()
    tmp = Filter(function(x) length(x) >= 3, tmp)
    tmp1 = lapply(tmp, function(x) x[[1]])
    tmp2 = lapply(tmp, function(x) x[[2]])
    tmp3 = lapply(tmp, function(x) x[[3]])
    tmp4 = lapply(tmp, function(x) x[[4]])
    tmp5 = lapply(tmp, function(x) x[[5]])
    tmp6 = lapply(tmp, function(x) x[[6]])
    tmp7 = lapply(tmp, function(x) x[[7]])
    tmp8 = lapply(tmp, function(x) x[[8]])
    tmp9 = lapply(tmp, function(x) x[[9]])
    tmp12 = lapply(tmp, function(x) x[[10]])

    result1 = c(result1, tmp1)
    result2 = c(result2, tmp2)
    result3 = c(result3, tmp3)
    result4 = c(result4, tmp4)
    result5 = c(result5, tmp5)
    result6 = c(result6, tmp6)
    result7 = c(result7, tmp7)
    result8 = c(result8, tmp8)
    result9 = c(result9, tmp9)
    result12 = c(result12, tmp12)
  }
}

Final = list("a0" = result1, "a0Var"= result2, "a1"=result3, "sigma1"=result4,
             "sigma2" = result5, "rho" = result6, "beta"=result7, "objVec"=result8,
             "Time" = result9, "ztest" = result10, "ztestb" = result11, "p.value" = result12, 
             "a0Var_ora" = a0Var_ora, "Power" = Power_th)

save(Final, file = savePath)


