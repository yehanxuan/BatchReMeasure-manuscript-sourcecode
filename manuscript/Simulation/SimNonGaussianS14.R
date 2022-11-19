## Main_NonGaussian ##
rm(list = ls())
library(snowfall)
library(RConics)
source("./code/ReMeasure_S1.R")
source("./code/Batch2Only.R")

args <- commandArgs(trailingOnly = TRUE)
method = args[1]
n = as.numeric(args[2])
n1 = as.numeric(args[3])
r1 = as.numeric(args[4])
r2 = as.numeric(args[5])
a0 = as.numeric(args[6])
a1 = as.numeric(args[7])
repID = as.numeric(args[8])
NoiseType = args[9]

method = "Batch2"
NoiseType = "t"
n = 100
n1 = 50
r1 = 2
r2 = 0.6
a0 = 0.25
a1 = 0.5
repID = 2



source(paste0("./Simulation/oneReplicate-NonGaussian-S1.R"))

savePath = paste0("./S1Data_NG/S1-method-", method, "-NoiseType-", NoiseType, "-repID-",repID,"-nc1-", n/2 ,"-nt2-", n/2,"-nc2-", n1 ,"-r1-",
                  r1,"-r2-", r2,"-a0-", a0, "-a1-", a1, ".RData")


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

if (method == "ReMeasure") {
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
    tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_ReMeasure_NG)
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
    result4 = c(result4, tmp4)
    result5 = c(result5, tmp5)
    result6 = c(result6, tmp6)
    result7 = c(result7, tmp7)
    result8 = c(result8, tmp8)
    result9 = c(result9, tmp9)
    a0Var_ora = c(a0Var_ora, tmp10)
    Power_th = c(Power_th, tmp11)
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
    tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_OnlyBatch2_NG)
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
    result1 = c(result1, tmp1)
    result2 = c(result2, tmp2)
    result3 = c(result3, tmp3)
    result4 = c(result4, tmp4)
    result5 = c(result5, tmp5)
    result6 = c(result6, tmp6)
    result7 = c(result7, tmp7)
    result8 = c(result8, tmp8)
  }
}

Final = list("a0" = result1, "a0Var"= result2, "a1"=result3, "sigma1"=result4,
             "sigma2" = result5, "rho" = result6, "beta"=result7, "objVec"=result8,
             "Time" = result9, "ztest" = result10, "ztestb" = result11, "p.value" = result12, 
             "a0Var_ora" = a0Var_ora, "Power" = Power_th)


save(Final, file=savePath)


