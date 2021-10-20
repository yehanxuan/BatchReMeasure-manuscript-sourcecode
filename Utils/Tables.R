# Compute sample mean and sample sd #
compute = function(resultT){
  res = Reduce("+", resultT)/length(resultT)
  resultS2 = lapply(resultT, function(tmp) tmp^2)
  res2 = Reduce("+", resultS2)/length(resultT)
  resSD = sqrt(res2 - res^2)/sqrt(length(resultT))
  return(list(res = res, resSD = resSD))
}

# computeMSE 
#resultBetaList = Final$beta 
#betaTrue = beta
MSEbeta = function(resultBetaList, betaTrue) {
  err = lapply(resultBetaList, function(x) (x - betaTrue)^2)
  Norm = lapply(err, function(x) sum(x)) # \|\beta - \betahat\|_2^2
  #lapply(err, col)
  res = Reduce("+", Norm)/length(Norm)
  resultS2 = lapply(Norm, function(tmp) tmp^2)
  res2 = Reduce("+", resultS2)/length(Norm)
  resSD = sqrt(res2 - res^2)/sqrt(length(Norm))
  return(list(res = res, resSD = resSD))
}

computeMSE = function(resultT, trueValue) {
  err = lapply(resultT, function(x) (x - trueValue)^2 )
  res = Reduce("+", err)/length(err)
  resultS2 = lapply(err, function(tmp) tmp^2)
  res2 = Reduce("+", resultS2)/length(err)
  resSD = sqrt(res2 - res^2)/sqrt(length(err))
  return(list(res = res, resSD = resSD))
}



#lapply(Final$beta, function)

#------ S1 Tables ----

allMethod = c("Ignore", "Batch2", "ReMeasure")
nc1 = 200
nt2 = nc1/2 
#nt2 = nc1 
nc2Seq = c(nc1/10, nc1/5, nc1/2) 

#nc2Seq = c(nc1/4, nc1/2, nc1*3/4)
beta = c(1, 1, 1, 1, 1)
s1 = 1
s2 = 2
r = 0.8
repID=2 
a0 = 1
a1 = 0.5
collectEstimates = function(allMethod, nc1, nt2, nc2Seq, s1, s2, r, a0, a1, repID, beta) {
  #FinalMean = matrix(0, 5, length(allMethod)*length(nc2Seq))
  #FinalSD = matrix(0, 5, length(allMethod)*length(nc2Seq))
  FinalMean = c()
  FinalSD = c()
  ValueList = list(a0, r, s1, s2, a1)
 for (c2 in 1:length(nc2Seq)){
   tmpMean = matrix(0, 6, length(allMethod))
   tmpSD = matrix(0, 6, length(allMethod))
   nc2 = nc2Seq[c2]
   for (m in 1:length(allMethod)) {
     method = allMethod[m]
   #  fileList = paste0("./simuDataSmall/S1-",method,"-",repID, "-", nc1, "-", nt2, "-",
    #                   nc2, "-", s1, "-", s2, "-", r, ".RData")
     fileList = paste0("./simuData/S1-method", method, "-repID-", 
     repID, "-nc1-", nc1, "-nt2-", nt2, "-nc2-", nc2, "-s1-", s1, "-s2-",
     s2, "-rho-", r, "-a0-", a0, "-a1-", a1, ".RData")
     load(fileList)
     k = 0
     for (i in c(1,6,4,5,3, 7)) {
       k = k + 1
       if (i == 7) {
        resultT = MSEbeta(Final$beta, beta) 
       } else {
         #resultT = compute(Final[[i]])
         resultT = computeMSE(Final[[i]], ValueList[[k]])
       }
       if (identical(resultT$res, numeric(0)) ) {
         tmpMean[k, m] = 999
         tmpSD[k, m] = 999
       } else {
         tmpMean[k, m] = resultT$res
         tmpSD[k, m] = resultT$resSD
       }
     }
    # compute(Final$beta)
     
   }
   FinalMean = cbind(FinalMean, tmpMean)
   FinalSD = cbind(FinalSD, tmpSD)
 }
  return(list(FinalMean, FinalSD))
}

FinalMean = collectEstimates(allMethod, nc1, nt2, nc2Seq, s1, s2, r, a0, a1, repID = 2, beta)[[1]]
FinalSD = collectEstimates(allMethod, nc1, nt2, nc2Seq, s1, s2, r, a0, a1, repID = 2, beta)[[2]]

#Names = names(Final)[c(1,6,4,5,3)]
Names = c("$\\|\\hat{a}_0 - a_0 \\|_2^2 $", "$\\| \\hat{\\rho} -\\rho \\|_2^2  $", "$\\| \\hat{\\sigma}_1 - \\sigma_1 \\|_2^2$", 
          "$\\| \\hat{\\sigma}_2 - \\sigma_2 \\|_2^2$", "$\\| \\hat{a}_1 - a_1 \\|_2^2$", "$\\|\\hat{\\beta} - \\beta_0 \\|_2^2$")
Output_Table = function(FinalMean, FinalSD, Names) {
  rk = 3
  cat("\\hline\n")
  
  for (rowI in 1:nrow(FinalMean)) {
    cat(Names[rowI])
    for (colJ in 1:9) {
      cat(" & ")
      if (FinalMean[rowI, colJ] == 999){
        cat("-")
      } else {
        cat(format(round( FinalMean[rowI, colJ],rk), nsmall = rk))
        cat("(",format(round(FinalSD[rowI, colJ],rk),nsmall = rk),")", sep="" )
      }
    }
    cat("\\\\", "\n")
  }
  cat("\\hline\n")
}

Output_Table(FinalMean, FinalSD, Names)



