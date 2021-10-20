### This section is used for Hypothesis testing ###
# For example 

Final$a0

Final$a0Var

a0List = Final$a0
a0VarList = Final$a0Var

Scaled_Statistics = function(a0List, a0VarList) {
  a0Scale = unlist(a0List)
  a0VarVec = unlist(a0VarList)
  Test_Stat = a0Scale/sqrt(a0VarVec)
  return(Test_Stat)
}

### 
quantile(a0Normal, 0.975)
quantile(a0Normal, 0.95)
quantile(a0Normal, 0.9)




#quantile(a0Normal, c(0.025, 0.975) )
# --- Create --- Quantile table ---#
allMethod
nc1 = 500
nt2 = nc1 
nc2Seq = c(nc1/4, nc1/2, nc1*3/4) 
s1 = 1
s2 = 2
r = 0.5
a0 = 2
CollectTestStatistics = function(){
  FinalQuantile = c()
  for (c2 in 1:length(nc2Seq)){
    tmpQuantile  = matrix(0, 2, length(allMethod))
    
    nc2 = nc2Seq[c2]
    for (m in 1:length(allMethod)) {
      method = allMethod[m]
      fileList = paste0("./simuData/S1-",method,"-",repID, "-", nc1, "-", nt2, "-",
                        nc2, "-", s1, "-", s2, "-", r, ".RData")
      load(fileList)
      a0List = Final$a0
      a0VarList = Final$a0Var
      a0Normal = Scaled_Statistics(a0List, a0VarList, a0)
      Q = quantile(a0Normal, c(0.025, 0.975))
      tmpQuantile[1, m] = Q[1]
      tmpQuantile[2, m] = Q[2]
    }
    FinalQuantile  = cbind(FinalQuantile, tmpQuantile)  
  }
  return(FinalQuantile)
}

QuantileTable = function() {
  rk = 3
  cat("\\hline\n")
  
  for (rowI in 1:nrow(FinalQuantile)) {
    if (rowI == 1) {
      cat("$q_{0.025}$")
    } else {
      cat("$q_{0.975}$")
    }
    for (colJ in 1:9) {
      cat(" & ")
        cat(format(round( FinalQuantile[rowI, colJ],rk), nsmall = rk))
    }
    cat("\\\\", "\n")
  }
  cat("\\hline\n")
}

#----Create Barplot----####
allMethod = c("Ignore", "Batch2", "ReMeasure")
levelSeq = c(0.025, 0.05, 0.1)
nc1 = 500
nt2 = 500
nc2 = 125
repID = 2
s1 = 1
s2 = 2
r = 0.5

TestPlot = function(allMethod, levelSeq, nc1, nt2, nc2, s1, s2, r, repID = 2) {
  BarData = data.frame()
  for (m in 1:length(allMethod)) {
    method = allMethod[m]
    fileList = paste0("./simuData/S1-",method,"-",repID, "-", nc1, "-", nt2, "-",
                      nc2, "-", s1, "-", s2, "-", r, ".RData")
    load(fileList)
    a0List = Final$a0
    a0VarList = Final$a0Var
    Test_Stat = Scaled_Statistics(a0List, a0VarList)
    tmp = data.frame()
    for (l in 1:length(levelSeq) ){
      level = levelSeq[l]
      # We need to use the quantile of standard normal distribution 
      p = qnorm(1-level, mean = 0, sd = 1)
      names(t) = NULL
     Ratio = sum(Test_Stat > p)/length(Test_Stat)
      #count = as.numeric(sum(Test_Stat > p))
      tmp = rbind(tmp, data.frame(method = method, Frequency = Ratio, level = level) )
    }
    BarData = rbind(BarData, tmp)
  }
}


BarData$method
BarData$level = factor(BarData$level)
library(ggplot2)
#aes( x = method, y=t, fill = level)
ggplot(data = BarData) + geom_boxplot(mapping = aes( x = method, y = Frequency, fill = level))

ggplot(data = BarData) + geom_bar(mapping = aes( x = count, fill = level), position = "dodge")

ggplot(data = BarData) + stat_summary(mapping = aes(x = method, y = Frequency, color = level), fun.min = 
                                        min, fun.max = max, fun = median)



