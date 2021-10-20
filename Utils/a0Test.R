### Test statistics ###
library(reshape2)
library(ggplot2)
collectTest = function(allMethod, levelSeq, nc1, nc2, s1, s2, r, a0, a1, repID = 2){
  TestData = data.frame()
    nt2 = nc1
    #nt2 = nc1/2
    for (m in 1:length(allMethod)) {
      method = allMethod[m]
      fileList = paste0("./simuData/S1-method", method, "-repID-", 
                        repID, "-nc1-", nc1, "-nt2-", nt2, "-nc2-", nc2, "-s1-", s1, "-s2-",
                        s2, "-rho-", r, "-a0-", a0, "-a1-", a1, ".RData")
      load(fileList)
      a0List = Final$a0
      a0VarList = Final$a0Var
      ## We need add an absolute value here 
      ztest = unlist(a0List)/sqrt(unlist(a0VarList) )
      
      tmp = data.frame()
      for (l in 1:length(levelSeq) ){
        level = levelSeq[l]
      # We need to use the quantile of standard normal distribution 
        p = qnorm(1-level/2, mean = 0, sd = 1)
        #names(t) = NULL
        rej = mean(abs(ztest) > p)
        ymax = rej + 1.96*sqrt(rej*(1 - rej)/length(ztest) ) 
        ymin = rej - 1.96*sqrt(rej*(1 - rej)/length(ztest) ) 
        ymin <- ifelse(ymin > 0, ymin, 0)
        tmp = rbind(tmp, data.frame(method = method, rate = rej, level = level, a0 = a0, ymax = ymax, ymin = ymin) )
      }
      TestData = rbind(TestData, tmp)
      TestData$level = factor(TestData$level)
    }
  return(TestData)
}

collectTest_all = function(allMethod, nc1Seq, levelSeq, s1, s2, r, a0, a1, repID ) {
  TestData = data.frame()
  for (i in 1:length(nc1Seq)){
    nc1 = nc1Seq[i]
    nc2Seq = c(nc1/4, nc1/2, nc1*3/4)
   # nc2Seq = c(nc1/2, nc1/5, nc1/10)
    nc2String = c("nc1/4", "nc1/2", "nc1*3/4")
   # nc2String = c("nc1/2", "nc1/5", "nc1/10")
    tmp1 = data.frame()
    for (j in 1:length(nc2Seq)) {
      nc2 = nc2Seq[j]
      tmpData = collectTest(allMethod, levelSeq, nc1, nc2, s1, s2, r, a0, a1, repID )
      tmpData$nc1 = nc1
      tmpData$nc2 = nc2String[j]
      tmp1 = rbind(tmp1, tmpData)
    }
    TestData = rbind(TestData, tmp1)
  }
  return(TestData)
}


#allMethod = c("Ignore", "Batch2", "ReMeasure")
#allMethod = c( "ReMeasure")
allMethod = c("Ignore", "ReMeasure")
levelSeq = c(0.005, 0.025, 0.05, 0.1)
#nc1Seq = c(100, 200, 500)
nc1Seq = c(100, 200, 500)
s1 = 1
s2 = 2
r = 0.5
a0 = 0 
a1 = 0.5 
repID = 2
#nc2 = 25

TestData_0 = collectTest_all(allMethod, nc1Seq, levelSeq, s1, s2, r, a0 = 0, a1, repID)
TestData_0.5 = collectTest_all(allMethod, nc1Seq, levelSeq, s1, s2, r, a0 = 0.5, a1, repID)
TestData_1 = collectTest_all(allMethod, nc1Seq, levelSeq, s1, s2, r, a0 = 1, a1, repID)
TestData = rbind(TestData_0, TestData_0.5, TestData_1)

TestData$nc1 = factor(TestData$nc1)
#pdf("./figs/Test-a0-2.pdf", width = 7, height = 6.18)
#dev.off()
nMeth = length(allMethod)
fill_cols <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7",  "#F0E442",   
          "steelblue", "orange", 'darkseagreen', 'firebrick2', 'gold2')[1:nMeth]  # Colorblind friendly
cols <- c("#56B4E9", "#009E73")
hline_data = data.frame(level = c(0.005, 0.025, 0.05, 0.1),
                        value = c(0.005, 0.025, 0.05, 0.1))
dodge <- position_dodge(width=0.92)
plot = ggplot(data = TestData, aes(x = nc2, y = rate,  fill = method, color = nc1) ) + 
  geom_bar(stat='identity',position=dodge) + 
  geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25, size = 0.25) +
  geom_hline(data = hline_data, aes(yintercept = value ), linetype = 2) +
  facet_grid(level~a0, scale = 'free', labeller = label_both)  + guides(fill=guide_legend(title="Method")) +
  labs(title = 'n2 = n1/2, rho = 0.5', x = "Re-measure Size", y ="Power" ) + 
  theme(plot.title = element_text(hjust = 0.5)) 
plot 
#+ scale_fill_manual(values = fill_cols) + scale_color_manual(values = cols)


