### Density histograms ### 
rm(list = ls())
Hist_Plot = function(method, nc1Seq, a0 = 2, s1 = 1, s2 = 2, r = 0.5, repID = 2) {
  Hist_data = data.frame()
  # The size we re-measured, we can change this if future setting is changed!
  nc2String = c("nc1/4", "nc1/2", "nc1*3/4")
  for (i in 1:length(nc1Seq)) {
    nc1 = nc1Seq[i]
    nc2Seq = c(nc1/4, nc1/2, nc1*3/4)
    nt2 = nc1 
    tmp = data.frame()
    for (j in 1:length(nc2Seq)) {
      nc2 = nc2Seq[j]
      fileList = paste0("./simuData/S1-",method,"-",repID, "-", nc1, "-", nt2, "-",
                        nc2, "-", s1, "-", s2, "-", r, ".RData")
      load(fileList)
      a0List = Final$a0
      a0VarList = Final$a0Var
      # This is the statistics, we expected it to be std normal 
      a0Scale = (unlist(Final$a0) - a0)/sqrt(unlist(a0VarList))
      Density_data = data.frame(x = a0Scale)
      Density_data$method = method
      Density_data$a0 = a0
      Density_data$nc1 = nc1
      Density_data$nt2 = nt2 
      Density_data$nc2 = nc2String[j]
      tmp = rbind(tmp, Density_data)
    }
    Hist_data = rbind(Hist_data, tmp)
  }
  return(Hist_data)
}

## For example 
nc1Seq = c(100, 200, 500)
#method = "Batch2"
method = "ReMeasure"
Hist_data = Hist_Plot(method, nc1Seq)
p = ggplot(data = Hist_data, aes(x = x))   + 
  geom_histogram(aes(y = ..density..), bins = 30, colour = "black", fill = "blue") + 
  geom_density(col=2) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1) ) + facet_grid(nc1~nc2, labeller = label_both)
# add the mean vertical line 
p + geom_vline(aes(xintercept = mean(x)), color = "red", linetype = "dashed", size = 1)
