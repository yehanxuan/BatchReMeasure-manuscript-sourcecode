#### Plot The MSE ####
rm(list = ls())
require(reshape)
require(ggplot2)
require(cowplot)
library(latex2exp)
# computeMSE 
computeMSE = function(resultT, trueValue) {
  err = lapply(resultT, function(x) (x - trueValue)^2 )
  res = Reduce("+", err)/length(err)
  resultS2 = lapply(err, function(tmp) tmp^2)
  res2 = Reduce("+", resultS2)/length(err)
  resSD = sqrt(res2 - res^2)/sqrt(length(err))
  return(list(res = res, resSD = resSD))
}

a0s = c(0.5)
a1s = c(0.5)
r1s = c(0.5, 1, 2)
r2s = c(0.3, 0.6, 0.9)
n1s = seq(5, 50, length = 10)
n = 100
methods = c("ReMeasure", "Batch2", "Ignore", "LS")
niter = 1000

array.a = array(NA, c(length(a0s), length(a1s), length(r1s), length(r2s), length(n1s), length(methods), niter), 
                   dimnames = list(TrueEffect = paste(a0s), LocationEffect = paste(a1s), ScaleEffect = paste(r1s), SNR = paste(r2s), 
                                   RemeasureNo = paste(n1s), Method = methods, Iter = paste(1: niter)) )

for (m in 1:length(methods) ) {
  method = methods[m]
  for (i in 1:length(a0s) ) {
    a0 = a0s[i]
    for (j in 1:length(a1s)) {
      a1 = a1s[j]
      for (k in 1:length(r1s)) {
        r1 = r1s[k]
        for (s in 1:length(r2s)) {
          r2 = r2s[s]
          for (c in 1:length(n1s)) {
            n1 = n1s[c]
            fileList = paste0("./S1Data_Revise/S1-method-", method, "-repID-", 
                              2, "-nc1-", n/2, "-nt2-", n/2, "-nc2-", n1, "-r1-", r1, 
                              "-r2-", r2, "-a0-", a0, "-a1-", a1, ".RData")
            load(fileList)
            for (iter in 1:niter) {
              diff_a = Final$a0[[iter]] - a0
              array.a[as.character(a0), as.character(a1), as.character(r1), as.character(r2),
                         as.character(n1), m, as.character(iter)] = diff_a
            }
          }
        }
      }
    }
  }
}


create.MSE = function(array.a) {
  a0.df = melt(array.a)
  colnames(a0.df)[ncol(a0.df)] = 'Value'
  error.info = aggregate(Value ~ TrueEffect + LocationEffect + ScaleEffect + SNR + RemeasureNo + Method,
                          a0.df, function(x) mean(is.na(x)), na.action = function (x) x)
  
  m = aggregate(Value ~ TrueEffect + LocationEffect + ScaleEffect+ SNR + RemeasureNo + Method,
                 a0.df, function(x)  mean(x^2) ) ## Modify later 
  
  sd  = aggregate(Value ~ TrueEffect + LocationEffect + ScaleEffect + SNR +  RemeasureNo + Method,
                   a0.df, function(x) {
                     sd(x^2)/sqrt(niter)
                   })
  
  MSE.df = cbind(m, ymax=m[, ncol(m)] +  sd[, ncol(m)], ymin=m[, ncol(m)] - sd[, ncol(m)],  
                 SD=sd[, ncol(sd)], ErrRate=error.info[, ncol(m)])
  
  MSE.df$TrueEffect = factor(MSE.df$TrueEffect, levels = paste(a0s))
  MSE.df$LocationEffect= factor(MSE.df$LocationEffect, levels = paste(a1s))
  levels(MSE.df$LocationEffect) = TeX( paste('$a_{10}$:', levels(MSE.df$LocationEffect)) )
  MSE.df$ScaleEffect = factor(MSE.df$ScaleEffect, levels = paste(r1s))
  levels(MSE.df$ScaleEffect) = TeX(paste('$\\sigma_{1}$:', levels(MSE.df$ScaleEffect)))
  MSE.df$SNR = factor(MSE.df$SNR, levels = paste(r2s))
  MSE.df$RemeasureNo = factor(MSE.df$RemeasureNo, levels = paste(n1s))
  
  return(list(error = error.info, MSE.df = MSE.df))
}


array.a
result = create.MSE(array.a)
MSE.df = result$MSE.df

MSE.df$Method  = factor(MSE.df$Method, levels = c("ReMeasure", "Batch2", "Ignore", "LS"))

nMeth = length(methods) 
obj.list = list()
i = 1

for (r2 in paste(r2s)) {
  for (a0 in paste(a0s)) {
    
    res3 = subset(MSE.df, TrueEffect %in% a0 & SNR %in% r2, drop=TRUE)	
    res3$RemeasureNo = as.numeric(as.character(res3$RemeasureNo))
    
    dodge = position_dodge(width=0.9)
    
    obj = ggplot(res3, aes(x=RemeasureNo, y=Value,  group=Method, 
                            col = Method, shape = Method, linetype = Method)) +
      geom_line() +
      geom_point(size = 2) +scale_fill_discrete( labels = methods ) +
    				geom_errorbar(aes(ymax=ymax, ymin=ymin)) +
      scale_color_manual(values = c("red", "#E69F00", "#56B4E9", "#009E73")) 
    
    if (a0 == '0') {
      obj = obj  + 
        ylab('MSE') +
        facet_wrap(~ScaleEffect, labeller=label_parsed)
    } else {
      obj = obj + 
        ylab('MSE') +
        facet_wrap(~ScaleEffect, labeller=label_parsed)
    }
    
    obj = obj +
      xlab("No. of remeasured samples") +
      ggtitle(TeX(paste0('$a_0$=', a0, ' $a_1$=', a1, ' $\\rho$=', r2)) ) +
      theme_bw(base_size = 22) +
      theme(legend.position="bottom") +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_y_log10()
    
    obj.list[[i]] = obj
    i = i + 1
    print(obj) 
  }
}


GraphName = paste0("./figure/Fig2.pdf")
pdf(GraphName, width = 15, height = 20)

obj = plot_grid( obj.list[[1]], obj.list[[2]], obj.list[[3]],
                  labels = "AUTO", ncol = 1, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()





 




