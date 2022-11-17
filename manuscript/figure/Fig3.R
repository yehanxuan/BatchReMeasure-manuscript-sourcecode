rm(list = ls())
require(reshape)
require(ggplot2)
require(cowplot)
require(latex2exp)
a0s <- c(0, 0.5, 0.8)
a1s = c(0.5)
r1s <- c(2)
r2s <- c(0.3, 0.6, 0.9)
n = 100
n1s <- seq(5, 50, len = 10)
methods = c("ReMeasure", "Batch2", "Ignore", "LS")
allMethods = c("ReMeasure", "Batch2", "Ignore", "LS")
niter = 1000

res.Sim1.a = array(NA, c(length(a0s), length(a1s), length(r1s), length(r2s), length(n1s), length(methods), niter), 
                    dimnames = list(TrueEffect = paste(a0s), LocationEffect = paste(a1s), ScaleEffect = paste(r1s), SNR = paste(r2s), 
                                    RemeasureNo = paste(n1s), Method = allMethods, Iter = paste(1: niter)) )
                              
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
                
                if (method == "ReMeasure" | method == "Batch2" | method == "Ignore" | method == "Gen") {
                  pv <- 2 * stats::pnorm(-abs(Final$a0[[iter]] / sqrt(Final$a0Var[[iter]]) ) )   
                } else if (method == "LS") {
                  pv <- Final$p.value[[iter]]
                  }else {
                 pv <- mean(abs(Final$ztestb[[iter]]) > abs(Final$ztest[[iter]]) )
                }
                res.Sim1.a[as.character(a0), as.character(a1), as.character(r1), as.character(r2),
                           as.character(n1), m, as.character(iter)] = pv
              }
            }
          }
        }
      }
    }
  }

create.data <- function (res.a) {
  # Calculate the standard error
  res.df <- reshape::melt(res.a)
  colnames(res.df)[ncol(res.df)] <- 'Value'
  # Error counts
  error.info <- aggregate(Value ~ TrueEffect + LocationEffect + ScaleEffect + SNR + RemeasureNo + Method,
                          res.df, function(x) mean(is.na(x)), na.action = function (x) x)
  
  m  <- aggregate(Value ~ TrueEffect + LocationEffect + ScaleEffect+ SNR + RemeasureNo + Method,
                  res.df, function(x) mean(x <= alpha))
  se  <- aggregate(Value ~ TrueEffect + LocationEffect + ScaleEffect + SNR +  RemeasureNo + Method,
                   res.df, function(x) {
                     p <- mean(x <= alpha)
                     sqrt(p * (1 - p) / length(x))
                   })
  
  
  
  res.df2 <- cbind(m, ymax=m[, ncol(m)] + 1.96 * se[, ncol(m)], ymin=m[, ncol(m)] - 1.96 * se[, ncol(m)],  
                   SE=se[, ncol(se)], ErrRate=error.info[, ncol(m)])
  res.df2$TrueEffect <- factor(res.df2$TrueEffect, levels = paste(a0s))
  res.df2$LocationEffect<- factor(res.df2$LocationEffect, levels = paste(a1s))
  levels(res.df2$LocationEffect) <- TeX(paste('$a_{1}$:', levels(res.df2$LocationEffect)) )
  res.df2$ScaleEffect <- factor(res.df2$ScaleEffect, levels = paste(r1s))
  levels(res.df2$ScaleEffect) <- TeX(paste('$\\sigma_{1}$:', levels(res.df2$ScaleEffect)))
  res.df2$SNR <- factor(res.df2$SNR, levels = paste(r2s))
  res.df2$RemeasureNo <- factor(res.df2$RemeasureNo, levels = paste(n1s))
  
  return(list(error = error.info, res.df = res.df2))
}

##############
res.Sim1.a
alpha <- 0.05
result <- create.data(res.Sim1.a)
res.df <- result$res.df
write.csv(result$error, file = 'Sim1.error.csv')
res.df$Method = factor(res.df$Method, levels = c("ReMeasure", "Batch2", "Ignore", "LS"))
obj.list <- list()
i <- 1
levels(res.df$SNR) = TeX(paste('$\\rho$:', levels(result$res.df$SNR)))

for (a0 in paste(a0s)) {
    res3 <- subset(res.df, TrueEffect %in% a0, drop=TRUE)	
    res3$RemeasureNo <- as.numeric(as.character(res3$RemeasureNo))
    
    dodge <- position_dodge(width=0.9)
    
    obj <- ggplot(res3, aes(x=RemeasureNo, y=Value,  group=Method, 
                            col = Method, shape = Method, linetype = Method)) +
      geom_line() +
      geom_point(size = 2) +scale_fill_discrete( labels = allMethods ) + 
      scale_color_manual(values = c("red", "#E69F00", "#56B4E9", "#009E73"))
    
    if (a0 == '0') {
      obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
        ylab('Type I error') +
        facet_grid( ~ SNR, labeller=label_parsed)
    } else {
      obj <- obj + 
        ylab('Power') +
        facet_grid(~ SNR, scales='free', labeller=label_parsed)
    }
    
    obj <- obj +
      xlab("No. of remeasured samples") + 
      ggtitle(TeX( paste0('$a_0$=', a0, ' $a_1$=', a1, ' $\\sigma_{1}$ =', r1) ) ) +
      theme_bw(base_size = 22) +
      theme(legend.position="bottom") +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) 
    
    obj.list[[i]] <- obj
    i <- i + 1
    print(obj) 
}


GraphName = "./manuscript/figure/Fig3.pdf"
pdf(GraphName, width = 15, height = 20)
obj <- plot_grid( obj.list[[1]], obj.list[[2]], obj.list[[3]],
                  labels = "AUTO", ncol = 1, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()

