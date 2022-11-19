### Theoretical curve ####
rm(list = ls())
library(ggplot2)
library(reshape)
require(cowplot)
library(latex2exp)
method = "ReMeasure"
ns = c(100, 200, 400)
a1s=a1= 0.5
a0s = c( 0.5, 0.8)
RecoverSeq = c(0.8, 0.9, 0.95)
r2s = seq(0, 0.9, len = 10)
r1s = c( 0.5, 1,2 )
alpha = 0.05

#load("./Oracle_S1/Ratio-ReMeasure-repID-2-a1-0.5-a0-0.5-0.8.RData")
load("./data/Ratio-ReMeasure-repID-2-a1-0.5-a0-0.5-0.8.RData")

sample.df = reshape2::melt(ratio.a)
colnames(sample.df)[ncol(sample.df)] <- 'ReMeasureNo'
sample.df$TrueEffect <- factor(sample.df$TrueEffect, levels = paste(a0s))
sample.df$LocationEffect<- factor(sample.df$LocationEffect, levels = paste(a1s))
sample.df$samples <- factor(2*sample.df$samples, levels = paste(ns))
levels(sample.df$LocationEffect) <- TeX( paste('$a_{10}$:', levels(sample.df$LocationEffect)) )
sample.df$ScaleEffect <- factor(sample.df$ScaleEffect, levels = paste(r1s))
levels(sample.df$samples) <- TeX( paste0('$n=$', levels(sample.df$samples)) )
sample.df$Ratio <- factor(sample.df$Ratio)

obj.list <- list()
i = 1
r1 = 2
for (a0 in paste(a0s)) {
    res <- subset(sample.df, TrueEffect %in% a0 & ScaleEffect %in% r1 , drop = TRUE)
    dodge <- position_dodge(width=0.9)
    obj = ggplot(res, aes(x = SNR, y=ReMeasureNo, group = Ratio, col = Ratio, 
                          shape = Ratio, linetype = Ratio)) +
      geom_point(size = 4) + geom_line()
    
    obj <- obj  +
      facet_grid(~ samples, labeller=label_parsed)
    
    obj <- obj +
      xlab("Correlation") + ylab("Proportion of samples remeasured") + 
      ggtitle(TeX( paste0('Cohen\'s d=', a0 ) ) ) +
      theme_bw(base_size = 22) +
      theme(legend.position="bottom") +
      theme(axis.text.x = element_text(angle = 45, hjust=1, size = 25),
            axis.text.y = element_text(size = 25))   +
      theme(legend.text = element_text(size = 20),
            legend.key.size = unit(2, "cm")
            ) + guides(col= guide_legend(title= "Relative power"),
                       shape =guide_legend(title= "Relative power"),
                       linetype = guide_legend(title= "Relative power"))
    
    obj.list[[i]] <- obj
    i <- i + 1
    print(obj) 
}

pdf("./figure/Fig4.pdf", width = 15, height = 16)
obj <- plot_grid( obj.list[[1]], obj.list[[2]], 
                  labels = "AUTO", ncol = 1, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()


