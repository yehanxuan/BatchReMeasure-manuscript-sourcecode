### Fig S16 ###
rm(list = ls())
library(ggplot2)
library(cowplot)
library(reshape)
library(stringi)
source("./code/OvarianCancer.R")
Datalist = c("C1C2C4C5_RNAseq", "C1C2C4C5_Agilent", "C1C2C4C5_Agilent")
allMethods = c("ReMeasure", "ComBat", "sva")
create.data = function(res.pv.a, procedure) {
  res.df <- reshape::melt(res.pv.a)
  colnames(res.df)[ncol(res.df)] <- 'Value'
  if (procedure == "none") {
    m <- aggregate(Value~Method + ReMeasureNo, res.df, function(x) sum(x <= alpha ))
  } else if (procedure == "BH") {
    m <- aggregate(Value~Method + ReMeasureNo, res.df, function(x) sum( p.adjust(x, method = "BH") < alpha ) )
  } else if (procedure == "bonferroni") {
    m <- aggregate(Value~Method + ReMeasureNo, res.df, function(x) sum( p.adjust(x, method = "bonferroni") < alpha ) )
  }
  m$ReMeasureNo = factor(m$ReMeasureNo)
  return(m)
}

obj.list = list()
procedure = "none"
nGene = 11861
Ind_all = c( rbind( ReMeasure_C1MES_name[1:10], ReMeasure_C2IMM_name[1:10],
                    ReMeasure_C4DIF_name[1:10], ReMeasure_C5PRO_name[1:10]) )
filepath = paste0("./data/IndependSvaComBatNULL-", paste(Datalist, collapse = '-'), ".RData")


load(filepath)
n1s = seq(10, length(Ind_all), by = 5)
res.pv.a = array(NA, c(length(allMethods), length(n1s), nGene), 
                 dimnames = list(Method = allMethods, ReMeasureNo = paste(n1s), 
                                 Gene = paste(1:nGene) ) )

for (m in 1:length(allMethods)) {
  method = allMethods[m]
  if (method == "ReMeasure") {
    pvMat = res$pv
  } else if (method == "ComBat") {
    pvMat = res$pv_par
  } else if (method == "sva") {
    pvMat = res$pv_sva
  } 
  
  for (num in 1:length(n1s)) {
    n1 = n1s[num]
    res.pv.a[method, as.character(n1), ] = pvMat[num, ] 
  }
  
}
error.info <- aggregate(value~Method + ReMeasureNo,
                        reshape::melt(res.pv.a), function(x) mean(is.na(x)), na.action = function (x) x)

alpha = 0.05
res.df = create.data(res.pv.a, procedure = procedure)
res.df$Method <- factor(res.df$Method, levels = c("ReMeasure", "ComBat", "sva", "LMM"))
str = stri_replace_all_regex(Datalist,
                             pattern =c("MES", "IMM", "PRO", "DIF", '_'),
                             replacement = c('','','','', ' '),
                             vectorize = FALSE)
str = sub("(C\\d)(C\\d)(C\\d)", "\\1+\\2+\\3+", str)
res.sub = subset(res.df, Method %in% c("ReMeasure", "ComBat", "sva", "LMM"))

GraphName = paste0("./figure/FigS16.pdf")
pdf(GraphName, width = 10, height = 8)
obj <- ggplot(res.sub, aes(x = ReMeasureNo, y = Value, group = Method,
                           col = Method, shape = Method,  linetype = Method)) + 
  geom_line(size = 1) + geom_point(size = 2) 
obj <- obj + xlab("No. of remeasured samples") + theme_bw(base_size = 22) + 
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1)) + 
  ggtitle(gsub("_", " ", paste(str[-3], collapse = ' vs. ') ) ) + 
  theme(plot.title = element_text(size = 15, hjust = 0.5)) + ylab("No. of discoveries") +
  theme(legend.title = element_blank(), panel.grid.minor = element_blank()) 
obj <- obj  + geom_hline(aes(yintercept=alpha*nGene), linetype=2)
print(obj)
dev.off()




