### real data Null hypothesis ###
rm(list = ls())
source("./code/OvarianCancer.R")
library(stringi)
library(ggplot2)
allMethods = c("ReMeasure", "Batch2", "Ignore", "LS")
Datalist = c("C1C2C4C5_RNAseq", "C1C2C4C5_Agilent", "C1C2C4C5_Agilent")
Ind_all = c( rbind( ReMeasure_C1MES_name[1:10], ReMeasure_C2IMM_name[1:10],
       ReMeasure_C4DIF_name[1:10], ReMeasure_C5PRO_name[1:10]) )


n1s <- seq(5, length(Ind_all), by = 5)
nGene = 11861
alpha = 0.05
procedure = "none"
create.data.discover = function(res.pv.a, procedure = "none") {
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


res.pv.a = array(NA, c(length(allMethods), length(n1s), nGene), 
                 dimnames = list(Method = allMethods, ReMeasureNo = paste(n1s), 
                                 Gene = paste(1:nGene) ) )

for (m in 1:length(allMethods)) {
  method = allMethods[m]
  result = read.csv(paste0( './data/order-', paste(Datalist, collapse = '-'), "-method-", method,
                            ".csv"), header = TRUE, check.names = FALSE)
  for (n1 in n1s) {
    for (gene in 1:nGene) {
      res.pv.a[method, as.character(n1), as.character(gene)] = result[gene, paste(n1)]
    }
  }
}

res.pv.a
res.df  = create.data.discover(res.pv.a, procedure)
res.df$ReMeasureNo <-  as.numeric(as.character(res.df$ReMeasureNo))
res.df$Method = factor(res.df$Method, levels = c("ReMeasure", "Batch2", "Ignore", "LS"))

str = stri_replace_all_regex(Datalist,
                             pattern =c("MES", "IMM", "PRO", "DIF", '_'),
                             replacement = c('','','','', ' '),
                             vectorize = FALSE)
str = sub("(C\\d)(C\\d)(C\\d)", "\\1+\\2+\\3+", str)

GraphName = paste0("./figure/Fig5.pdf")
pdf(GraphName, width = 10, height = 8)
obj <- ggplot(res.df, aes(x = ReMeasureNo, y = Value, group = Method, 
                          col = Method, shape = Method,  linetype = Method)) +
  geom_line(size = 1) + geom_point(size = 2) 

obj <- obj + geom_hline(aes(yintercept=alpha*nGene), linetype=2)
obj <- obj + xlab("No. of remeasured samples") + theme_bw(base_size = 22) + 
  theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1)) + 
  ggtitle(gsub("_", " ", paste(str[-3], collapse = ' vs. ') ) ) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + ylab("No. of discoveries") 
print(obj)
dev.off()






