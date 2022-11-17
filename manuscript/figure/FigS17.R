### Fig S17 ###
rm(list = ls())
library(ggplot2)
library(cowplot)
library(reshape)
library(stringi)
source("./code/OvarianCancer.R")
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


nGene = 11861
DList = list(
  c("C1C2C5_RNAseq", "C4DIF_Agilent", "C1C2C5_Agilent"),
  c("C1C4C5_RNAseq", "C2IMM_Agilent", "C1C4C5_Agilent"),
  c("C1C2C4_RNAseq", "C5PRO_Agilent", "C1C2C4_Agilent"),
  c("C2C4C5_RNAseq", "C1MES_Agilent", "C2C4C5_Agilent")
)

obj.list = list()
procedure = "BH"
for (l in 1:length(DList)) { 
  Datalist = DList[[l]]
  ds1 = get(Datalist[1])
  nGene = nrow(ds1)
  if ( isTRUE(grepl("C1MES", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C1MES_name
  } else if ( isTRUE(grepl("C2IMM", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C2IMM_name
  } else if ( isTRUE(grepl("C4DIF", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C4DIF_name
  } else if ( isTRUE(grepl("C5PRO", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C5PRO_name
  } else if ( isTRUE(grepl("C1C2C4", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C1C2C4_name
  } else if ( isTRUE(grepl("C1C2C5", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C1C2C5_name
  } else if ( isTRUE(grepl("C1C4C5", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C1C4C5_name
  } else if (isTRUE(grepl("C2C4C5", Datalist[1], fix = TRUE)) ) {
    Ind_all = ReMeasure_C2C4C5_name
  }
  
  filepath = paste0("./data/IndependSvaComBat-", paste(Datalist, collapse = '-'), ".RData")
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
  
  alpha = 0.05
  res.df = create.data(res.pv.a, procedure = procedure)
  res.df$Method <- factor(res.df$Method, levels = c("ReMeasure", "ComBat", "sva", "LMM"))
  str = stri_replace_all_regex(Datalist,
                               pattern =c("MES", "IMM", "PRO", "DIF", '_'),
                               replacement = c('','','','', ' '),
                               vectorize = FALSE)
  str = sub("(C\\d)(C\\d)(C\\d)", "\\1+\\2+\\3+", str)
  obj <- ggplot(res.df, aes(x = ReMeasureNo, y = Value+1, group = Method,
                            col = Method, shape = Method,  linetype = Method)) + 
    geom_line(size = 1) + geom_point(size = 2) 
  obj <- obj + xlab("No. of remeasured samples") + theme_bw(base_size = 16) + 
    theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1)) + 
    ggtitle(gsub("_", " ", paste(str[-3], collapse = ' vs. ') ) ) + 
    theme(plot.title = element_text(size = 15, hjust = 0.5)) + ylab("No. of discoveries") +
    theme(legend.title = element_blank(), panel.grid.minor = element_blank()) 
  obj.list[[l]] <- obj 
}

GraphName = paste0("./figure/FigS17.pdf")
pdf(GraphName, width = 12, height = 10)
plot_grid(obj.list[[1]], obj.list[[2]], 
          obj.list[[3]], obj.list[[4]], ncol = 2)
dev.off()



