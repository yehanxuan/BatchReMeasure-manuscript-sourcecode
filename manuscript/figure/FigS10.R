### Figure S10 ###
rm(list = ls())
source("./code/OvarianCancer.R")
library(ggplot2)
library(reshape)
library(cowplot)
require(latex2exp)
library(stringi)
DList = list(
  c("C1C2C5_RNAseq", "C4DIF_Agilent", "C1C2C5_Agilent"),
  c("C1C4C5_RNAseq", "C2IMM_Agilent", "C1C4C5_Agilent"),
  c("C1C2C4_RNAseq", "C5PRO_Agilent", "C1C2C4_Agilent"),
  c("C2C4C5_RNAseq", "C1MES_Agilent", "C2C4C5_Agilent")
)

allMethods = c("ReMeasure", "Batch2")
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


strength = "Strong"
nGene = length( get( paste0("Gene_", strength)) )
obj.list = list()
  for (l in 1:length(DList)) {
    Datalist = DList[[l]]
    if (isTRUE( grepl("recur", Datalist[1], fix=TRUE) ) ) {
      Ind_all = ReMeasure_recur_name
    } else if ( isTRUE(grepl("non",Datalist[1], fix = TRUE ) ) ) {
      Ind_all = ReMeasure_non_name
    } else if ( isTRUE(grepl("C1MES", Datalist[1], fix = TRUE)) ) {
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
    
    n1s = seq(10, length(Ind_all), by = 5)
    res.pv.a = array(NA, c(length(allMethods), length(n1s), nGene), 
                     dimnames = list(Method = allMethods, ReMeasureNo = paste(n1s), 
                                     Gene = paste(1:nGene) ) )
    
    for (m in 1:length(allMethods)) {
      method = allMethods[m]
      result = read.csv(paste0( './data/strength-', strength, "-", paste(Datalist, collapse = '-'), "-method-", method,
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
    res.df$Method = factor(res.df$Method, levels = c("ReMeasure", "Batch2"))
    
    str = stri_replace_all_regex(Datalist,
                                 pattern =c("MES", "IMM", "PRO", "DIF", '_'),
                                 replacement = c('','','','', ' '),
                                 vectorize = FALSE)
    str = sub("(C\\d)(C\\d)", "\\1+\\2+", str)
    
    obj <- ggplot(res.df, aes(x = ReMeasureNo, y = Value, group = Method, 
                              col = Method, shape = Method,  linetype = Method)) +
      geom_line(size = 1) + geom_point(size = 2) 
    obj <- obj + xlab("No. of remeasured samples") + theme_bw(base_size = 16) + 
      theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1)) + 
      ggtitle(gsub("_", " ", paste(str[-3], collapse = ' vs. ') ) ) + 
      theme(plot.title = element_text(size = 15, hjust = 0.5)) + ylab("No. of discoveries") + 
      scale_color_manual(values = c("#F8766D", "#7CAE00"))
    
    obj.list[[l]] <- obj 
    
  }



GraphName = paste0("./figure/FigS10.pdf")
pdf(GraphName, width = 12, height = 10)
plot_row <- plot_grid(obj.list[[1]], obj.list[[2]], 
          obj.list[[3]], obj.list[[4]], ncol = 2)
# Joint plot titles
title <- ggdraw() + draw_label(paste(strength, "correlation"), fontface = 'bold',
                               x = 0.1,
                               hjust = 0, size = 20) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1, 1))
dev.off()


