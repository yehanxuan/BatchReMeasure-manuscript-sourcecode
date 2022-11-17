# Plot the NULL 
rm(list = ls())
library(ggplot2)
library(cowplot)
library(reshape)
library(stringi)
source("./code/OvarianCancer.R")
allMethods = c("ReMeasure", "ComBat", "sva")
Datalist = c("C1C2C4C5_RNAseq", "C1C2C4C5_Agilent", "C1C2C4C5_Agilent")
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

filepath = paste0("./data/IndependSvaComBatNULL-", paste(Datalist, collapse = '-'), ".RData")
load(filepath)

Ind_all = c( rbind( ReMeasure_C1MES_name[1:10], ReMeasure_C2IMM_name[1:10],
                    ReMeasure_C4DIF_name[1:10], ReMeasure_C5PRO_name[1:10]) )
n1s = seq(10, length(Ind_all), by = 5)
# Provide the histogram of p-value 
res_Re = res$pv
rownames(res_Re) <- n1s 
Re_frame= data.frame( melt(res_Re), "Methods" = "ReMeasure")
colnames(Re_frame) <- c("Size", "Gene", "Pvalue", "Methods")
res_Com = res$pv_par
rownames(res_Com) <- n1s 
Com_frame = data.frame( melt(res_Com), "Methods" = "ComBat")
colnames(Com_frame) <- c("Size", "Gene", "Pvalue", "Methods")
res_sva = res$pv_sva 
rownames(res_sva) <- n1s 
sva_frame = data.frame( melt(res_sva), "Methods" = "sva")
colnames(sva_frame) <- c("Size", "Gene", "Pvalue", "Methods")
hist_frame = rbind(Re_frame, Com_frame, sva_frame)

str = stri_replace_all_regex(Datalist,
                             pattern =c("MES", "IMM", "PRO", "DIF", '_'),
                             replacement = c('','','','', ' '),
                             vectorize = FALSE)
str = sub("(C\\d)(C\\d)(C\\d)", "\\1+\\2+\\3+", str)

l = 3
n1 = n1s[l]
df.hist = subset(hist_frame, Size %in% n1)
df.hist$Methods <- factor(df.hist$Methods, levels = c("ReMeasure", "ComBat", "sva", "LMM"))
pobj = ggplot() + geom_histogram(data = df.hist , aes(x = Pvalue),breaks = seq(0,1,length.out = 30)) +  
  theme_bw(base_size = 16) + xlab("P-value") + ylab("Count") + 
  theme(panel.grid = element_blank()) + facet_wrap(~Methods, scales='free', nrow = 1) + 
  theme(strip.background =  element_blank()) + 
  ggtitle(gsub("_", " ", paste(str[-3], collapse = ' vs. ') ) ) + 
  theme(plot.title = element_text(size = 15, hjust = 0.5))

GraphName = paste0("./figure/FigS18a.pdf")
pdf(GraphName, width = 12, height = 4)
pobj
dev.off()

### non-NULL ###
DList = list(
  c("C1C2C5_RNAseq", "C4DIF_Agilent", "C1C2C5_Agilent"),
  c("C1C4C5_RNAseq", "C2IMM_Agilent", "C1C4C5_Agilent"),
  c("C1C2C4_RNAseq", "C5PRO_Agilent", "C1C2C4_Agilent"),
  c("C2C4C5_RNAseq", "C1MES_Agilent", "C2C4C5_Agilent")
)


  Datalist = DList[[4]]
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
 
  
res_Re = res$pv
rownames(res_Re) <- n1s 
Re_frame= data.frame( melt(res_Re), "Methods" = "ReMeasure")
colnames(Re_frame) <- c("Size", "Gene", "Pvalue", "Methods")
res_Com = res$pv_par
rownames(res_Com) <- n1s 
Com_frame = data.frame( melt(res_Com), "Methods" = "ComBat")
colnames(Com_frame) <- c("Size", "Gene", "Pvalue", "Methods")
res_sva = res$pv_sva 
rownames(res_sva) <- n1s 
sva_frame = data.frame( melt(res_sva), "Methods" = "sva")
colnames(sva_frame) <- c("Size", "Gene", "Pvalue", "Methods")
hist_frame = rbind(Re_frame, Com_frame, sva_frame)

l = 6 
n1 = n1s[l]
df.hist = subset(hist_frame, Size %in% n1)
df.hist$Methods <- factor(df.hist$Methods, levels = c("ReMeasure", "ComBat", "sva", "LMM"))
pobj = ggplot() + geom_histogram(data = df.hist , aes(x = Pvalue),breaks = seq(0,1,length.out = 30)) +  
  theme_bw(base_size = 16) + xlab("P-value") + ylab("Count") + 
  theme(panel.grid = element_blank()) + facet_wrap(~Methods, scales='free', nrow = 1) + 
  theme(strip.background =  element_blank()) + 
  ggtitle(gsub("_", " ", paste(str[-3], collapse = ' vs. ') ) ) + 
  theme(plot.title = element_text(size = 15, hjust = 0.5))

GraphName = paste0("./figure/FigS18b.pdf")
pdf(GraphName, width = 12, height = 4)
pobj
dev.off()

