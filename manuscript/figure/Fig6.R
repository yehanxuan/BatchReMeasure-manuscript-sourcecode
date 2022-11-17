### plot Ovarian Rank ###
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

pRank.list = list()
allMethods = c("ReMeasure", "Batch2", "Ignore", "LS")
nGene = 11861
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

  ### Golden Standard ###
  T_pvVec = NULL
  for (i in 1:nrow(recur_Agilent)) {
    if ( isTRUE(grepl("recur", Datalist[1], fix = TRUE)) ) {
      xData = recur_Agilent
      yData = non_Agilent
    } else if ( isTRUE(grepl("non", Datalist[1], fix = TRUE)) ) {
      xData = non_Agilent
      yData = recur_Agilent
    } else if ( isTRUE(grepl("C1MES", Datalist[1], fix = TRUE)) ) {
      xData = C1MES_Agilent
      yData = C2C4C5_Agilent
    }  else if ( isTRUE(grepl("C2IMM", Datalist[1], fix = TRUE)) ) {
      xData = C2IMM_Agilent
      yData = C1C4C5_Agilent
    } else if ( isTRUE(grepl("C4DIF", Datalist[1], fix = TRUE)) ) {
      xData = C4DIF_Agilent
      yData = C1C2C5_Agilent
    } else if ( isTRUE(grepl("C5PRO", Datalist[1], fix = TRUE)) ) {
      xData = C5PRO_Agilent
      yData = C1C2C4_Agilent
    } else if ( isTRUE(grepl("C1C2C4", Datalist[1], fix = TRUE)) ) {
      xData = C1C2C4_Agilent
      yData = C5PRO_Agilent
    } else if ( isTRUE(grepl("C1C2C5", Datalist[1], fix = TRUE)) ) {
      xData = C1C2C5_Agilent
      yData = C4DIF_Agilent
    } else if ( isTRUE(grepl("C1C4C5", Datalist[1], fix = TRUE)) ) {
      xData = C1C4C5_Agilent
      yData = C2IMM_Agilent
    } else if (isTRUE(grepl("C2C4C5", Datalist[1], fix = TRUE)) ) {
      xData = C2C4C5_Agilent
      yData = C1MES_Agilent
    }
    x = xData[i, ]
    y = yData[i, ]
    if ( (0.95 < var(x)/var(y) ) && ( var(x)/var(y) < 1.05 ) ) {
      Equal = TRUE
    } else {
      Equal = FALSE
    }

    Ttest = t.test(x, y, mu = 0, alternative = "two.sided", var.equal = Equal, paired = F)
    T_pvVec = c(T_pvVec, Ttest$p.value)
  }
  alpha_t = 0.05
  procedure = "BH"
  Rej_t_Ind = (p.adjust(T_pvVec, method = procedure) < alpha_t)
  sum(Rej_t_Ind)
  # significant genes
  Sig_Ind = which(Rej_t_Ind == TRUE)
  cat(length(Sig_Ind))
  n1s = seq(10, length(Ind_all), by = 5)
  res.Rank.a = array(NA, c(length(allMethods), length(n1s)), 
                                      dimnames = list(Method = allMethods, ReMeasureNo = paste(n1s)
                                                       ) )
  
  for (m in 1:length(allMethods)) {
    method = allMethods[m]
    result = read.csv(paste0( './Data_real/', paste(Datalist, collapse = '-'), "-method-", method,
                              ".csv"), header = TRUE, check.names = FALSE)
    
    for (n1 in n1s) {
        sort_result = sort(result[, paste(n1)], decreasing =TRUE, index.return = TRUE)$ix
        res.Rank.a[method, as.character(n1)] = mean( which( sort_result %in% Sig_Ind) )
    }
  }
  df.rank = melt(res.Rank.a)
  str = stri_replace_all_regex(Datalist,
                               pattern =c("MES", "IMM", "PRO", "DIF", '_'),
                               replacement = c('','','','', ' '),
                               vectorize = FALSE)
  str = sub("(C\\d)(C\\d)", "\\1+\\2+", str)
  
  df.rank$Method = factor(df.rank$Method, levels = c("ReMeasure", "Batch2", "Ignore", "LS"))
  df.rank$ReMeasureNo  <-  as.numeric(as.character(df.rank$ReMeasureNo))
  pRank = ggplot(df.rank, aes(x = ReMeasureNo, y = value, group = Method, 
                              col = Method, shape = Method,  linetype = Method)) + 
    geom_line(size = 1) + geom_point(size = 2) + xlab("No. of remeasured samples") + theme_bw(base_size = 16) +
    theme() +
    theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1)) + 
    ggtitle( gsub("_", " ", paste(str[-3], collapse = ' vs. ') )  ) + 
    theme(plot.title = element_text(size = 15, hjust = 0.5)) + ylab("Average rank")
  pRank.list[[l]] <- pRank
}


GraphName = paste0("./figure/Fig6.pdf")
pdf(GraphName, width = 12, height = 10)
plot_grid(pRank.list[[1]], pRank.list[[2]], 
          pRank.list[[3]], pRank.list[[4]], ncol = 2)
dev.off()


