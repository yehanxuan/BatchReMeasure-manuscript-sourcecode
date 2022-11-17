## Real data ##
rm(list = ls())
load("./data/Ovarian.Rdata")
Agilent.dat # Agilent batch 2 
RNAseq.dat # RNAseq batch 1

Agilent.recur
RNAseq.recur
remeasure.recur
unique(Agilent.recur)
sum(Agilent.recur == "yes")
sum(Agilent.recur == "no")
sum(RNAseq.recur == "yes")
sum(RNAseq.recur == "no")
sum(remeasure.recur == "yes")
sum(remeasure.recur == "no")
unique(RNAseq.recur)
unique(remeasure.recur)

Agilent.subtype
RNAseq.subtype
remeasure.subtype
unique(Agilent.subtype)
unique(RNAseq.subtype)
unique(remeasure.subtype)
table(Agilent.subtype)
table(RNAseq.subtype)
table(remeasure.subtype)

dim(RNAseq.dat)
dim(Agilent.dat)

colnames(RNAseq.dat)
rownames(RNAseq.recur)
rownames(RNAseq.subtype)

# recurrence data 
Agilent_recur_name =  colnames(Agilent.dat[ ,which( (Agilent.recur == "yes") %in% TRUE) ])
RNAseq_recur_name = colnames(RNAseq.dat[, which( (RNAseq.recur == "yes") %in% TRUE)])
ReMeasure_recur_name = rownames(remeasure.recur[which( (remeasure.recur == "yes") %in% TRUE), , drop = F])
recur_Agilent = Agilent.dat[, Agilent_recur_name]
recur_RNAseq = RNAseq.dat[ , RNAseq_recur_name]
dim(recur_Agilent)
dim(recur_RNAseq)
# non-recurrence data 
Agilent_non_name = colnames(Agilent.dat[ ,which( (Agilent.recur == "no") %in% TRUE) ])
RNAseq_non_name = colnames(RNAseq.dat[, which( (RNAseq.recur == "no") %in% TRUE)])
ReMeasure_non_name = rownames(remeasure.recur[which( (remeasure.recur == "no") %in% TRUE), , drop = F])
non_Agilent = Agilent.dat[ , Agilent_non_name]
non_RNAseq = RNAseq.dat[ , RNAseq_non_name]
dim(non_Agilent)
dim(non_RNAseq)
# C1-MES
Agilent_C1MES_name = colnames(Agilent.dat[ ,which( (Agilent.subtype == "C1-MES") %in% TRUE) ])
RNAseq_C1MES_name = colnames( RNAseq.dat[, which( (RNAseq.subtype == "C1-MES") %in% TRUE)])
ReMeasure_C1MES_name = rownames(remeasure.subtype[which( (remeasure.subtype == "C1-MES") %in% TRUE), , drop = F])
C1MES_Agilent = Agilent.dat[ , Agilent_C1MES_name]
C1MES_RNAseq = RNAseq.dat[ , RNAseq_C1MES_name]
dim(C1MES_Agilent)
dim(C1MES_RNAseq)
# C2-IMM
Agilent_C2IMM_name = colnames(Agilent.dat[ ,which( (Agilent.subtype == "C2-IMM") %in% TRUE) ])
RNAseq_C2IMM_name = colnames(RNAseq.dat[, which( (RNAseq.subtype == "C2-IMM") %in% TRUE)])
ReMeasure_C2IMM_name = rownames(remeasure.subtype[which( (remeasure.subtype == "C2-IMM") %in% TRUE), , drop = F])
C2IMM_Agilent = Agilent.dat[ , Agilent_C2IMM_name]
C2IMM_RNAseq =  RNAseq.dat[ , RNAseq_C2IMM_name]
dim(C2IMM_Agilent)
dim(C2IMM_RNAseq)
# C4-DIF 
Agilent_C4DIF_name = colnames(Agilent.dat[ ,which( (Agilent.subtype == "C4-DIF") %in% TRUE) ])
RNAseq_C4DIF_name = colnames(RNAseq.dat[, which( (RNAseq.subtype == "C4-DIF") %in% TRUE)])
ReMeasure_C4DIF_name = rownames(remeasure.subtype[which( (remeasure.subtype == "C4-DIF") %in% TRUE), , drop = F])
C4DIF_Agilent =  Agilent.dat[ , Agilent_C4DIF_name]
C4DIF_RNAseq = RNAseq.dat[ , RNAseq_C4DIF_name]
# C5-PRO 
Agilent_C5PRO_name = colnames(Agilent.dat[ , which( (Agilent.subtype == "C5-PRO") %in% TRUE) ])
RNAseq_C5PRO_name = colnames(RNAseq.dat[, which( (RNAseq.subtype == "C5-PRO") %in% TRUE)])
ReMeasure_C5PRO_name = rownames(remeasure.subtype[which( (remeasure.subtype == "C5-PRO") %in% TRUE), , drop = F])
C5PRO_Agilent = Agilent.dat[ , Agilent_C5PRO_name]
C5PRO_RNAseq = RNAseq.dat[ , RNAseq_C5PRO_name]


C2C4C5_Agilent = cbind(C2IMM_Agilent, C4DIF_Agilent, C5PRO_Agilent)
C2C4C5_RNAseq = cbind(C2IMM_RNAseq, C4DIF_RNAseq, C5PRO_RNAseq)
dim(C2C4C5_Agilent)
dim(C2C4C5_RNAseq)
ReMeasure_C2C4C5_name = c(ReMeasure_C2IMM_name, ReMeasure_C4DIF_name, ReMeasure_C5PRO_name)
length(ReMeasure_C2C4C5_name)

C1C2C5_Agilent = cbind(C1MES_Agilent, C2IMM_Agilent, C5PRO_Agilent)
C1C2C5_RNAseq = cbind(C1MES_RNAseq, C2IMM_RNAseq, C5PRO_RNAseq)
dim(C1C2C5_Agilent)
dim(C1C2C5_RNAseq)
ReMeasure_C1C2C5_name = c(ReMeasure_C1MES_name, ReMeasure_C2IMM_name, ReMeasure_C5PRO_name)
length(ReMeasure_C1C2C5_name)

C1C2C4_Agilent = cbind(C1MES_Agilent, C2IMM_Agilent, C4DIF_Agilent)
C1C2C4_RNAseq = cbind(C1MES_RNAseq, C2IMM_RNAseq, C4DIF_RNAseq)
dim(C1C2C4_Agilent)
dim(C1C2C4_RNAseq)
ReMeasure_C1C2C4_name = c(ReMeasure_C1MES_name, ReMeasure_C2IMM_name, ReMeasure_C4DIF_name)
length(ReMeasure_C1C2C4_name)

C1C4C5_Agilent = cbind(C1MES_Agilent, C4DIF_Agilent, C5PRO_Agilent)
C1C4C5_RNAseq = cbind(C1MES_RNAseq, C4DIF_RNAseq, C5PRO_RNAseq)
dim(C1C4C5_Agilent)
dim(C1C4C5_RNAseq)
ReMeasure_C1C4C5_name = c(ReMeasure_C1MES_name, ReMeasure_C4DIF_name, ReMeasure_C5PRO_name)
length(ReMeasure_C1C4C5_name)


Agilent_Re = Agilent.dat[ , rownames(remeasure.recur)]
RNAseq_Re = RNAseq.dat[ , rownames(remeasure.recur)]
dim(Agilent_Re)
dim(RNAseq_Re)
CorSeq = rep(0, nrow(Agilent_Re))
for (g in 1:nrow(Agilent_Re)) {
  CorSeq[g] = cor(Agilent_Re[g, ], RNAseq_Re[g, ])  
}
CorSeq_abs = abs(CorSeq)
q_low = quantile(CorSeq_abs, 1/4)
q_high = quantile(CorSeq_abs, 3/4)
Gene_Weak = which(CorSeq_abs < q_low)
Gene_Moderate = which( (CorSeq_abs < q_high) & (CorSeq_abs >= q_low) )
Gene_Strong = which(CorSeq_abs >= q_high)


