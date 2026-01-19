#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: apTransition.r <inputTable> <outputTable>", call.=FALSE)
}

#load library
library(dplyr)

# Load arguments
inputTable <- args[1]
outputTable <- args[2]

## Load tmb scores
transi=read.table(inputTable, header=TRUE, sep= "\t")

#Select only SNV
transi_snv=transi[which(transi$type =="snp"),]

#Add column with squashed transitions
transi_snv$transi=paste(transi_snv$ref,transi_snv$alt, sep = "_")

#Compute number of transition & indels

table_transi=table(factor(transi_snv$transi,levels = c("A_C","A_G","A_T","C_A","C_G","C_T","G_A","G_C","G_T","T_A","T_C","T_G")))
indels=length(which(transi$type == "indel"))

table_transi=data.frame(cbind(ID=as.character(transi_snv$ID[1]),t(as.matrix(table_transi)),indels))

# Squash transitions:
table_transi=table_transi%>%
  mutate(AG_TC=as.numeric(as.character(table_transi$A_G))+as.numeric(as.character(table_transi$T_C)))%>%
  mutate(CT_GA=as.numeric(as.character(table_transi$C_T))+as.numeric(as.character(table_transi$G_A)))%>%
  mutate(AC_TG=as.numeric(as.character(table_transi$A_C))+as.numeric(as.character(table_transi$T_G)))%>%
  mutate(CA_GT=as.numeric(as.character(table_transi$C_A))+as.numeric(as.character(table_transi$G_T)))%>%
  mutate(AT_TA=as.numeric(as.character(table_transi$A_T))+as.numeric(as.character(table_transi$T_A)))%>%
  mutate(CG_GC=as.numeric(as.character(table_transi$C_G))+as.numeric(as.character(table_transi$G_C)))

table_transi=table_transi[,c("ID","AG_TC","CT_GA","AC_TG","CA_GT","AT_TA","CG_GC","indels")]
write.table(table_transi,outputTable, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, )
