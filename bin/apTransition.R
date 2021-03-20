#!/usr/bin/env Rscript

## Load tmb scores
transi=read.table("transi.tsv", header=T, sep= "\t")

transi_snv=transi[which(transi$type =="snp"),]
# Add column with squashed transitions
transi_snv$transi=paste(transi_snv$ref,transi_snv$alt, sep = "_")

test_table = NULL
indels = NULL

for(Barcode in unique(transi_snv$ID)){
    table_transi=table(transi_snv[which(transi_snv$ID == Barcode),"transi"])
    print(table_transi)
    indels=c(indels,length(which(transi$type == "indel" & transi$ID == Barcode)))
    test_table = rbind(test_table,table_transi)
}

test_table=cbind(test_table,indels)
rownames(test_table) = unique(transi_snv$ID)

write.csv(test_table,"transiTable.tsv", row.names = TRUE, col.names = TRUE)