#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: compute_clust.r <inputTable> <outputDir> <figureName>", call.=FALSE)
}

# Load arguments
inputTable <- args[1]
outputDir<-args[2]
figureName<-args[3]

#inputTable <- "/bioinfo/users/tgutman/Documents/Tom/SeqOIA/SeqOIA_KDI/analysis/test_clustering/bam/clust_matrix.tsv"
#outputDir<-"/data/users/tgutman/Documents/Tom/SeqOIA/SeqOIA_KDI/analysis/test_clustering/bam"
#figureName<-"test_clust"

# Handle path & Output Names
outputFile=paste(outputDir,"/",figureName,sep = "")
outputPlot=paste(outputDir,"/",figureName,".png",sep = "")
outputMatrix=paste(outputDir,"/",figureName,"_identito.csv",sep = "")

# Load libraries
library(vegan)
library(pheatmap)
library(RColorBrewer)

# Saving function
save_pheatmap_png <- function(x, filename, width=1000, height=1000) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    png(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

#Load data
clust_mat=read.table(inputTable, header=T, stringsAsFactors = F)
#clust_mat[,-1]=apply(clust_mat[,-1],2,function(x){x[which(x=="NEC")]=0;return(as.numeric(x))})
rownames(clust_mat)=clust_mat$rs
clust_mat=clust_mat[,-1]

temp_colnames=colnames(clust_mat)
temp_rownames=rownames(clust_mat)

# Filter matrix
#small_clust_mat=clust_mat[, !sapply(clust_mat, is.character)]

# Handle NA:
clust_mat=as.matrix(clust_mat)
clust_mat[which(clust_mat == "NEC")]=NA
clust_mat=matrix(as.numeric(clust_mat),ncol=ncol(clust_mat))

colnames(clust_mat)=temp_colnames
rownames(clust_mat)=temp_rownames

# Handle sample full of NA
for (i in 1:nrow(clust_mat)){
    if (sum(is.na(clust_mat[i,])) == ncol(clust_mat)){
        warning('One sample contains only NAs, this sample will be removed')
        clust_mat=clust_mat[-i,]
        print(clust_mat)
    }
}

# Compute distance
jacdist <- vegdist(clust_mat,method="jaccard",na.rm=TRUE)
jacdist_mat = as.matrix(jacdist)

## reduce matrix
#jacdist <- vegdist(small_clust_mat,method="jaccard",na.rm=TRUE)
#jacdist_mat = as.matrix(jacdist)

## 1-dist matrix:
corr_mat= 1-jacdist_mat

# Create heatmap
clust=pheatmap(corr_mat,
         method = "Ward.D",
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         display_numbers = TRUE,
         legend = TRUE ,
         fontsize = 20,
         silent = TRUE)

 clust_order=clust$tree_col$order

 # Save figure:
 save_pheatmap_png(clust, outputPlot)
 write.csv(corr_mat[clust_order,clust_order], outputMatrix, quote=TRUE, row.names=TRUE)

# End message + Path:
print("script ended successfully")
print(paste("plot can be found here:", outputFile))
