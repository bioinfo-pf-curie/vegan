#!/usr/bin/env Rscript
.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

suppressMessages(library(facets))
library(optparse)
set.seed(1234)

option_list <- list(make_option(c("-i", "--input"), type="character", default=NULL, help="Facets pileup input file"),
                    make_option(c("-n", "--name"), type="character", default=NULL, help="Sample name"),
                    make_option(c("-o", "--outDir"), type="character", default='./', help="Output directory"),
                    make_option(c("-a", "--assembly"), type="character", default=NULL, help="Genome assembly, Must be hg18, hg19, hg38, mm9, mm10"),
                    make_option(c("--normalDepth"), type="numeric", default=25, help="Minimum normal sample depth"),
                    make_option(c("--maxDepth"), type="numeric", default=1000, help="Maximum coverage to consider a SNP"),
                    make_option(c("--windowSize"), type="numeric", default=250, help="Window size"),
                    make_option(c("--cval"), type="numeric", default=150, help="Critical value for segmentation"),
                    make_option(c("--cvalPreproc"), type="numeric", default=25, help="Criticial segmentation value for preprocessing"),
                    make_option(c("--ampCopy"), type="numeric", default=5, help="Copy number to call an amplification"),
                    #make_option(c("--maxCluster"), type="numeric", default=5, help=""),
                    make_option(c("--hetThres"), type="numeric", default=0.25, help="VAF value to call a SNP heterozygous"),
                    make_option(c("--unmatch"), action="store_true", default=FALSE, help="Normal sample is not matched with tumor. In this case, heterogygous SNPs are called using tumor reads only and logOR calculations are different. Use het.thresh = 0.1 or lower in that case.")
                    )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
name <- opt$name
outDir <- opt$outDir

#########################
##
## Functions
##
#########################

facetsPlot <- function (x, emfit = NULL, clustered = FALSE, chromlevels = c(1:22,"X"), plot.type = c("em", "naive", "both", "none"), sname = NULL) {
    plot.type <- match.arg(plot.type)
    if (plot.type == "none")
        layout(matrix(1:2, ncol = 1))
    if (plot.type == "em")
        layout(matrix(rep(1:4, c(9, 9, 6, 1)), ncol = 1))
    if (plot.type == "naive")
        layout(matrix(rep(1:4, c(9, 9, 6, 1)), ncol = 1))
    if (plot.type == "both")
        layout(matrix(rep(1:6, c(9, 9, 6, 1, 6, 1)), ncol = 1))
    par(mar = c(0.25, 3, 0.25, 1), mgp = c(1.75, 0.6, 0), oma = c(3,0, 1.25, 0))
    jseg <- x$jointseg
    if (missing(emfit)) {
        out <- x$out
        if (plot.type == "em" | plot.type == "both") {
            warning("emfit is missing; plot.type set to naive")
            plot.type <- "naive"
        }
    }
    else {
        out <- emfit$cncf
    }
    if (clustered) {
        cnlr.median <- out$cnlr.median.clust
        mafR <- out$mafR.clust
        mafR[is.na(mafR)] <- out$mafR[is.na(mafR)]
    }
    else {
        cnlr.median <- out$cnlr.median
        mafR <- out$mafR
    }
    mafR <- abs(mafR)
    chrcol <- 1 + rep(out$chr - 2 * floor(out$chr/2), out$num.mark)
    nn <- cumsum(table(jseg$chrom[is.finite(jseg$cnlr)]))
    segbdry <- cumsum(c(0, out$num.mark))
    segstart <- segbdry[-length(segbdry)]
    segend <- segbdry[-1]
    plot(jseg$cnlr[is.finite(jseg$cnlr)], pch = ".", cex = 2,  col = c("grey", "lightblue", "azure4", "slateblue")[chrcol],  ylab = "log-ratio", xaxt = "n")
    abline(h = median(jseg$cnlr, na.rm = TRUE), col = "green2")
    abline(h = x$dipLogR, col = "magenta4")
    segments(segstart, cnlr.median, segend, cnlr.median, lwd = 1.75, col = 2)
    abline(v=cumsum(tapply(out$num.mark,out$chr,sum)),col="grey80",lty=2)
    
    plot(jseg$valor[is.finite(jseg$cnlr)], pch = ".", cex = 2.5,  col = c("grey", "lightblue", "azure4", "slateblue")[chrcol], ylab = "log-odds-ratio", ylim = c(-4, 4), xaxt = "n")
    segments(segstart, sqrt(mafR), segend, sqrt(mafR), lwd = 1.75,  col = 2)
    segments(segstart, -sqrt(mafR), segend, -sqrt(mafR), lwd = 1.75,  col = 2)
    abline(v=cumsum(tapply(out$num.mark,out$chr,sum)),col="grey80",lty=2)
    cfpalette <- c(colorRampPalette(c("white", "steelblue"))(10), "bisque2")
    if (plot.type == "naive" | plot.type == "both") {
        out$tcn[out$tcn > 10] <- 9 + log10(out$tcn[out$tcn >10])
        ii <- which(out$lcn > 5)
        if (length(ii) > 0)
            out$lcn[ii] <- 5 + log10(out$lcn[ii])
        plot(c(0, length(jseg$cnlr)), c(0, max(out$tcn)), type = "n",  ylab = "copy number (nv)", xaxt = "n")
        segments(segstart, out$tcn, segend, out$tcn, lwd = 1.75, col = 1)
        segments(segstart, out$lcn, segend, out$lcn, lwd = 1.75, col = 2)
        abline(v=cumsum(tapply(out$num.mark,out$chr,sum)),col="grey80",lty=2)
        plot(c(0, length(jseg$cnlr)), 0:1, type = "n", ylab = "",  xaxt = "n", yaxt = "n")
        mtext("cf-nv", side = 2, at = 0.5, line = 0.3, las = 2,  cex = 0.75)
        cfcol <- cfpalette[round(10 * out$cf + 0.501)]
        rect(segstart, 0, segend, 1, col = cfcol, border = NA)
        abline(v=cumsum(tapply(out$num.mark,out$chr,sum)),col="grey65",lty=2)
    }
    if (plot.type == "em" | plot.type == "both") {
        out$tcn.em[out$tcn.em > 10] <- 9 + log10(out$tcn.em[out$tcn.em > 10])
        ii <- which(out$lcn.em > 5)
        if (length(ii) > 0)
            out$lcn.em[ii] <- 5 + log10(out$lcn.em[ii])
        plot(c(0, length(jseg$cnlr)), c(0, max(out$tcn.em)), type = "n", ylab = "copy number (em)", xaxt = "n")
        segments(segstart, out$tcn.em, segend, out$tcn.em, lwd = 1.75, col = 1)
        segments(segstart, out$lcn.em, segend, out$lcn.em, lwd = 1.75,col = 2)
        abline(v=cumsum(tapply(out$num.mark,out$chr,sum)),col="grey80",lty=2)
        plot(c(0, length(jseg$cnlr)), 0:1, type = "n", ylab = "", xaxt = "n", yaxt = "n")
        mtext("cf-em", side = 2, at = 0.5, line = 0.2, las = 2,cex = 0.75)
        cfcol <- cfpalette[round(10 * out$cf.em + 0.501)]
        rect(segstart, 0, segend, 1, col = cfcol, border = NA)
        abline(v=cumsum(tapply(out$num.mark,out$chr,sum)),col="grey50",lty=2)
        
    }
    if (missing(chromlevels)) {
        if (length(nn) == 23) {
            axis(labels = c(1:22, "X"), side = 1, at = (nn + c(0, nn[-23]))/2, cex = 0.65)
        }
        else {
            axis(labels = 1:22, side = 1, at = (nn + c(0, nn[-22]))/2,cex = 0.65)
        }
    }
    else {
        axis(labels = chromlevels, side = 1, at = (nn + c(0, nn[-length(nn)]))/2, cex = 0.65)
    }
    mtext(side = 1, line = 1.75, "Chromosome", cex = 0.8)
    pos= (nn + c(0, nn[-length(nn)]))/2
    cellularfraction=round(out$cf[which(duplicated(out$cf)==F)],2)
    cellularfractioncol=cfcol[which(duplicated(out$cf)==F)]
    mtext(side=1,line=1.75,at=pos[1],col="black",cex=0.8,font=2,text="cell frac:")
    lc=lapply(1:length(cellularfraction),function(i){mtext(side = 1, line = 1.75, at=pos[i+1], cellularfraction[i], cex = 0.8,col=cellularfractioncol[i],font=2)})
    
    if (!missing(sname))
        mtext(sname, side = 3, line = 0, outer = TRUE, cex = 0.8)
}


#####################
##
## Facets
##
#####################

## Load pileup
message("Loading pileup ...")
rcmat <- readSnpMatrix(opt$input)
rcmat

xx <- preProcSample(rcmat, ndepth=opt$normalDepth, het.thresh=opt$hetThres,
                    snp.nbhd=opt$windowSize, cval=opt$cvalPreproc, deltaCN=0,
                    gbuild=opt$assembly, hetscale=TRUE, unmatched=opt$unmatch,
                    ndepthmax=opt$maxDepth)

pc_het_snp <- round(nrow(xx$pmat[which(xx$pmat$het==1),])*100/nrow(xx$pmat),2)

oo <- suppressWarnings(procSample(xx, cval=opt$cval,
                                  min.nhet=15, dipLogR=NULL))

diplogR <- oo$dipLogR
flags <- oo$flags

## Estimate copy number and cellular fraction
fit2 <- suppressWarnings(emcncf(oo, trace=FALSE, unif=FALSE,
                                 min.nhet=15, maxiter=50, eps=1e-3))

fit2$cncf$purity <- fit2$purity
fit2$cncf$ploidy <- fit2$ploidy

purity <- round(unique(fit2$cncf$purity),2)
ploidy <- round(unique(fit2$cncf$ploidy),2)
    
if(!is.null(fit2$emflags)){
    title <- paste0("Allele-specific facets CNV profile of ",name, " : cellularity=",purity, ", ploidy=",ploidy, " , Warning:", fit2$emflags)
} else {
    title <- paste0("Allele-specific facets CNV profile of ",name, " : cellularity=", purity, ", ploidy=", ploidy)
}

write.table(fit2$cncf,paste0(outDir,"/",name,"_subclonal_allele_spe_cnv_", purity, "cellularity_", ploidy,"ploidy.txt"),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

pdf(paste0(outDir,"/",name,"_subclonal_allele_spe_cnv_", purity, "cellularity_", ploidy, "ploidy.pdf"), width=16, height=8)
facetsPlot(x=oo, emfit=fit2, sname=title, chromlevels=unique(xx$pmat$chrom))
invisible(dev.off())

png(paste(outDir,"/",name,"_subclonal_allele_spe_cnv.png",sep=""), width=1200, height=800)
facetsPlot(x=oo, emfit=fit2, sname=title, chromlevels=unique(xx$pmat$chrom))
invisible(dev.off())

## Summarize results in a table
tab <- fit2$cncf
tab$ID <- name
tab$ploidy <- unique(round(tab$ploidy))
tab[is.na(tab$lcn.em),"lcn.em"] <- 0
tab$A <- tab$tcn.em-tab$lcn.em
tab$B <- tab$lcn.em
tab$Ag <- apply(tab[,c("A","B")],1,function(x){if(x[1]!=0){paste(rep("A",x[1]),collapse="")}})
tab$Bg <- apply(tab[,c("A","B")],1,function(x){if(x[2]!=0){paste(rep("B",x[2]),collapse="")}})
tab$Geno <- paste(tab$Ag,tab$Bg,sep="")
tab$Geno <- gsub("NULL","",tab$Geno)
tab <- tab[,c("ID","chrom","start","end","tcn.em","Geno","cnlr.median.clust","ploidy")]
colnames(tab) <- c("ID","chrom","loc.start","loc.end","CNt","Geno","logratio","ploidy")

tab$call="NEUTRAL"
if(length(which(tab$CNt>floor(unique(tab$ploidy))))>0)
    tab[which(tab$CNt>floor(unique(tab$ploidy))),"call"]="GAIN"
if(length(which(tab$CNt>(opt$ampCopy+unique(tab$ploidy))))>0)
    tab[which(tab$CNt>(opt$ampCopy+unique(tab$ploidy))),"call"]="AMP"
if(length(which(tab$CNt<floor(unique(tab$ploidy))))>0)
    tab[which(tab$CNt<floor(unique(tab$ploidy))),"call"]="LOSS"
if(length(which(tab$CNt==0))>0)
    tab[which(tab$CNt==0),"call"]="DEL"
tab$LOH <- 0
tab[grep("B",tab$Geno,invert=T),"LOH"] <- 1
tab$chrom <- as.character(tab$chrom)

if(opt$assembly=="mm10" | opt$assembly=="mm9")
    tab[which(tab$chrom=="20"),"chrom"]="X"
if(opt$assembly=="hg19" | opt$assembly=="hg18")
    tab[which(tab$chrom=="23"),"chrom"]="X"

write.table(tab, paste0(outDir,"/",name,"_subclonal_allele_spe_cnv_",round(unique(fit2$cncf$purity),2),"cellularity_",unique(tab$ploidy),"ploidy.transformed.txt"),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

write.table(fit2$emflags, paste0(outDir,"/",name,"_subclonal_allele_spe_cnv_",round(unique(fit2$cncf$purity),2),"cellularity_",unique(tab$ploidy),"ploidy.flags.txt"),
            quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

pp<-unique(data.frame(name=name,ploidy=fit2$cncf$ploidy, purity=fit2$cncf$purity))
write.table(pp, paste0(outDir,"/",name,"_subclonal_allele_spe_cnv_",round(unique(fit2$cncf$purity),2),"cellularity_",unique(tab$ploidy),"ploidy.optimal.txt",sep=""),
            quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")


