#!/usr/bin/env Rscript
.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

library(data.table)
library(optparse)
library(ggplot2)

option_list = list(
    make_option("--cov", type="character", default=NULL, help="Gene coverage file (bed)", metavar="character"),
    make_option("--oprefix", type="character", default=NULL, help="Output file", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$cov)){
    print_help(opt_parser)
    stop("Please provide a bed file with the gene coverage.", call.=FALSE)
}
if (is.null(opt$output)){
    opt$output = "covdensity.tsv"

}

x <- read.table(opt$cov, sep="\t")
dt <- data.table(x)
gene.cov <- dt[, .(mean=mean(V5)), by=list(V4)]
if (max(gene.cov$mean>1500)){
    dens <- hist(gene.cov$mean, breaks=c(seq(0,100,by=10),seq(150,1500,by=50),ceiling(max(gene.cov$mean))))
}else{
    dens <- hist(gene.cov$mean, breaks=c(seq(0,100,by=10),seq(150,1500,by=50)))
}
    
dens.df <- data.frame(mean=dens$breaks, counts=c(1,1-cumsum(dens$counts/sum(dens$counts)))*100)


## plot
pdf(paste0(opt$oprefix, ".pdf"))
ggplot(dens.df) + geom_line(aes(y=counts,x=mean)) + theme_classic() + ylim("Cumumative Gene Fraction (%)") + xlim("Coverage (X)")
dev.off()

## save
write.table(dens.df, file=paste0(opt$oprefix, ".mqc"), sep="\t", row.names=FALSE, col.names=FALSE)
