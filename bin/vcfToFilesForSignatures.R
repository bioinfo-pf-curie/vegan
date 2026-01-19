#!/usr/bin/env Rscript

# rm(list=ls())

# if (!suppressPackageStartupMessages(require("tidyverse"))) {install.packages("tidyverse", dependencies=TRUE, repos='http://cloud.r-project.org/')}
# suppressPackageStartupMessages(require(tidyverse))
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE) 
indir <- args[1] ; stopifnot(!is.na(indir)) 
outdir <- args[2] ; stopifnot(!is.na(outdir)) 
pairid <- args[3] ; stopifnot(!is.na(pairid))

file <- dir(path = indir, pattern = "*.tsv", recursive=TRUE,full.names=TRUE)

res <- read_delim(file, col_names = TRUE, delim = "\t",show_col_types = FALSE) %>% 
      mutate(ID = pairid) %>%
      mutate(chr = gsub('chr','',CHROM)) %>%
      rename(pos=POS, ref=REF, alt=ALT) %>%
      select(chr, pos, ID, ref, alt)

write.table(res, paste(outdir,"/",pairid,"_sigprofiler.vcf",sep=""), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")