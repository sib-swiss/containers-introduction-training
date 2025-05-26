#!/usr/bin/env Rscript

# load packages required for this script
write("Loading packages required for this script", stderr())
suppressPackageStartupMessages({
    library(DESeq2)
    library(optparse)
})

# workaround for issue 112: https://github.com/thelovelab/DESeq2/issues/112
# this can probably be removed in the future
setOldClass("ExpData")

# load dependency packages for testing installations
write("Loading dependency packages for testing installations", stderr())
suppressPackageStartupMessages({
    library(apeglm)
    library(IHW)
    library(limma)
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(RColorBrewer)
    library(scales)
    library(stringr)
})

# parse options with optparse
option_list <- list(
    make_option(c("--rows"),
        type = "integer",
        help = "Number of rows in dummy matrix [default = %default]",
        default = 100
    )
)

opt_parser <- OptionParser(
    option_list = option_list,
    description = "Runs DESeq2 on dummy data"
)
opt <- parse_args(opt_parser)

# create a random dummy count matrix
cnts <- matrix(rnbinom(n = opt$row * 10, mu = 100, size = 1 / 0.5), ncol = 10)
cond <- factor(rep(1:2, each = 5))

# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~cond)

# standard analysis
dds <- DESeq(dds)
res <- results(dds)

# print results to stdout
print(res)
