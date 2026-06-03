#!/usr/bin/env Rscript

# Load packages required for this script
write("Loading packages required for this script", stderr())
suppressPackageStartupMessages({
    library(DESeq2)
    library(optparse)
})

# Workaround for issue 112: https://github.com/thelovelab/DESeq2/issues/112
# This can probably be removed in the future
setOldClass("ExpData")

# Load dependency packages for testing installations
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

# Create parsing options list
option_list <- list(
    make_option(c("--rows"),
        type = "integer",
        help = "Number of rows in dummy matrix [default = %default]",
        default = 100
    )
)

# Implement parser with optparse
opt_parser <- OptionParser(
    option_list = option_list,
    description = "Runs DESeq2 on dummy data"
)

# Parse options with optparse
opt <- parse_args(opt_parser)

# Create a random dummy count matrix
cnts <- matrix(rnbinom(n = opt$row * 10, mu = 100, size = 1 / 0.5), ncol = 10)
cond <- factor(rep(1:2, each = 5))

# Object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~cond)

# Standard analysis
dds <- DESeq(dds)
res <- results(dds)

# Print results to stdout
print(res)
