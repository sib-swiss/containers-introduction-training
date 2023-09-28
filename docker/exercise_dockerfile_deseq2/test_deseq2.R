#!/usr/bin/env Rscript

library(DESeq2)
library(optparse)

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



# see vignette for suggestions on generating
# count tables from RNA-Seq data
cnts <- matrix(rnbinom(n = opt$row * 10, mu = 100, size = 1 / 0.5), ncol = 10)
cond <- factor(rep(1:2, each = 5))

# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~cond)

# standard analysis
dds <- DESeq(dds)
res <- results(dds)

print(res)
