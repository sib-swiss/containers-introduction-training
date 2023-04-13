#!/usr/bin/env Rscript

library(biomaRt)
library(optparse)

option_list <- list(
    make_option(c("--pattern"),
        type = "character",
        help = "Search pattern [default = %default]",
        default = "mouse"
    )
)

opt_parser <- OptionParser(
    option_list = option_list,
    description = "Searches biomaRt ensembl datasets"
)
opt <- parse_args(opt_parser)

ensembl <- useEnsembl(biomart = "ensembl")
searchDatasets(mart = ensembl, pattern = opt$pattern)
