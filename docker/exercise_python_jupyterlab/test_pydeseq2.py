#!/usr/bin/env python3

# Import libraries
import sys

# Load packages required for this script
print('Loading packages required for this script', file=sys.stderr)

# ruff: disable[E402]
import argparse
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
# ruff: enable[E402]

# Create argument parser
parser = argparse.ArgumentParser(description='Runs PyDESeq2 on dummy data')
parser.add_argument(
    '--rows',
    type=int,
    default=100,
    help='Number of rows (genes) in dummy matrix [default: %(default)s]',
)

# Parse arguments
print('Parsing arguments', file=sys.stderr)
args = parser.parse_args()

# Create a random dummy count matrix (samples x genes)
print('Creating dummy count matrix', file=sys.stderr)
counts_df = pd.DataFrame(
    np.random.negative_binomial(n=2, p=0.02, size=(10, args.rows)),
    index=[f'sample{i + 1}' for i in range(10)],
    columns=[f'gene{i + 1}' for i in range(args.rows)],
)

# Metadata / condition factor
metadata = pd.DataFrame(
    {'condition': ['A'] * 5 + ['B'] * 5}, index=counts_df.index
)

# Run DESeq2 analysis
print('Running PyDESeq2', file=sys.stderr)
dds = DeseqDataSet(
    counts=counts_df, metadata=metadata, design_factors='condition'
)
dds.deseq2()

# Get results
print('Computing results...', file=sys.stderr)
stat_res = DeseqStats(dds, contrast=["condition", "B", "A"])
stat_res.summary()

# Print all results (commented out by default to avoid flooding stdout)
# print('Results:', file=sys.stderr)
# print(stat_res.results_df.to_string())
