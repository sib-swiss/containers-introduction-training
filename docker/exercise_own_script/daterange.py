#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = "Get a daterange")

parser.add_argument('-d', '--date', type=str, required=True, 
                    help='Date. Format: [YYYYMMDD]')

args = parser.parse_args()

dates = pd.date_range(args.date, periods=7)

for d in dates:
    print(d)
