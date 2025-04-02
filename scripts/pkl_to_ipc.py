#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import polars as pl

aparser = argparse.ArgumentParser(description=(" "))
aparser.add_argument("-i", "--inputf", help="Input file")
aparser.add_argument('-o', "--outf", help="Output file")
args = aparser.parse_args()

if args.inputf.endswith(".pkl"):
    inputdata = pl.from_pandas(pd.read_pickle(args.inputf))
    inputdata.write_ipc(args.outf)
elif args.inputf.endswith(".ipc"):
    inputdata = pl.read_ipc(args.inputf).to_pandas()
    inputdata.to_pickle(args.outf)
else:
    raise ValueError
