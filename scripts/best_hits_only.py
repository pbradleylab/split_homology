#!/usr/bin/env python3


#!/usr/bin/env python3

import os
import sys
import argparse
import polars as pl

aparser = argparse.ArgumentParser(description=("Keep only the best hit per bacterial protein"))
aparser.add_argument("-d", "--dataframe", help="Input IPC dataframe containing hits data")
aparser.add_argument('-o', "--out", help="IPC dataframe file receiving location values")

args = aparser.parse_args()

if args.dataframe.endswith(".ipc"):
    hits = pl.scan_ipc(args.dataframe)
elif args.dataframe.endswith(".tsv"):
    hits = pl.scan_csv(args.dataframe, separator="\t", infer_schema_length=13100999)
else:
    raise ValueError

# Sort in descending order of bitscore and then take the top human protein hit by bacterial protein/lineage
consolidated = hits.sort(
    ["bac_protein", "Lineage", "bitscore"],
    descending = [False, False, True]
).group_by(
    ["bac_protein", "Lineage"], maintain_order=True
).first()

consolidated.collect(streaming=True).write_ipc(args.out)
