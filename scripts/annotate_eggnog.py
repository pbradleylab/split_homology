#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import polars as pl


aparser = argparse.ArgumentParser(description="Assign EGGnog annotations to UHGP IDs in supplied search hits.")
aparser.add_argument('-d', "--dataframe", help="Pickled or IPC data frame containing filtered search hits.")
aparser.add_argument('-e', "--eggnog", help="Eggnog annotations TSV file.")
aparser.add_argument('-o', "--out", help="Output file receiving pickled or IPC data frame with augmented search hits.")

args = aparser.parse_args()

print("Reading dataframe...")
if args.dataframe.endswith(".ipc"):
    hits = pl.scan_ipc(args.dataframe)
elif args.dataframe.endswith(".pkl"):
    hits = pl.from_pandas(pd.read_pickle(args.dataframe))
else:
    raise ValueError

print("Reading and extracing EGGnog data...")
#eggnog = pd.read_csv(args.eggnog, sep='\t', names=range(0,22), header=None)
#eggnog = eggnog[[0,6,8,9,10,11,12,13,21]]
#eggnog.columns = ["bac_protein", "EN1", "EN2", "EN3", "EN4", "EN5", "EN6", "EN7", "desc"]
#eggnog = pl.from_pandas(eggnog).lazy()

eggsel = [f"column_{i+1}" for i in [0, 6, 8, 9, 10, 11, 12, 13, 21]]
eggreplace = ["bac_protein"] + [f"EN{i}" for i in range(1,8)] + ["desc"]

eggnog = pl.scan_csv(
    args.eggnog, separator='\t', has_header=False
).select(
    eggsel
).rename(
    dict(zip(eggsel, eggreplace))
)


print("merging")
hits2 = hits.join(eggnog, how="left", on="bac_protein")

if args.out.endswith(".pkl"):
    hits2.collect().to_pandas().to_pickle(args.out)
elif args.out.endswith(".ipc"):
    hits2.collect(streaming=True).write_ipc(args.out)
elif args.out.endswith(".tsv"):
    hits2.collect(streaming=True).write_csv(args.out, separator='\t')
else:
    raise ValueError
