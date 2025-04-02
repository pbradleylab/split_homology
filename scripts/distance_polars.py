#!/usr/bin/env python3

# Polars rewrite of distance using shifts (much faster)

import os
import sys
import argparse
import pandas as pd
import polars as pl

LargeInt = 999999

aparser = argparse.ArgumentParser(description=(" "))
aparser.add_argument("-d", "--dataframe", help="Pickled or IPC dataframe containing hits data")
aparser.add_argument('-o', "--out", help="IPC dataframe file receiving location values")
aparser.add_argument('-m', "--maxdist", help="Maximum distance to retain (default: 10000)", default=10000, type=int)
args = aparser.parse_args()

if args.dataframe.endswith(".ipc"):
    hits = pl.read_ipc(args.dataframe)
elif args.dataframe.endswith(".pkl"):
    hits = pl.from_pandas(pd.read_pickle(args.dataframe))
else:
    raise ValueError

# Extract the genome type from the data frame columns.  Determined by samgenome step.
dfcols = hits.columns
if "src_genome" in dfcols:
    print("Source genomes found")
    genometype = "src_genome"
elif "rep_genome" in dfcols:
    print("Representative genomes found")
    genometype = "rep_genome"
elif "g" in dfcols:
    print("Generic genomes found")
    genometype = "g"
else:
    print("No genome type column found in data frame.")
    sys.exit(1)

# remove redundant entries (keep best bitscore)
hits = hits.sort(
    ["hum_protein",genometype,"bac_protein","contig","strand","featnum","proteinID","bitscore"],
    descending=[False,False,False,False,False,False,False,True]
)
hits = hits.group_by(
    ["hum_protein",genometype,"bac_protein","contig","strand","featnum"]
).first(
).sort(
    ["hum_protein",genometype,"bac_protein","contig","strand","featnum","proteinID","bitscore"],
    descending=[False,False,False,False,False,False,False,True]
)

filtered_hits = hits.with_columns(
    lastfn = pl.col("featnum").shift().over(["hum_protein", genometype, "contig", "strand"]),
    nextfn = pl.col("featnum").shift(-1).over(["hum_protein", genometype, "contig", "strand"])
).with_columns(
    dist_next = abs(pl.col("featnum") - pl.col("lastfn")),
    dist_last = abs(pl.col("featnum") - pl.col("nextfn"))
).with_columns(
    dist_nearest = pl.when(
        pl.col("dist_next").is_null()
    ).then(
        pl.col("dist_last")
    ).when(
        pl.col("dist_last").is_null()
    ).then(
        pl.col("dist_next")
    ).otherwise(
        pl.min_horizontal("dist_next", "dist_last")
    )
).filter(
    ~(pl.col("dist_nearest").is_null())
).filter(
    pl.col("dist_nearest") <= args.maxdist
).select(
    pl.exclude(["lastfn", "nextfn", "dist_next", "dist_last"])
).sort(
    ["hum_protein", genometype, "contig", "strand", "featnum"]
)

# Add group numbers. If the feature number changes too much within a
# hum_protein+genome+contig+strand group, then it's part of a separate putative
# operon, so we increment.

print("Adding group numbers...")
filtered_hits = filtered_hits.with_columns(
    lastfn = pl.col("featnum").shift(1, fill_value=LargeInt).over(
        ["hum_protein", genometype, "contig", "strand"]
    )
).with_columns(
    skip = abs(pl.col("featnum") - pl.col("lastfn")) > args.maxdist
).with_columns(
    groupnum = pl.col("skip").cum_sum()
).select(
    pl.exclude(["lastfn", "skip"])
)

print("Writing to disk...")
filtered_hits.write_ipc(args.out)
