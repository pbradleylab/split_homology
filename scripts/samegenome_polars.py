#!/usr/bin/env python3
import math
import psutil
import sys
import os
import argparse
import pickle
from datetime import datetime
import pandas as pd
import multiprocessing
from multiprocessing import Pool
import shelve
#import pyarrow as pa
#import pyarrow.dataset as ds
import polars as pl

from collections import namedtuple


process = psutil.Process()

aparser = argparse.ArgumentParser(description=("Select hits for proteins seen "
        "together in same genome."))
aparser.add_argument('-d', "--dataframe",
        help="Pickled data frame containing filtered search hits.")
aparser.add_argument('-a', "--arrowmap",
        help="cluster_map UHGP arrow dataset")
aparser.add_argument('-G', "--genometype",
        required=True,
        help="Either 'src' for source "
        "(GUT_GENOME...) genomes or 'rep' for representative (MGYG...) genomes")
aparser.add_argument('-T', "--test",
                     default=False,
                     action='store_true',
        help="Do an abridged version for testing")
aparser.add_argument('-t', "--hits_threshold",
        required=False,
        type=int,
        help="Threshold of number of blast hits in a given human protein group above which "
        "the human protein is excluded from the results. If not supplied, all hits are "
        "used in this script (may impact processing time significantly).")
aparser.add_argument('-m', "--minimum_per_genome",
        default=1,
        type=int,
        help="Minimum number of hits per genome to retain")
aparser.add_argument('-o', "--outfile",
        help="Output file (IPC format)")
args = aparser.parse_args()

# Create output directory if it does not exist

results = pd.read_pickle(args.dataframe)

# add clusterID (by removing text from bac_protein)
resultsP = pl.from_pandas(results).with_columns(clusterID = pl.col("bac_protein").str.replace("GUT_GENOME", "").str.replace("_","").cast(pl.UInt64) )

# if testing, make this finish faster...
if args.test:
    this_n_rows = 250000
else:
    this_n_rows = None

# syntax to read partitioned dataset is a little different (*/*);
# join to the results as we read in
dset = pl.read_ipc(os.path.join(args.arrowmap, "clusterG*/*"),
                   n_rows=this_n_rows).with_columns(
    g = (pl.col("proteinID")/100000).floor().cast(pl.UInt32())
).join(
    resultsP, on="clusterID", how="inner"
)

print(f" -- Dataset read in...")
print(f" -- Memory usage: {process.memory_info().rss}")

print(f" -- Filtering and writing to {args.outfile}...")
# Now filter only cases with at least one bac_protein hit per genome and write to disk
# (Note, we will later filter for cases where there are multiple hits)
dset.select(
    pl.all()
    .sort_by("bitscore", descending=True)
    .over(["hum_protein", "g"])
).with_columns(
    count = pl.col("bac_protein")
    .count()
    .over(["hum_protein", "g"])
).filter(
    pl.col("count") >= args.minimum_per_genome
).write_ipc(args.outfile)

print(f" -- Finished!")
print(f" -- Memory usage: {process.memory_info().rss}")
sys.stdout.flush()
