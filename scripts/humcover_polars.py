#!/usr/bin/env python3

import sys
import os
import argparse
import pickle
import pandas as pd
import polars as pl
from pathlib import Path
import gffutils
from multiprocessing import Pool

from Bio import SeqIO


aparser = argparse.ArgumentParser(description="Filter to keep hits where most of human protein is covered.")
aparser.add_argument('-d', "--dataframe", help="Pickled data frame containing filtered search hits")
aparser.add_argument('-f', "--fasta_humprot", help="FASTA file containing human protein sequences")
aparser.add_argument('-o', "--outfile", help="Output file receiving data frame with augmented search hits")
aparser.add_argument('-c', "--coverage_full",
                     help=("Select full coverage of human sequences. Defaults to partial coverage using internal"
                     " threshold."),
                     action="store_true")
aparser.add_argument('-I', "--IPC",
                     help=("Dataset is in IPC mode, not pickled"),
                     action="store_true")
aparser.add_argument('-G', "--group_mode",
                     help=("Calculate hum_cover over groups"),
                     action="store_true")
aparser.add_argument('-p', "--partpercent",
                     type=int,
                     help="Coverage fraction of partial human protein matches")
aparser.add_argument('-C', "--covdiff",
                     type=int,
                     default=20,
                     help="Percent humcover difference necessary to be counted")
aparser.add_argument('-m', "--minbaclength",
                     type=int,
                     default=80,
                     help="Minimum length of a bacterial protein to be counted")
args = aparser.parse_args()

mincount = 1
if args.coverage_full:
    fullcover = True
else:
    fullcover = False
    mincount = 2
if args.partpercent:
    cov_fraction = args.partpercent / 100
    print(f"cov_fraction = {cov_fraction}")
else:
    print("A --partpercent value must be provided.")
    sys.exit(1)

if args.IPC:
    results = pl.read_ipc(args.dataframe)
else:
    results = pl.from_pandas(pd.read_pickle(args.dataframe))

covdiff = args.covdiff / 100

# Extract the genome type from the data frame columns.  Determined by samegenome step.
dfcols = results.columns
if "src_genome" in dfcols:
    print("Source genomes found")
    genometype = "src_genome"
elif "rep_genome" in dfcols:
    print("Representative genomes found")
    genometype = "rep_genome"
elif "g" in dfcols:
    print("Generic genomes found")
    genometype = "g"
elif "Lineage" in dfcols:
    print("Using lineage as substitute for genomes [must be already merged]")
    genometype = "Lineage"
else:
    print("No genome type column found in data frame.")
    sys.exit(1)


# Figure out which columns to retain

result_selection = ["hum_protein", genometype, "sstart", "send"]
if args.group_mode:
    result_selection.append("groupnum")

# For these purposes, we only care about the intervals covering human proteins,
# so we can collapse down all the bacterial proteins for now. Also take this
# time to get rid of very small bacterial proteins and cases where there are too
# few remaining counts

print("Filtering results...")
resultsF = results.filter(
    pl.col("length_bac") >= args.minbaclength
).with_columns(
    count = pl.col("bac_protein").count().over(["hum_protein", genometype])
).filter(
    pl.col("count") >= mincount
).select(
    result_selection
).unique()

# Sort by starting and ending location within each hum_protein/genome pair.
# Then, shift the ending location down one so we can find the running maximum.
# If an sstart is ahead of the last running maximum, it's a new interval.
# Getting the cumulative sum over the booleans of "new interval" gives you an
# interval group per hum_protein/genome pair.
print("Computing intervals...")
intervals = resultsF.sort("hum_protein", genometype, "sstart","send").with_columns(
    last_send=pl.col("send").shift(fill_value=0).over(["hum_protein", genometype])
).with_columns(
    running_max = pl.max_horizontal("send","last_send")
).with_columns(
    last_running_max=pl.col("running_max").shift(fill_value=0).over(["hum_protein", genometype])
).with_columns(
    interval_group = (pl.col("sstart") > pl.col("running_max").shift(fill_value=0)).cum_sum().over(["hum_protein", genometype])
)

# Now summarize first each interval group (collapsing to min and max), then each
# distinct interval (since we now know these don't overlap) to get the total length
print("Tallying intervals...")

# If this is the second pass, we've organized our results into groups
if args.group_mode:
    totlength_groups = ["hum_protein", genometype, "groupnum"]
else:
    totlength_groups = ["hum_protein", genometype]

total_length = intervals.group_by(totlength_groups + ["interval_group"]).agg(
    length=pl.col("send").max() - pl.col("sstart").min()
).group_by(totlength_groups).agg(
    tot_length=pl.col("length").sum()
)


# Read in FASTA file to get human protein lengths
print("Reading in FASTA...")
ids = []
lengths = []
hprot_seqs = args.fasta_humprot
for rec in SeqIO.parse(hprot_seqs, "fasta"):
    ids.append(rec.id)
    lengths.append(len(rec.seq))
hseqs = pl.DataFrame({"hum_protein": ids, "length_hum": lengths}, schema={"hum_protein": pl.Utf8, "length_hum": pl.UInt64})

# Join together with results
print("Computing totals...")
total_humcover = total_length.join(hseqs, on="hum_protein", how="left").with_columns(
    tot_humcover = pl.col("tot_length") / pl.col("length_hum")
)

# Now get "best" individual humcover

print("Computing individual totals...")
indiv_humcover = resultsF.with_columns(
    indiv_length = pl.col("send") - pl.col("sstart")
).join(
    hseqs, on="hum_protein", how="left"
).with_columns(
    indiv_humcover = pl.col("indiv_length") / pl.col("length_hum")
).group_by(["hum_protein", genometype]).agg(
    best_indiv_humcover = pl.col("indiv_humcover").max()
)

overall = indiv_humcover.join(total_humcover, on=["hum_protein", genometype], how="inner")

print("Getting final results...")
if fullcover:
    overallF = overall.filter(
        pl.col("best_indiv_humcover") >= cov_fraction
    )
else:
    # Merge two analyses, then filter for cases where total humcover is over
    # covdiff% more than best_indiv_humcover, humcover is at least cov_fraction, and
    # best_indiv_humcover is less than cov_fraction
    overallF = overall.filter(
        pl.col("best_indiv_humcover") < cov_fraction,
        pl.col("tot_humcover") >= cov_fraction,
        (pl.col("tot_humcover") - covdiff) > pl.col("best_indiv_humcover")
    )


# Prune original results to only keep entries in overallF
print("Writing results...")
if args.group_mode:
    final_results = results.join(overallF, on=["hum_protein", genometype, "groupnum"], how="inner")
else:
    final_results = results.join(overallF, on=["hum_protein", genometype], how="inner")
final_results.write_ipc(args.outfile)
