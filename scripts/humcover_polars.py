#!/usr/bin/env python3

import sys
import os
import shutil
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
aparser.add_argument('-P', "--prewrite",
                     help=("Write first to a temporary directory then copy (sometimes helps with network storage)"),
                     default="")
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
aparser.add_argument('-v', "--verbose",
                     help="Increase verbosity (helpful for debugging)",
                     action="store_true")
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

result_selection = ["hum_protein", genometype, "sstart", "send", "indiv_humcover"]
if args.group_mode:
    result_selection.append("groupnum")

# Read in FASTA file to get human protein lengths
print("Reading in FASTA...")
ids = []
lengths = []
hprot_seqs = args.fasta_humprot
for rec in SeqIO.parse(hprot_seqs, "fasta"):
    ids.append(rec.id)
    lengths.append(len(rec.seq))
hseqs = pl.DataFrame({"hum_protein": ids, "length_hum": lengths}, schema={"hum_protein": pl.Utf8, "length_hum": pl.UInt64})

# Calculate human coverage per individual hit

print("Computing individual totals...")
results=results.with_columns(
    indiv_length = pl.col("send") - pl.col("sstart")
).join(
    hseqs, on="hum_protein", how="left"
).with_columns(
    indiv_humcover = pl.col("indiv_length") / pl.col("length_hum")
)
if (args.verbose):  print(results)

# For these purposes, we only care about the intervals covering human proteins,
# so we can collapse down all the bacterial proteins for now. Also take this
# time to get rid of very small bacterial proteins and cases where there are too
# few remaining counts. Finally, remove entries that are individually "too long"

print("Filtering results...")
minbaclength=args.minbaclength

if fullcover:
    final_results = results.filter(
        pl.col("length_bac") >= minbaclength
    ).filter(
        pl.col("indiv_humcover") >= cov_fraction
    ).with_columns(
        best_indiv_humcover=pl.col("indiv_humcover"),
        tot_humcover=pl.col("indiv_humcover"),
        tot_length = pl.col("indiv_length")
    ).unique()
    if (args.verbose):  print(final_results)
else: # Partial coverage
    resultsF = results.filter(
        pl.col("length_bac") >= minbaclength
    ).filter(
        pl.col("indiv_humcover") < cov_fraction
    ).with_columns(
        count = pl.col("bac_protein").count().over(["hum_protein", genometype])
    ).filter(
        pl.col("count") >= mincount
    ).select(
        result_selection
    ).unique()
    if (args.verbose):  print(resultsF)
    # Sort by starting and ending location within each hum_protein/genome pair.
    # Then, find the running maximum of the end ("send") locations.
    # If an sstart is ahead of the last running maximum, it's a new interval.
    # Getting the cumulative sum over the booleans of "new interval?" gives you an
    # interval group per hum_protein/genome pair.
    print("Computing intervals...")
    intervals = resultsF.sort("hum_protein", genometype, "sstart","send").with_columns(
            running_max = pl.col("send").cum_max().over(["hum_protein", genometype])
    ).with_columns(
        interval_group = (pl.col("sstart") > pl.col("running_max").shift(fill_value=0)).cum_sum().over(["hum_protein", genometype])
    )
    if (args.verbose):  print(intervals)
    # Now summarize first each interval group (collapsing to min and max), then each
    # distinct interval (since we now know these don't overlap) to get the total length
    print("Tallying intervals...")
    # If this is the second pass, we've organized our results into groups that we want to preserve
    if args.group_mode:
        totlength_groups = ["hum_protein", genometype, "groupnum"]
    else:
        totlength_groups = ["hum_protein", genometype]
    total_length = intervals.group_by(totlength_groups + ["interval_group"]).agg(
        length=pl.col("send").max() - pl.col("sstart").min()
    ).group_by(totlength_groups).agg(
        tot_length=pl.col("length").sum()
    )
    if (args.verbose):  print(total_length)
    total_humcover = total_length.join(hseqs, on="hum_protein", how="left").with_columns(
        tot_humcover = pl.col("tot_length") / pl.col("length_hum")
    )
    if (args.verbose):  print(total_humcover)
    # Now get "best" individual humcover
    indiv_humcover = resultsF.group_by(["hum_protein", genometype]).agg(
        best_indiv_humcover = pl.col("indiv_humcover").max()
    )
    if (args.verbose):  print(indiv_humcover)
    overall = indiv_humcover.join(total_humcover, on=["hum_protein", genometype], how="inner")
    if (args.verbose):  print(overall)
    if (args.verbose):  print(overall.with_columns(indiv_covdiff=pl.col("tot_humcover")-pl.col("best_indiv_humcover")).select(pl.max("indiv_covdiff")))
    if (args.verbose):  print(overall.select(pl.min("best_indiv_humcover")))
    if (args.verbose):  print(overall.select(pl.max("tot_humcover")))

    print("Getting final results...")
    # Merge two analyses, then filter for cases where total humcover is over
    # covdiff% more than best_indiv_humcover, humcover is at least cov_fraction, and
    # best_indiv_humcover is less than cov_fraction
    overallF = overall.filter(
        pl.col("best_indiv_humcover") < cov_fraction,
        pl.col("tot_humcover") >= cov_fraction,
        (pl.col("tot_humcover") - covdiff) > pl.col("best_indiv_humcover")
    )

    if (args.verbose):  print(overallF)
    # Keep only non-redundant columns or columns used for matching
    overall_selection = ["hum_protein", genometype, "tot_humcover", "tot_length", "best_indiv_humcover"]
    if args.group_mode:
        overall_selection.append("groupnum")
    overallF = overallF.select(
        overall_selection
    )
    if (args.verbose):  print(overallF)

    # Prune original results to only keep entries in overallF
    if args.group_mode:
        final_results = results.join(overallF, on=["hum_protein", genometype, "groupnum"], how="inner")
    else:
        final_results = results.join(overallF, on=["hum_protein", genometype], how="inner")
    if (args.verbose):  print(final_results)

# Write output

print("Writing results...")
if args.prewrite=="":
    final_results.write_ipc(args.outfile)
else:
    os.makedirs(args.prewrite, exist_ok=True)
    intermediate_file=os.path.join(args.prewrite, os.path.basename(args.outfile))
    final_results.write_ipc(intermediate_file)
    shutil.copy(intermediate_file, args.outfile)
    os.remove(intermediate_file)

# finally, test that output saved properly
test = pl.scan_ipc(args.outfile)
test_len = test.select(pl.count()).collect().item()
actual_len = final_results.select(pl.count()).item()
if (test_len != actual_len): raise Exception("Something went wrong writing to disk; results may not be trustworthy")
