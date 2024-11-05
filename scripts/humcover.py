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
aparser.add_argument('-p', "--partpercent",
                     type=int,
                     help="Coverage fraction of partial human protein matches")
args = aparser.parse_args()

if args.coverage_full:
    fullcover = True
else:
    fullcover = False
    if args.partpercent:
        cov_fraction = args.partpercent / 100
        print(f"cov_fraction = {cov_fraction}")
    else:
        print("If --coverage_full not supplied, a --partpercent value must be provided.")
        sys.exit(1)

if args.IPC:
    results = pl.read_ipc(args.dataframe).to_pandas()
else:
    results = pd.read_pickle(args.dataframe)

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
else:
    print("No genome type column found in data frame.")
    sys.exit(1)


# 3) Keep only hits where most of human protein is covered
#----------------------------------------------------------
# Get human protein lengths from drug targets FASTA file
#  add as column value for each human protein
# Compute (humprot_hit_length / humprot_length) and add as column  'humprot_coverage' to hits data frame.
def merge_intervals( intervals ):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


ids = []
lengths = []
hprot_seqs = args.fasta_humprot
for rec in SeqIO.parse(hprot_seqs, "fasta"):
    ids.append(rec.id)
    lengths.append(len(rec.seq))
hseqs = pd.DataFrame(list(zip(ids, lengths)), columns=["hum_protein", "length_hum"])

print("Associating human sequence lengths...")
hcov_hits = results.merge(hseqs,
                          how="left",
                          on="hum_protein")

print("Calculating human sequence coverage fraction...")
# Keep only sets of hits from a given genome where the set provides coverage of
# most of the human protein, but where individual bacterial protein hits do not.


def coverage(hphits):
    khits = pd.DataFrame()
    hprot = hphits["hum_protein"].iloc[0]
    print(hprot)
    sys.stdout.flush()
    for ggrp in hphits.groupby(genometype):

        # Keep only hits with bacterial sequence length > 100
        hits = ggrp[1].copy()
        hits = hits[ hits["length_bac"] > 100 ]

        hits['humcover_frac'] = ((hits["send"] - hits["sstart"]) / hits["length_hum"])
        hits = hits[ hits['humcover_frac'] < 0.75 ]

        # If only a single hit remains after the above filtering, skip this genome entirely.
        if hits.shape[0] < 2:
            continue

        # compose list of all intervals
        intervals = list( hits[["sstart", "send"]].apply(tuple, axis=1) )

        # merge all intervals into new list of non-overlapping intervals
        merged_intvs = merge_intervals(intervals)

        # sum up lengths of all intervals (numerator)
        accum = 0
        for intv in merged_intvs:
            accum += (intv[1] - intv[0])

        # calculate human protein coverage fraction
        hprot_coverage = accum / hits["length_hum"].iloc[0]

        # Compare collective human protein coverage of these hits with
        # the maximum individual humcover_frac found in this set.
        # If the difference between these is less than 0.2, ignore this set
        # of hits.
        max_humcover_frac = hits['humcover_frac'].max()
        covdiff = hprot_coverage - max_humcover_frac
        if covdiff < 0.2:
            continue

        # if human protein coverage fraction is > 0.8, keep all hits from this genome.
        # Q: All hits from this genome, or a selection of hits from this genome?
        #if hprot_coverage >= 0.8:   # TODO: Examine: Too stringent?
        if hprot_coverage >= cov_fraction:
            khits = khits.append(hits)

    #shortprot = hprot.split('|')[2].split('_')[0]
    #khits.to_pickle(f"humcover_{shortprot}.pkl")
    return khits


if fullcover:
    keephits = hcov_hits[ ((hcov_hits["send"] - hcov_hits["sstart"]) / hcov_hits["length_hum"]) > 0.75 ]
else:
    keephits = pd.DataFrame()
    hprot_grps = hcov_hits.groupby("hum_protein")
    groups = [hprot_grps.get_group(x) for x in hprot_grps.groups]
    print(f"Number of hum_protein groups: {len(groups)}")
    pool = Pool()
    hitsets = pool.map(coverage, groups)
    pool.close()

    print("Aggragating results...")
    for i, hitset in enumerate(hitsets):
       keephits = keephits.append(hitset)
    print(keephits)


keephits = keephits.reset_index(drop=True)
keephits.to_pickle(args.outfile)
