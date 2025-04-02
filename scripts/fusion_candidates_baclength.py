#!/usr/bin/env python3

import argparse
import pandas as pd


# Read BLAST search output into a data frame.
# Run first filter stage used to identify possible bacterial operon fusions into human protein coding genes.
#    Keep only hits that are full-length bacterial sequence matches
#      'full-length' being 80% in this case.

colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "quend", "sstart", "send", "evalue", "bitscore"]

pd.options.display.width = 1200
pd.set_option('display.max_columns', None)

aparser = argparse.ArgumentParser(description="Filter BLAST results for horizontal gene transfer project.")
aparser.add_argument("-i", "--infile", help="BLAST-6 TSV input file containing search hits.")
aparser.add_argument("-s", "--seq_lengths", help="Pickled data frame containing bacterial sequence lengths.")
aparser.add_argument("-b", "--bac_overlap", type=float, help="Bacterial sequence overlap fraction.")
aparser.add_argument('-o', "--out", help="Output file receiving pickled data frame with augmented search hits.")
args = aparser.parse_args()

results = pd.read_csv(args.infile, sep='\t')
results.columns = colnames


# Determine which column represents the human genes, as this changes depending upon how the BLAST
# search was structured. i.e.
#  search for druggable targets in UHGP-50 database
#   - or -
#  search for UHGP-50 genes within the druggable targets database
# These two approaches provide different results because the behavior of the BLAST tool is
# such that multiple hits may be returned for each query sequence and the number of query
# sequences differs depending on which side of this comparison is used as the database.
if "HUMAN" in results["qseqid"][0]:
    results = results.rename(columns={"qseqid":"hum_protein", "sseqid":"bac_protein"})
    results = results.rename(columns={"sstart":"bacstart", "send":"bacend"})
else:
    results = results.rename(columns={"qseqid":"bac_protein", "sseqid":"hum_protein"})
    results = results.rename(columns={"qstart":"bacstart", "quend":"bacend"})


# 1) Hits must be full-length on bacterial side
seqs = pd.read_pickle(args.seq_lengths)


# Perform a merge to associate the collected bacterial sequence length in seqs dataframe
# with the corresponding hit row in the search results so that overlap fraction
# may be easily computed.
# NOTE: Time-consuming.
print("Merging in bacterial protein sequence lengths...")
results = results.merge(seqs,
                        how='left',
                        left_on="bac_protein",
                        right_index=True,
                        suffixes=[None, "_bac"])

overlap = args.bac_overlap
print(f"Calculating bacterial sequence overlap fraction. Removing hits with overlap < {overlap:0.0%}.")
# TODO: change to overlap_frac_bac, to distinguish it.
results['overlap_frac'] = (results["bacend"] - results["bacstart"]) / results["length_bac"]
results = results[ results['overlap_frac'] >= overlap ]
num_genes = len(results['hum_protein'].unique())
#print(f"Human genes with hits having >0.8 bacterial overlap fraction: {len(hashits)}/{num_genes}")
results.to_pickle(args.out)

