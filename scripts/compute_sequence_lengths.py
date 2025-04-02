#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO


aparser = argparse.ArgumentParser(description="Calculate sequence lengths from the provided FASTA file"
                                    " and save a pickled data frame of the protein:length mapping.")
aparser.add_argument('-i', "--infile", required=True, help="FASTA file containing sequence data")
aparser.add_argument('-o', "--outfile", required=True, help="Pickled data frame output filename.")
args = aparser.parse_args()


# Collect all bacterial protein sequence lengths for use in next stage.
#   Do this a single time and read the cached results thereafter.
print("Calculating bacterial sequence lengths...")
ids = []
lengths = []
for rec in SeqIO.parse(args.infile, "fasta"):
    ids.append(rec.id)
    lengths.append(len(rec.seq))
seqs = pd.DataFrame(ids, columns=["genomeid"])
seqs['length'] = lengths
seqs = seqs.set_index("genomeid")
seqs.to_pickle(args.outfile)

