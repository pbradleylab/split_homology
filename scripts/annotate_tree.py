#!/usr/bin/env python3

import argparse
import pandas as pd
import polars as pl
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace

aparser = argparse.ArgumentParser(description="Annotate the leaf nodes of the given tree of genomes with taxonomic names.")

aparser.add_argument('-t', "--treefile", help="Tree file in Newick format.")
aparser.add_argument('-d', "--dataframe", help="Pickled data frame containing filtered search hits")
aparser.add_argument('-o', "--outfile", help="Output tree file. Will be written in Newick format.")

args = aparser.parse_args()

if args.dataframe.endswith(".ipc"):
    df = pl.read_ipc(args.dataframe).to_pandas()
elif args.dataframe.endswith(".pkl"):
    df = pd.read_pickle(args.dataframe)
else:
    raise ValueError

# Extract the genome type from the data frame columns.  Determined by samgenome step.
dfcols = df.columns
if "src_genome" in dfcols:
    print("Source genomes found")
    genometype = "src_genome"
elif "rep_genome" in dfcols:
    print("Representative genomes found")
    genometype = "rep_genome"
else:
    print("No genome type column found in data frame.")
    sys.exit(1)

tree = Tree(args.treefile)
ts = TreeStyle()

for t in tree.traverse():
    if "GUT_GENOME" in t.name:
        genome = t.name[:16]
        t.name = f"{genome} - {df[ df[genometype] == genome]['Lineage'].iloc[0]}"

tree.write(outfile=args.outfile)


