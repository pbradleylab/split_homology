#!/usr/bin/env python3

import argparse
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace

aparser = argparse.ArgumentParser(description="Annotate the leaf nodes of the given tree of genomes with taxonomic names.")

aparser.add_argument('-t', "--treefile", help="Tree file in Newick format.")
aparser.add_argument('-T', "--title", help="Title to place on tree plot")
aparser.add_argument('-o', "--outfile", help="Output tree plot PDF file.")

args = aparser.parse_args()

tree = Tree(args.treefile)
ts = TreeStyle()

ts.title.add_face(TextFace(f"      {args.title}"), column=0)
tree.render(args.outfile, tree_style=ts)
