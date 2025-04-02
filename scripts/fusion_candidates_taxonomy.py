#!/usr/bin/env python3

import argparse
import pandas as pd
import polars as pl

aparser = argparse.ArgumentParser(description="Assign taxonomic ID for MGYG... organism IDs taken from the UHGG paper's supplental materials.")
aparser.add_argument('-d', "--dataframe", help="Pickled data frame containing filtered search hits.")
aparser.add_argument('-t', "--taxonomy", help="Path to genomes-all_metadata.tsv file from UHGP project")
aparser.add_argument('-o', "--outfile", help="Output file receiving pickled data frame with augmented search hits.")
aparser.add_argument('-C', "--collapse", help="Collapse and count identical hits (useful for full list)", action='store_true')
args = aparser.parse_args()

results = pl.scan_ipc(args.dataframe)

# Read taxonomy info for MGYG genome IDs and associate with genomes.
#tax = pd.read_excel(args.taxonomy,
#                    header=1,
#                    skiprows=[2])
tax = pl.scan_csv(args.taxonomy,
                  separator='\t')

# Extract the genome type from the data frame columns.  Determined by samgenome step.
dfcols = results.columns
if "src_genome" in dfcols:
    print("Source genomes found")
    genometype = "src_genome"
    tax = tax.select(["Genome", "Lineage"])
    tax = tax.rename({"Genome":"src_genome"})
elif "rep_genome" in dfcols:
    print("Representative genomes found")
    genometype = "rep_genome"
    tax = tax.select(["MGnify_accession", "Lineage"])
    tax = tax.unique()
    tax = tax.rename({"MGnify_accession":"rep_genome"})
elif "g" in dfcols:
    print("Generic genomes found")
    genometype = "g"
    tax = tax.select(["Genome", "Lineage"]).with_columns(
        g = pl.col("Genome").str.replace("GUT_GENOME","").str.replace("_", "").cast(pl.UInt32)
    )

else:
    print("No genome type column found in data frame.")
    sys.exit(1)

# Subset and adjust column names.
#tax = tax[["Species representative", "MGnify accession", "Taxonomy lineage (GTDB)"]]
#tax = tax.rename(columns = {"MGnify accession":"rep_genome"})
#tax = tax.rename(columns = {"Taxonomy lineage (GTDB)":"Taxonomy-GTDB"})
#results = results.merge(tax[["rep_genome","Taxonomy-GTDB"]], how="left", on="rep_genome")
results = results.join(tax, how="left", on=genometype)

if args.collapse:
    results = results.group_by(
        pl.exclude(["Genome"], genometype, "proteinID")
    ).agg(
        pl.exclude(["Genome"], genometype, "proteinID").unique(),
        nGenomes = pl.col("Genome").count()
    )

if args.outfile.endswith(".pkl"):
    results.collect().to_pandas().to_pickle(args.outfile)
elif args.outfile.endswith(".ipc"):
    results.collect(streaming=True).write_ipc(args.outfile)
