#!/usr/bin/env python3

import os
import subprocess as sp
import argparse
import csv
import lmdb

aparser = argparse.ArgumentParser(description=("Create a LMDB database containing "
                " _compact_ UHGP-xx protein map data.\n NOTE: Time-consuming. Performance "
                "can be improved by moving the input TSV file to a ramdrive (i.e. "
                "/dev/shm) first. The IDs will all be reduced to only the numerical"
                "portion to avoid storing redundant ID prefix strings."),
                formatter_class=argparse.RawTextHelpFormatter)
aparser.add_argument('-i', "--infile", required=True, help="UHGP-xx cluster map TSV file.")
aparser.add_argument('-o', "--outname", required=True, help=("Location of "
                        "directory to create to hold the database."))

args = aparser.parse_args()

# For a ~27GB input set, this map_size works.
env = lmdb.open(args.outname, map_size=int(1E12))
txn = env.begin(write=True)

print("Counting input file lines...")
proc = sp.run(f"wc -l {args.infile}".split(),
                stdout=sp.PIPE,
                stderr=sp.PIPE,
                universal_newlines=True)
totrows = int(proc.stdout.split()[0])
rowcount = 1

# Store each record as (protein ID, cluster ID) such that
# a query against the DB for a protein ID will produce the
# ID of the cluster into which it was grouped.
# Input file has pairs ordered (UHGP cluster ID, protein ID).
with open(args.infile, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        # To minimize sizse, store cluster and protein ID values as
        # 11-character strings of digits.
        #  ......*****
        #  First 6 digits (.) are the genome ID
        #  Last  5 digits (*) are the protein ID suffix
        protid = row[1].replace("GUT_GENOME","").replace("_","").encode()
        u90id = row[0].replace("GUT_GENOME","").replace("_","").encode()
        txn.put(protid, u90id)
        if rowcount % 1000000 == 0:
            print(f"{rowcount}/{totrows}  {(rowcount/totrows)*100:0.4}%")
        rowcount += 1

print("Committing transaction...")
txn.commit()
env.close()

