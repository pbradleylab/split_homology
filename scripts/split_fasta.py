#!/usr/bin/env python3

import os
from math import ceil
import argparse
import subprocess as sp
from Bio import SeqIO

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


aparser = argparse.ArgumentParser(description=("Split a FASTA file into 'num_chunks' "
            "chunks for parallel processing."))
aparser.add_argument('-n', "--num_chunks", help="Number of chunks into which the input file will be split")
aparser.add_argument('-i', "--input", help="FASTA file to split")
aparser.add_argument('-o', "--outdir", help="Output directory to receive the files from the splitting")
args = aparser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
nprefix = os.path.basename(args.input).replace(".fasta", "").replace(".faa", "")
format_type = "fasta"
num_records = int(sp.getoutput(f"grep '>' {args.input} | wc -l"))
chunksize = ceil(num_records / int(args.num_chunks))
record_iter = SeqIO.parse(open(args.input), format_type) 
print(f"{num_records} records found")
for i, batch in enumerate(batch_iterator(record_iter, chunksize)):
    filename = os.path.join(args.outdir, f"{nprefix}_{i:03}.fasta")
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, format_type)
    print(f"Wrote {count} records to {filename}")

