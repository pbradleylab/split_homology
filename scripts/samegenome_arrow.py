#!/usr/bin/env python3
import math
import psutil
import sys
import os
import argparse
import pickle
from datetime import datetime
import pandas as pd
import multiprocessing
from multiprocessing import Pool, Value, Manager
import shelve
import pyarrow as pa
import pyarrow.dataset as ds

from collections import namedtuple


process = psutil.Process()

aparser = argparse.ArgumentParser(description=("Select hits for proteins seen "
        "together in same genome."))
aparser.add_argument('-d', "--dataframe",
        help="Pickled data frame containing filtered search hits.")
aparser.add_argument('-p', "--processes",
                     default=1,
        help="Number of processes")
aparser.add_argument('-a', "--arrowmap",
        help="cluster_map UHGP arrow dataset")
aparser.add_argument('-G', "--genometype",
        required=True,
        help="Either 'src' for source "
        "(GUT_GENOME...) genomes or 'rep' for representative (MGYG...) genomes")
aparser.add_argument('-T', "--test",
                     default=False,
                     action='store_true',
        help="Do an abridged version for testing")
aparser.add_argument('-t', "--hits_threshold",
        required=False,
        type=int,
        help="Threshold of number of blast hits in a given human protein group above which "
        "the human protein is excluded from the results. If not supplied, all hits are "
        "used in this script (may impact processing time significantly).")
aparser.add_argument('-o', "--outdir",
        help="Output directory receiving multiple pickled dataframe files, one for each "
             "human protein group.")
args = aparser.parse_args()

# Create output directory if it does not exist
if os.path.exists(args.outdir):
    if not os.path.isdir(args.outdir):
        print(f"{args.outdir} is not a directory. Please check arguments.")
else:
    os.makedirs(args.outdir)

results = pd.read_pickle(args.dataframe)

# Cluster association named tuple for composing cluster_map_per_genome.
Classoc = namedtuple('Classoc', ['protein_val', 'UHGP_val'])


genometype = args.genometype
if genometype == 'src':
    genome_col = 'src_genome'
if genometype == 'rep':
    genome_col = 'rep_genome'


# TODO: generalize for 'src' as well as 'rep' genomes.
# Utility functions for converting the strings stored in the cluster map database
# to full protein/cluster ID names, and vice versa.
def expand(prot_val):
    return f"GUT_GENOME{prot_val[0:6]}_{prot_val[6:]}"

def compact(prot_id):
    return prot_id.replace("GUT_GENOME","").replace('_','')



# print("Reading cluster map list into dictionary...")
print("Reading cluster map list into list...")
sys.stdout.flush()
before = datetime.now()

# uhgp_ids_per_genome = pd.read_pickle(args.arrowmap)
cm_schema = pa.schema([
    ('clusterID', pa.uint64()),
    ('proteinID', pa.uint64())
])

# uhgp_ids_per_genome = []
# read into a list for lower memory usage, but still use the dict to add to lists efficiently
# by_genome_dict = dict()
dset = ds.dataset(args.arrowmap, format="arrow", schema=cm_schema)
dt = dset.to_table()

reader = dset.to_batches()
#
#for (i, batch) in enumerate(reader):
#    if args.test:
#        if i == 10000: break
#    bd = batch.to_pydict()
#    for pID, cID in zip(bd['proteinID'], bd['clusterID']):
#        gn = math.floor(pID / 100000)
#        g = f"GUT_GENOME{gn:06}"
#        p = f"{pID:011}"
#        c = f"{cID:011}"
#        classoc = Classoc(protein_val = p, UHGP_val = c)
#        if g in by_genome_dict.keys():
#            uhgp_ids_per_genome[by_genome_dict[g]][1].append(classoc)
#        else:
#            uhgp_ids_per_genome.append((g, [classoc]))
#            by_genome_dict[g] = len(uhgp_ids_per_genome) - 1
#    del bd
#del reader
#del dset
#del by_genome_dict # no longer need this once read in
#after = datetime.now()
#reader = dset.to_batches()
print(f" -- Memory usage (1): {process.memory_info().rss}")
udict = dict()
for (i, batch) in enumerate(reader):
    if args.test:
        if i == 10000: break
    bd = batch.to_pydict()
    for pID, cID in zip(bd['proteinID'], bd['clusterID']):
        gn = math.floor(pID / 100000)
        g = f"GUT_GENOME{gn:06}"
        p = f"{pID:011}"
        c = f"{cID:011}"
        classoc = Classoc(protein_val = p, UHGP_val = c)
        if g in udict.keys():
            #uhgp_ids_per_genome[by_genome_dict[g]][1].append(classoc)
            udict[g].append(classoc)
        else:
            udict[g] = [classoc]
            #uhgp_ids_per_genome.append((g, [classoc]))
            #by_genome_dict[g] = len(uhgp_ids_per_genome) - 1
    del bd
del reader
del dset
print(f" -- Memory usage (2): {process.memory_info().rss}")
uhgp_ids_per_genome = udict.items() # convert to list of tuples representation
del udict
print(f" -- Memory usage (3): {process.memory_info().rss}")
after = datetime.now()

print(f"Done.  Read data in {after-before}")
sys.stdout.flush()


def intersection_hits(hphits):
    '''Accepts group of hits against a particular human protein'''

    start = datetime.now()
    best = pd.DataFrame()
    humprot = hphits.iloc[0]['hum_protein']
    numhits = len(hphits)

    # Skip human protein groups that have too many hits in order to save some processing time.
    if args.hits_threshold and numhits > args.hits_threshold:
        stop = datetime.now()
        print(f"{humprot} - {numhits} hits - {stop-start} - Number of hits exceeds threshold. Skipping.")
        sys.stdout.flush()
        return

    protein_sets_Ax = []  # List of U90 protein sets matching BLAST hits
    genomes_Bx = []       # The genomes from which each of those sets come.
    unclust_Cx = []       # list of dicts
  
    hit_compact_ids = set(hphits['bac_protein'].apply(compact))

    for (key, classocs) in uhgp_ids_per_genome:
        #classocs = uhgp_ids_per_genome[key]
        genome_uhgp_compact_ids = set([cl[1] for cl in classocs])

        intersect = hit_compact_ids & genome_uhgp_compact_ids

        if len(intersect) > 1:
            # if args.test: print(f"Found hits for {humprot} in {key}")
            expanded_set = [expand(x) for x in intersect]
            protein_sets_Ax.append([expand(x) for x in intersect])
            #genomeid = expand(item[0])
            genomes_Bx.append(key)
            del expanded_set

        del genome_uhgp_compact_ids
        del intersect

    ## Keep the hits with U90 IDs that appear in each intersection and associate the
    ## Bx genome in question with each group of such hits.
    best_list = []
    for i, genomeid in enumerate(genomes_Bx):

        # Keep the hits with U90 IDs that appear in each intersection and associate the
        # Bx genome in question with each group of such hits.
        filt = hphits[ hphits["bac_protein"].isin(protein_sets_Ax[i]) ]
        filt = filt.reset_index(drop=True)
        filt[genome_col] = genomeid

        if filt.shape[0] < 2:   # Skip genomes with fewer than two total hits.
            continue

        # If multiple bacterial protein hits remain for a given genome and
        # human gene, keep only the best hit.
        filt = filt.iloc[ filt.groupby('bac_protein')['bitscore'].agg(pd.Series.idxmax) ]
        best_list.append(filt)#, ignore_index=True)
        del filt

    best = pd.concat(best_list)
    if args.test: print(best)
    del protein_sets_Ax
    del best_list

    if len(best) == 0:
        stop = datetime.now()
        print(f"{humprot} - {numhits} hits - {stop-start} - no hits after filtering. skipped.")
        return

    # Instead of returning data frame to starmap for concatenation and writing,
    # dump the result to disk ASAP and worry about putting it all back together later.
    short_humprot = humprot.split('|')[-1]
    best.to_pickle( os.path.join(args.outdir, f"{short_humprot}.pkl" ))
    stop = datetime.now()
    print(f"{humprot} - {numhits} hits - {stop-start} - wrote output")
    sys.stdout.flush()
    return(0)


hprot_grps = results.groupby("hum_protein")
orig_groups = [hprot_grps.get_group(x) for x in hprot_grps.groups]


# Ignore all groups with only a single hit
culled_groups = []
for group in orig_groups:
    if len(group) != 1:
        culled_groups.append(group)


# Sort groups by number of hits
print("Sorting groups...")
sys.stdout.flush()
sorted_groups = sorted(culled_groups, key = lambda l: l.shape[0], reverse=True)
if args.test: sorted_groups = sorted_groups[0:50]

numgroups = len(sorted_groups)
print(f"Number of hum prots: {numgroups}")
print(f" -- Memory usage (4): {process.memory_info().rss}")
print("Begin parallel processing...")
sys.stdout.flush()

if __name__ == '__main__':

    process2 = psutil.Process()
    with multiprocessing.get_context('fork').Pool(processes=int(args.processes), maxtasksperchild=1) as pool:
    #with multiprocessing.get_context('forkserver').Pool() as pool:
        pool.map(intersection_hits, sorted_groups)
        pool.close()
        pool.join()
    print(f" -- Memory usage (final): {process2.memory_info().rss}")
