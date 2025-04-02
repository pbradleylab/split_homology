#!/usr/bin/env python3

from p_tqdm import p_map
import sys
import os
import argparse
import pandas as pd
from glob import glob
from multiprocessing import set_start_method
#set_start_method("spawn", force=True)
#set_start_method("fork", force=True)
from multiprocessing import Pool
import gffutils
import pickle
import shelve

from collections import namedtuple
import lmdb


aparser = argparse.ArgumentParser(description=("Generate cluster map per-genome "
    "for the given genome type (source or representative). "
    ""
    "Cluster map is dictionary-like:"
    "   - genome_value is the key"
    "   - a list of tuples is the value"
    "     tuples contain"
    "        - unclustered protein ('protein_val' - compact ID)"
    "        - cluster to which it belongs ('UHGP_val' - compact ID)"
    ))
aparser.add_argument('-G', "--genometype",
        required=True,
        help="Either 'src' for source "
        "(GUT_GENOME...) genomes or 'rep' for representative (MGYG...) genomes")
aparser.add_argument('-g', "--gffdir",
        help="Root of directory tree containing representative genome (MGYG...) GFF files.")
aparser.add_argument('-u', "--u90db",
        help=("Database location containing UHGP-90 protein cluster mapping for use"
        " when dealing with raw source genomes."))
aparser.add_argument('-p', "--processes", type=int,
        help=("Number of processes to spawn in parallel."))
aparser.add_argument('-o', "--outfile",
        help="Output file name to receive pickled cluster_map_per_genome list.")
aparser.add_argument('-c', "--copydb", action='store_true',
        help="If true, cache lmdb in storage under $TMPDIR")
args = aparser.parse_args()


num_processes = args.processes

genometype = args.genometype
if genometype == 'src':
    genome_col = 'src_genome'
    genid_prefix = "GUT_GENOME"
if genometype == 'rep':
    genome_col = 'rep_genome'
    genid_prefix = "MGYG-HGUT-"
    if not args.gffdir:
        print("Argument -g,--gffdir required for 'rep' genometype specified.")
        sys.exit(1)
    gffdir = args.gffdir


# Utility functions for converting the strings stored in the cluster map database
# to full protein/cluster ID names, and vice versa.
def expand(prot_val):
    return f"GUT_GENOME{prot_val[0:6]}_{prot_val[6:]}"

def compact(prot_id):
    return prot_id.replace("GUT_GENOME","").replace('_','')


# Cluster association named tuple for composing cluster_map_per_genome.
Classoc = namedtuple('Classoc', ['protein_val', 'UHGP_val'])

cluster_map_per_genome = []


# retain only entries from the DB where the full UHGP90 ID happens to also be in the blast results
# only original protein IDs that map to full UHGP90 IDs from the hits
# first step 
# subset the DB based on the content of the blast results
#    keep pairs where the protein ID matches the U90 ID from the hits
#   UHGP90 ID is in bac_protein (full length ID)
#   appears verbatim in the bac_protein column
#   yields a reduced list
#     protein   UHGP90 ID
#   take the stems of the protein ID column from this reduced DB set
#  
#
#  want to obtain every genome that has a protein in the
#   
# - Get list of unclustered protein IDs for each UHGP values that appear in the incoming hits table
# - Get list of unique genomes from this list of unclustered protein IDs.
#   "every genome that has a protein in these U90 clusters"
# - This list of genomes is then used in the existing logic below to make
#   the cluster_map_per_genome
# - Preserve the initial unclustered protein ID, at intersection phase
#   

# Below, 'vals' are compact protein identifiers "GGGGGGPPPPP" as strings
#        'ids' are full protein IDs  "GUT_GENOMEGGGGGGPPPPP"


# Open UHGP-90 inverse mapping database  UHGP_val key -> multiple protein vals
#ienv = lmdb.open("/fs/project/bradley.720/db/uhgp/uhgp-90/uhgp90_inv_compact.lmdb",
#        max_dbs=3,
#        map_size=int(1E12))
#db = ienv.open_db("multivalue".encode(), dupsort=True)
#itxn = ienv.begin(db=db)
#cursor = itxn.cursor()


# Open UHGP-90 cluster map database
if (args.copydb):
  tmpdir = os.environ["TMPDIR"]
  if not (os.path.exists(tmpdir)): raise IOError(f"Temporary directory {tmpdir} not found")
  print(f"Copying database to {tmpdir}...")
  os.system(f"cp -r {args.u90db} {tmpdir}")
  db_loc = os.path.join(tmpdir, os.path.basename(args.u90db))
else:
  db_loc = args.u90db
print(f"Loading database from {db_loc}...")
env = lmdb.open(db_loc, map_size=int(1E12), readonly=True)


if genometype == "src":

    # Generate genome ID values for all UHGG genomes.
    unique_genvals = [f"{x:06}" for x in range(1,286998)]

    numgenomes = len(unique_genvals)
    print(f"Number of unique genomes: {numgenomes}")


    def src_lookup(genome_val):
        '''Compose a dictionary with the genome value as key
        and a list of tuples as the value
        tuples contain
            - unclustered protein value (compact ID)
            - cluster value to which it belongs (compact ID)'''

        #print(f"{genome_val}  tot_genomes:{numgenomes}")
        with env.begin() as txn: # Open read-only LMDB transaction.
          clust_ids = []
          idx = 1
          while idx < 10000:  # No GFF files have protein IDs greater than this value.
              protval = f"{genome_val}{idx:05d}"
              idx += 1
              clustval = txn.get(protval.encode())
              if not clustval: # Account for potential gaps in protein ID sequence.
                  continue
              clustval = clustval.decode()
              classoc = Classoc(protein_val = protval, UHGP_val = clustval)
              clust_ids.append(classoc)
        return genome_val, clust_ids


    #pool = Pool(processes=num_processes)
    gcmaps = p_map(src_lookup, unique_genvals, num_cpus=num_processes)

else:
    # For each representative genome ID, iterate over all bacterial protein ID values
    # found in the rep genome's GFF file obtained from the EBI site archive and look up
    # the UHGP-90 cluster ID associated with each. Store these in the
    # cluster_map_per_genome dict that holds the protein ID and its associated UHGP ID
    # as a tuple.

    def rep_lookup(gfile):
        genome_id = os.path.basename(gfile).split(".gff")[0]
        genome_val = genome_id.split('-')[-1]
        print(genome_id)
        with env.begin() as txn: # Open read-only LMDB transaction.
          clust_ids = []
          gfdb = gffutils.create_db(str(gfile), ":memory:")
          for feat in gfdb.all_features():
              protval = compact(feat.id)
              if feat.id not in clust_ids:
                  clustval = txn.get(protval.encode())
                  if clustval:
                      clustval = clustval.decode()
                      classoc = Classoc(protein_val = protval, UHGP_val = clustval)
                      clust_ids.append(classoc)
        return genome_val, clust_ids 


    #pool = Pool(processes=num_processes)
    gcmaps = p_map(rep_lookup, glob(f"{args.gffdir}/{genid_prefix}*.gff"), num_cpus=num_processes)


#pool.close()


# genome ID  (key)
#     [ (unclust_compact_ID, UHGP_compact_ID), (unclust_compact_ID, UHGP_compact_ID), ... ]
print("Composing cluster_map_per_genome...")

num_maps = len(gcmaps)
for i, gcmap in enumerate(gcmaps):
    genomeid = gcmap[0]
    classoc_list = gcmap[1]
    cluster_map_per_genome.append( [genomeid, classoc_list] )
print(f"Finished creating cluster map per genome.")

print("Writing cluster_map_per_genome to disk...")
with open(args.outfile, 'wb') as f:
    pickle.dump(cluster_map_per_genome, f, protocol=pickle.HIGHEST_PROTOCOL)
