# Pipeline Data Inputs

## Inputs that must be manually acquired

Where possible, crucial data is acquired via rules within the pipeline workflow
definition. 

The human target sets need to be downloaded manually due to the way in which they
are selected from their source.

### Human Protein Target Sequences

#### Update: run download.sh in data/raw/ to get the full set of UniProt sequences.
 
## Automatically Acquired Inputs

### Bacterial Protein Sequences

The pre-generated clustering identity threshold selected for this investigation
is 90%, thus the corresponding UHGP-90 data set is used.

This UHGP-90 data is acquired by the pipeline rule get_uhgp90_data, but may
take some time as the download pulls several tens of GB of data.

This input data will be placed in data/raw/uhgp-90/uhgp-90.faa
                                  data/raw/uhgp-90/uhgp-90.tsv
                                  data/raw/uhgp-90/uhgp-90_eggNOG.tsv
by the "get_uhgp90_data" workflow rule.

### Bacterial Genomes

Bacterial protein data is taken from version 1.0 of Unified Human
Gastrointestinal Genome Project (UHGP) data set to serve as inputs for the
bacterial side of the search and filtering. 

Two such UHGG groupings are acquired. The first, termed "representative
genomes", is where the protein IDs are clustered into 4644 unique species
genomes.  The second, termed "source genomes" is the set of all bacterial
genome IDs in the set; 286,997 separate genomes, in this case.

These input data sets will be placed in data/raw/UHGG_repgenomes_GFFs
                                        data/raw/UHGG_srcgenomes_GFFs
by the "get_uhgg_gffs" workflow rule.

### UHGG Metadata

For taxonomic associations to be made in the pipeline, some metadata for the UHGG data needs to be obtained.

These input data will be placed in data/raw/genomes-all_metadata.tsv
by the "get_uhgg_metadata" workflow rule.


