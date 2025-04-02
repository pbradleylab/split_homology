#!/bin/bash

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
gunzip UP*gz
mv UP000005640_9606.fasta HumanUPR_targets.fasta
