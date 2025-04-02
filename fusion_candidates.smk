import os

from snakemake.utils import min_version
min_version("6.15.5")

himem_slurm_partition = "hugemem"
medmem_slurm_partition = "cpu"

# These rules will run on the same node from which the workflow is invoked
# (the login node) since compute nodes cannot access the internet to get data sets.
localrules:  get_uhgg_gffs, get_rep_to_src_map, get_uhgp90_data, get_uhgg_metadata


#TARGETS = ['HPAFDA', 'HPAmetab']
TARGETS = ['HumanUPR']
BAC_OVERLAPS = ['0.5', '0.67', '0.75']
#GENOME_TYPES = ['src', 'rep']  # 'source' vs 'representative' UHGG genomes
GENOME_TYPES = ['src']
#COVERAGES = ['part70', 'part80', 'full']
PARTCOVERAGES = ['part60', 'part70', 'part80']
FULLCOVERAGES = ['full60', 'full70', 'full80']
COVERAGES = PARTCOVERAGES + FULLCOVERAGES
#HITS_THRESHOLDS = ['all', 20000]   # Applied at 'samegenome' step to limit runtime
HITS_THRESHOLDS = [20000]   # Applied at 'samegenome' step to limit runtime


# OVERRIDE
num_nodes = 20
max_target_seqs = 5


# Only make phylogenies from the partial group
rule all:
    input:
        expand("data/processed/export/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_eggnog.tsv",
            bac_overlap=BAC_OVERLAPS,
            targetset=TARGETS,
            genometype=GENOME_TYPES,
            thresh=HITS_THRESHOLDS,
            length=[x for x in PARTCOVERAGES]),

        expand("data/processed/summary/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/tree_summary",
            bac_overlap=BAC_OVERLAPS,
            targetset=TARGETS,
            genometype=GENOME_TYPES,
            thresh=HITS_THRESHOLDS,
            length=[x for x in PARTCOVERAGES]),

        expand("data/processed/export/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_eggnog.tsv",
            bac_overlap=BAC_OVERLAPS,
            targetset=TARGETS,
            genometype=GENOME_TYPES,
            thresh=HITS_THRESHOLDS,
            length=[x for x in COVERAGES]),

        expand("data/processed/humcover3/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/part_humcover3.ipc",
            bac_overlap=BAC_OVERLAPS,
            targetset=TARGETS,
            genometype=GENOME_TYPES,
            thresh=HITS_THRESHOLDS,
            length=[x[-2:] for x in PARTCOVERAGES]),
        expand("data/processed/humcover3/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/full_humcover3.ipc",
            bac_overlap=BAC_OVERLAPS,
            targetset=TARGETS,
            genometype=GENOME_TYPES,
            thresh=HITS_THRESHOLDS,
            length=[x[-2:] for x in FULLCOVERAGES]),


# Locally executed data download rules

rule get_uhgp90_data:
    output: "data/raw/uhgp-90/uhgp-90.faa",
            "data/raw/uhgp-90/uhgp-90.tsv",
            "data/raw/uhgp-90/uhgp-90_eggNOG.tsv",
    params: url = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgp_catalogue/uhgp-90.tar.gz"
    shell:
        "cd data/raw; "
        "curl -O {params.url}; "
        "tar zxvf uhgp-90.tar.gz"


# Locally executed rule
# Note: this can be time-consuming.
# TODO: needs end-to-end validation (test full-download and flattening)
rule get_uhgg_gffs:
    output: directory("data/raw/UHGG_{genometype}genomes_GFFs"),
    shell: "scripts/get_uhgg_gffs {wildcards.genometype} {output}"


rule get_uhgg_metadata:
    output: "data/raw/genomes-all_metadata.tsv"
    shell:
        "cd data/raw; "
        "curl -O https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv"


# Data processing rules


rule targets_db:
    input: "data/raw/{targetset}/{targetset}_targets.fasta"
    output: directory("data/raw/{targetset}/database"),
    conda: "envs/blast.yaml",
    shell: "makeblastdb -in {input} -dbtype prot -out {output}/{wildcards.targetset}"


# Split single FASTA file containing sequence queries into multiple similarly-sized
# chunks to allow for better parallelization across cluster nodes.
checkpoint split_search_queries:
    input: ancient("data/raw/uhgp-90/uhgp-90.faa"),
    output: directory("data/processed/split_search_queries/{targetset}/UHGP-90_queries"),
    shell: "scripts/split_fasta.py -i {input} -o {output} -n {num_nodes}"


# intermediate rule
rule search_in_targets:
    input:
        dbdir = "data/raw/{targetset}/database",
        query_chunk = "data/processed/split_search_queries/{targetset}/UHGP-90_queries/uhgp-90_{chunk}.fasta",
    output: "data/processed/search_in_targets/{targetset}/blast_{chunk}.tsv"
    conda: "envs/blast.yaml",
    threads: 24
    shell:
        ("blastp -outfmt 6 -num_threads {threads} "
          "-db {input.dbdir}/{wildcards.targetset} "
          "-query {input.query_chunk} -max_target_seqs 5 "
          "-evalue 1e-4 -out {output}")


def aggregate_blast(wildcards):
    checkpoint_out = checkpoints.split_search_queries.get(**wildcards).output[0]
    return expand("data/processed/search_in_targets/{targetset}/blast_{chunk}.tsv",
                    targetset = wildcards.targetset,
                    chunk = glob_wildcards(os.path.join(checkpoint_out, "uhgp-90_{chunk}.fasta")).chunk)


# Concat all the BLAST output files together into a single TSV file.
rule merged_search:
    input: aggregate_blast,
    output: f"data/processed/merged_search/{{targetset}}/U90_in_{{targetset}}_{max_target_seqs}tseqs.tsv"
    shell: "cat $(dirname {input[0]})/*.tsv > {output}"


## Extract UHGP-90 bacterial sequence lengths from FASTA file and store them in
## a pickled data frame for use by a downstream filter stage.
rule baclengths:
    input: ancient("data/raw/uhgp-90/uhgp-90.faa"),
    output: "data/processed/baclengths/uhgp-90_lengths.pkl",
    shell:
        "scripts/compute_sequence_lengths.py --infile {input} -o {output}"



# Create an IPC database mapping all UHGP proteins to UHGP-90 cluster IDs.
rule prot_cluster_arrowdb:
    input: ancient("data/raw/uhgp-90/uhgp-90.tsv"),
    output: directory("data/processed/prot_cluster_arrowdb/"),
    resources:
        time = "48:00:00",
        # time = "4:00:00",
        # # slurm_partition = himem_slurm_partition,
        mem_mb = 64000
    shell: "scripts/uhgp_to_arrow.R -i {input} -o {output}"

#------------------------------------------------------------------------------

# Select full-length search hits on the bacterial side.
rule fullbac:
    input:
        data = f"data/processed/merged_search/{{targetset}}/U90_in_{{targetset}}_{max_target_seqs}tseqs.tsv", # IDENTICAL
        lengths = "data/processed/baclengths/uhgp-90_lengths.pkl",   # IDENTICAL
    output: "data/processed/fullbac/{targetset}/{targetset}_{bac_overlap}_fullbac.pkl",
    shell:
        "scripts/fusion_candidates_baclength.py "
        "-s {input.lengths} "
        "-i {input.data} "
        "-b {wildcards.bac_overlap} "
        "-o {output}"

# Dev note: adapted to remove the clustermap step
# Determine which genomes encode distinct bacterial proteins matching input set
rule samegenome:
    input: data = "data/processed/fullbac/{targetset}/{targetset}_{bac_overlap}_fullbac.pkl",
           arrowmap = "data/processed/prot_cluster_arrowdb/"
    output: "data/processed/samegenome/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}/samegenome.ipc"
    threads: 20
    resources:
        # partition = himem_slurm_partition,
        time = "1:00:00",
        mem_mb = 224000
    run:
        thresh_flag = ""
        if wildcards.thresh != "all":
            thresh_flag = f"-t {wildcards.thresh} "
        shell("POLARS_MAX_THREADS={threads} "
              "scripts/samegenome_polars.py "
              "-d {input.data} "
              "--arrowmap {input.arrowmap} "
              "--genometype {wildcards.genometype} "
              f"{thresh_flag}"
              "-o {output}")

# Keep hits where coverage of human gene is covered by a collection of bacterial protein hits
# from the same genome or a single bacterial sequence, as requested.
rule humcover:
    input:
        data = "data/processed/samegenome/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}/samegenome.ipc",
        seqs = "data/raw/{targetset}/{targetset}_targets.fasta",
    output: "data/processed/humcover/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_humcover.ipc",
    threads: 20
    resources:
        mem_mb = 192000
    run:
        coverage_flag = ""
        part_percent_flag = f"-p {wildcards.length[-2:]} " # changed logic to allow customization of threshold
        if wildcards.length.startswith("full"):
            coverage_flag = "--coverage_full "
        shell("POLARS_MAX_THREADS={threads} "
              "scripts/humcover_polars.py "
              "-d {input.data} "
              "-f {input.seqs} "
              "-I "
              f"{coverage_flag}"
              f"{part_percent_flag}"
              "-o {output}")

# Next parts of the pipeline fork for full vs. part, because for full-length we skip localizing the proteins

# Integrate information from GFFs to get feature numbers, strand, and contig info
rule location_part:
    input:
        data = "data/processed/humcover/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{length}_humcover.ipc",
        gffdir = ancient(f"data/raw/UHGG_srcgenomes_GFFs"),
    output:
        "data/processed/location/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{length}_location.ipc",
    resources:
        mem_mb = 64000
    threads: 8
    run:
        shell("POLARS_MAX_THREADS={threads} scripts/localize_hits -d {input.data} -g {input.gffdir} -o {output}")

# For each remaining bacterial protein "hit" to a human protein, get its
# distance to the closest bacterial gene on the same contig on the same strand.
# Lose any bacterial hits with no other hits on the same contig/strand. Also
# filter by maximum distance (typically 3, can turn this off by setting -m very
# large if desired) and organize data into "feature groups" for next step
rule distance_part:
    input:
        "data/processed/location/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{length}_location.ipc",
    output:
        "data/processed/distance/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{length}_distance.ipc",
    resources:
        mem_mb = 64000,
        time = "1:00:00"
    threads: 8
    run:
        shell("POLARS_MAX_THREADS={threads} scripts/distance_polars.py -d {input} -o {output} -m 3")

# Since we've now lost some entries after this distance filtering, repeat
# humcover step on groups
rule humcover2:
    input:
        data = "data/processed/distance/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{length}_distance.ipc",
        seqs = "data/raw/{targetset}/{targetset}_targets.fasta",
    output: "data/processed/humcover2/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{length}_humcover.ipc",
    threads: 20
    resources:
        mem_mb = 128000
    run:
        coverage_flag = ""
        part_percent_flag = f"-p {wildcards.length[-2:]} " # changed logic to allow customization of threshold
        if wildcards.length.startswith("full"):
            coverage_flag = "--coverage_full "
        shell("POLARS_MAX_THREADS={threads} "
              "scripts/humcover_polars.py "
              "-d {input.data} "
              "-f {input.seqs} "
              "-I "
              "-G "    # grouped this time into feature sets
              f"{coverage_flag}"
              f"{part_percent_flag}"
              "-o {output}")

### Associate taxonomy for representative genomes per the UHGG archive's 'all_metadata' file.
rule taxonomy_part:
    input:
        data ="data/processed/humcover2/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{length}_humcover.ipc",
        tax = "data/raw/genomes-all_metadata.tsv",
    output: "data/processed/taxonomy/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{length}_taxonomy.ipc",
    resources:
        mem_mb = 32000,
        time = "1:00:00"
    shell:
       "scripts/fusion_candidates_taxonomy.py -d {input.data} -t {input.tax} -o {output}"

### If we're doing full-length, we can skip all the stuff about them needing to
### be a certain distance, and we can also collapse duplicates which should save
### a lot of space
rule taxonomy_full:
    input:
        data = "data/processed/humcover/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_full{length}_humcover.ipc",
        tax = "data/raw/genomes-all_metadata.tsv",
    output: "data/processed/taxonomy/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_full{length}_taxonomy.ipc",
    resources:
        mem_mb = 172000,
        time = "1:00:00"
    shell:
       "scripts/fusion_candidates_taxonomy.py -C -d {input.data} -t {input.tax} -o {output}"

# Full/part streams now merge

# Annotate search hits with eggNOG orthology IDs taken from the UHGP
# database file for the given cluster threshold.
rule eggnog:
    input:
        data = "data/processed/taxonomy/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_taxonomy.ipc",
        eggnog = ancient("data/raw/uhgp-90/uhgp-90_eggNOG.tsv"),
    output: "data/processed/eggnog/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_eggnog.ipc",
    threads: 6
    resources:
        mem_mb = 172000
    shell:
        "POLARS_MAX_THREADS={threads} scripts/annotate_eggnog.py -d {input.data} -e {input.eggnog} -o {output}"


# Export the final data frames as TSV for use by other tooling.
rule export:
    input: "data/processed/eggnog/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_eggnog.ipc",
    output: "data/processed/export/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_eggnog.tsv",
    resources:
        mem_mb = 150000
    shell: "POLARS_MAX_THREADS=1 scripts/convert_df -i {input} -f tsv -o {output}"

# Select only best hits for full, after filtering out "suspiciously good" hits (90%+ ID)
rule preprocess:
    input:
        full="data/processed/export/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_full{ln}_eggnog.tsv",
        part="data/processed/export/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_part{ln}_eggnog.tsv"
    output:
        full="data/processed/preprocess/{targetset}_{bac_overlap}_{genometype}_{thresh}_{ln}/full_no_contaminants.ipc",
        part="data/processed/preprocess/{targetset}_{bac_overlap}_{genometype}_{thresh}_{ln}/part_no_contaminants.ipc"
    resources:
        mem_mb = 150000
    shell: "scripts/preprocess-large-files.R {input.full} {input.part} $(dirname {output.part})"

rule best_hits_only:
    input: "data/processed/preprocess/{targetset}_{bac_overlap}_{genometype}_{thresh}_{ln}/full_no_contaminants.ipc"
    output: "data/processed/bho/{targetset}_{bac_overlap}_{genometype}_{thresh}_{ln}/full_no_contaminants.ipc"
    resources: mem_mb = 150000
    shell: "POLARS_MAX_THREADS=1 scripts/best_hits_only.py -d {input} -o {output}"

rule best_hits_only_dummy:
    input: "data/processed/preprocess/{targetset}_{bac_overlap}_{genometype}_{thresh}_{ln}/part_no_contaminants.ipc"
    output: "data/processed/bho/{targetset}_{bac_overlap}_{genometype}_{thresh}_{ln}/part_no_contaminants.ipc"
    resources: mem_mb = 150000
    shell: "ln -s $(realpath {input}) $(realpath {output})"

rule humcover3:
    input:
        data = "data/processed/bho/{targetset}_{bac_overlap}_{genometype}_{thresh}_{ln}/{pf}_no_contaminants.ipc",
        seqs = "data/raw/{targetset}/{targetset}_targets.fasta"
    output: "data/processed/humcover3/{targetset}_{bac_overlap}_{genometype}_{thresh}_{ln}/{pf}_humcover3.ipc"
    threads: 20
    resources:
        mem_mb = 128000
    run:
        coverage_flag = ""
        part_percent_flag = f"-p {wildcards.ln[-2:]} " # changed logic to allow customization of threshold
        if wildcards.pf.startswith("full"):
            coverage_flag = "--coverage_full "
        shell("POLARS_MAX_THREADS={threads} "
              "scripts/humcover_polars.py "
              "-d {input.data} "
              "-f {input.seqs} "
              "-I "
              f"{coverage_flag}"
              f"{part_percent_flag}"
              "-o {output}")

# Splits into multiple output files within the output directory
checkpoint prep_align:
    resources:
        mem_mb = 32000
    input: data="data/processed/eggnog/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_eggnog.ipc",
           bac_seqs=ancient("data/raw/uhgp-90/uhgp-90.faa"),
           target_seqs="data/raw/{targetset}/{targetset}_targets.fasta",
    output: directory("data/processed/prep_align/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}"),
    shell: "scripts/prepare_sorted_clustalo_inputs "
           "-d {input.data} "
           "-b {input.bac_seqs} "
           "-g {input.target_seqs} "
           "-o {output}"


# Perform multi-sequence alignment
rule align:
    input: "data/processed/prep_align/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}_and_opcands.fasta",
    output: "data/processed/align/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}.aln",
    conda: "envs/clustalo.yaml",
    threads: 20
    resources:
        mem_mb = 32000
    shell:
        "clustalo -i {input} > {output}"


# Filter alignments
rule filt_align:
    input: "data/processed/align/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}.aln",
    output: "data/processed/filt_align/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}.filt"
    shell: "scripts/filter_alignment "
           "-a {input} "
           "--coverage 0.7 "
           "-o {output}"


# Trim results of MSA according to the provided trim fraction.
rule trim:
    input: "data/processed/filt_align/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}.filt"
    output: "data/processed/trim/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}.trim"
    conda: "envs/clustalo.yaml",
    shell: "AMAS.py trim -i {input} -f fasta -d aa -t 0.6 -o {output}"


rule tree:
    input: "data/processed/trim/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}.trim"
    output: "data/processed/tree/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}.treefile",
    conda: "envs/trees.yaml",
    shell: "cd data/processed/tree/{wildcards.targetset}/{wildcards.targetset}_{wildcards.bac_overlap}_{wildcards.genometype}_{wildcards.thresh}_{wildcards.length}; "
           "ln -sf ../../../trim/{wildcards.targetset}/{wildcards.targetset}_{wildcards.bac_overlap}_{wildcards.genometype}_{wildcards.thresh}_{wildcards.length}/{wildcards.humprot_id}.trim .; "
           "iqtree -s {wildcards.humprot_id}.trim -T AUTO -m TEST -B 1000 -pre {wildcards.humprot_id}"


rule annotate_tree:
    input: tree="data/processed/tree/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}.treefile",
           data="data/processed/eggnog/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}_eggnog.ipc"
    output: "data/processed/annotate_tree/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}-tax.treefile",
    conda: "envs/trees.yaml",
    threads: 1
    shell: "scripts/annotate_tree.py -t {input.tree} -d {input.data} -o {output}"


rule render_tree:
    input: tree="data/processed/annotate_tree/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}-tax.treefile",
    output: "data/processed/render_tree/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{humprot_id}-tree.pdf",
    conda: "envs/trees.yaml",
    threads: 1
    shell: "QT_QPA_PLATFORM=offscreen  "
           "scripts/render_tree.py "
           "-t {input.tree} "
           "-T \"{wildcards.humprot_id} - {wildcards.targetset}_{wildcards.bac_overlap}_{wildcards.genometype}_{wildcards.thresh}_{wildcards.length}\" "
           "-o {output}"


def aggregate_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the scatter step
    '''
    checkpoint_output = checkpoints.prep_align.get(**wildcards).output[0]
    return expand('data/processed/render_tree/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/{i}-tree.pdf',
                    targetset = wildcards.targetset,
                    bac_overlap = wildcards.bac_overlap,
                    genometype = wildcards.genometype,
                    thresh = wildcards.thresh,
                    length = wildcards.length,
                    i=glob_wildcards(os.path.join(checkpoint_output, '{i}_and_opcands.fasta')).i)


rule tree_summary:
    input: aggregate_input,
    output: "data/processed/summary/{targetset}/{targetset}_{bac_overlap}_{genometype}_{thresh}_{length}/tree_summary"
    shell: "echo {input} > {output}"

