#!/usr/bin/env Rscript

# Take UHGG pangenomes and convert them to an Arrow database that gives which
# genomes have which UHGP-90 clusters (for samegenome step)

library(tidyverse)
library(arrow)
library(optparse)

main <- function(opt) {
  set_cpu_count(opt$threads)
  message(sprintf("Opening arrow dataset %s", opt$u90))
  u90 <- open_dataset(opt$u90, format="arrow")
  message(sprintf("Reading pangenome file %s", opt$pangenome))
  pg <- read_csv(opt$pangenome, col_types=cols())
                                        # Convert into "long" data frame
  message(" -> pivoting...")
  pg_justg <- pg %>%
    select(-(Gene:`Avg group size nuc`))
  # Drop empty cells when converting to long to save on space
  pg_long <- pg_justg %>%
    mutate(across(everything(), function(.) replace(., .=="", NA))) %>%
    pivot_longer(cols = everything(),
                 names_to = "genome",
                 values_to = "proteins",
                 values_drop_na = TRUE) %>%
    filter(values != "") %>%
    separate_longer_delim(proteins, "\t") %>%
    filter(!is.na(proteins)) %>%
    mutate(tag = opt$tag)
  message(sprintf(" -> pivoted with %d genomes found", ncol(pg_justg)))
                                        # Prepare data frame for Arrow merge
  pg_arrow <- pg_long %>%
    mutate(proteinID = as.numeric(
             gsub("_", "",
                  gsub("GUT_GENOME", "", proteins)))) %>%
    arrow_table() %>%
    mutate(proteinID = cast(proteinID, uint64()))
  pg_mapped <- right_join(u90, pg_arrow, by="proteinID") %>%
    mutate(clusterID = cast(clusterID, uint64())) %>%
    collect()
                                        # Clean up results and count clusterIDs (some may appear multiple times/genome)
                                        # Actually, skip this since we want to preserve proteinID
                                        #  pg_trim_count <- pg_mapped %>%
                                        #    mutate(clusterID = cast(clusterID, uint64())) %>%
                                        #    select(genome, clusterID, tag) %>%
                                        #    group_by(genome, clusterID, tag) %>%
                                        #    count() %>%
                                        #    arrange(genome, -n, clusterID) %>%
                                        #    ungroup() %>%
                                        #    mutate(species = opt$species) %>%
                                        #    group_by(tag, species) %>%
                                        #    arrange(genome)

                                        # Write to appropriate directory, making sure each MGYG gets its own non-
                                        # overwritten file and that they get partitioned by the "tag" for
                                        # performance reasons
  output_dir <- opt$output
  message(sprintf(" -> writing to directory %s...", opt$output))
  write_dataset(pg_mapped,
                opt$output,
                format="arrow",
                existing_data_behavior="overwrite",
                basename_template=paste0("part-{i}.arrow"))
}

if (sys.nframe() == 0) {
  option_list = list(
    make_option(c("-a", "--u90"), type="character",
                default=NULL,
                help="UHGP-90 arrow file",
                metavar="character"),
    make_option(c("-p", "--pangenome"), type="character",
                default=NULL,
                help="path to pan-genome",
                metavar="character"),
    make_option(c("-g", "--tag"), type="character",
                default="TAG",
                help="string to set tag column to (default: TAG)",
                metavar="character"),
    make_option(c("-s", "--species"), type="character",
                default="SPECIES",
                help="string to set species column to (default: SPECIES)",
                metavar="character"),
    make_option(c("-o", "--output"), type="character",
                default=NULL,
                help="output Arrow file directory",
                metavar="character"),
    make_option(c("-n", "--firstn"), type="integer",
                default=13L,
                help="how many characters to take for the partition string",
                metavar="integer"),
    make_option(c("-t", "--threads"), type="integer",
                default=1,
                help="how many threads to use",
                metavar="integer"),
    make_option(c("-T", "--test"), type="logical",
                default=FALSE,
                help="limit to randomly sampled 200 species for testing?",
                metavar="logical")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  main(opt)
}
