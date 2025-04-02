#!/usr/bin/env Rscript

# Take UHGG pangenomes and convert them to an Arrow database that gives which
# genomes have which UHGP-90 clusters (for samegenome step)

library(tidyverse)
library(arrow)
library(optparse)

main <- function(opt) {

  set_cpu_count(opt$threads)
  cmap <- open_dataset(opt$cmap, format="arrow")
  pangenomes <- list.files(path=opt$pangenomes,
                           pattern="genes_presence-absence_locus.csv$",
                           recursive=TRUE)
  # expected structure is "raw/UHGG_srcgenomes_GFFs/MGYG-HGUT-000/MGYG-HGUT-00001/pan-genome/genes...csv"
  names(pangenomes) <- basename(dirname(dirname(pangenomes)))
  # for running smaller tests

  # iterate one "tag" (MGYG-HGUT-???) at a time for efficiency
  tags <- basename(dirname(dirname(dirname(pangenomes))))
  ut <- unique(tags)
  if (opt$test) {
    ut = sample(ut, 5)
  }
  # Doing this instead of a "for" in the hope that it improves memory release
  lapply(ut, function(u) {

    message(sprintf("Tag %s (%d of %d)", u, which(ut==u)[1], length(ut)))
    indices <- which(tags == u)
    gc()

    tag <- u

    arrow_tbls <- map(indices, function(gi) {
      gn <- names(pangenomes)[gi]
      pangenomeF <- pangenomes[gi]
      whichN <- which(indices == gi)
      message(sprintf("  -> %07d: %s", whichN, gn))
      pg <- read_csv(file.path(opt$pangenomes, pangenomeF), col_types=cols())
                                        # Convert into "long" data frame
      pg_justg <- pg %>%
        select(-(Gene:`Avg group size nuc`))
      pg_long <- pg_justg %>%
        pivot_longer(cols = everything(),
                     names_to = "genome",
                     values_to = "proteins") %>%
        separate_longer_delim(proteins, "\t") %>%
        filter(!is.na(proteins)) %>%
        mutate(tag = tag)
      message(sprintf("    -> %d genomes found", ncol(pg_justg)))
                                        # Prepare data frame for Arrow merge
      pg_arrow <- pg_long %>%
        mutate(proteinID = as.numeric(
                 gsub("_", "",
                      gsub("GUT_GENOME", "", proteins)))) %>%
        arrow_table() %>%
        mutate(proteinID = cast(proteinID, uint64()))
      pg_mapped <- right_join(cmap, pg_arrow, by="proteinID") %>% collect()
                                        # Clean up results and count clusterIDs (some may appear multiple times/genome)
      pg_trim_count <- pg_mapped %>%
        select(genome, clusterID, tag) %>%
        group_by(genome, clusterID, tag) %>%
        count() %>%
        arrange(genome, -n, clusterID) %>%
        ungroup() %>%
        arrange(genome)
      return(pg_trim_count)
    })
                                        # Write to appropriate directory, making sure each MGYG tag gets its own file
    arrow_bound <- bind_rows(arrow_tbls) %>%
      group_by(tag)
    write_dataset(arrow_bound,
                  opt$output,
                  format="arrow",
                  existing_data_behavior="overwrite",
                  basename_template=paste0("part-{i}-", tag, ".arrow"))
  })
}

if (sys.nframe() == 0) {
  option_list = list(
    make_option(c("-c", "--cmap"), type="character",
                default=NULL,
                help="UHGP arrow file",
                metavar="character"),
    make_option(c("-p", "--pangenomes"), type="character",
                default=NULL,
                help="directory for pan-genomes",
                metavar="character"),
    make_option(c("-o", "--output"), type="character",
                default=NULL,
                help="output directory for Arrow database",
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
