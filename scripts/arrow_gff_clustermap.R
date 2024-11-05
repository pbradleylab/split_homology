#!/usr/bin/env Rscript

library(tidyverse)
library(arrow)
library(optparse)
library(ape)

read_gff <- function(gffF) {
  # grab ID parameter
  gff <- ape::read.gff(gffF)
  gff %>% mutate(ID = gsub(".*ID=([^;]+);.*", "\\1", attributes))
}

main <- function(opt) {
  cmap <- open_dataset(opt$cmap)
  gffs <- list.files(path=opt$gffs, pattern=".*.gff$", recursive=TRUE)
  names(gffs) <- basename(gffs)
  for (gi in 1:length(gffs)) {
    gn <- names(gffs)[gi]
    message(sprintf("%07d: %s", gi, gn))
    tag <- substr(gn, 1, opt$firstn)
    n <- gsub("\.gff$", "", gn)
    gff <- read_gff(gffs[gn])
    gff2 <- select(gff, ID) %>%
      mutate(ID = gsub("GUT_GENOME", "", ID)) %>%
      mutate(ID = gsub("_", "", ID)) %>%
      mutate(tag = tag, genome = n) %>%
      select(tag, genome, ID) %>%
      rename(proteinID = ID)
    gffM <- left_join(gff2, cmap, by="proteinID") %>%
      collect() %>%
      group_by(tag)
    write_dataset(gffM,
                  opt$output,
                  format="arrow",
                  existing_data_behavior = "overwrite",
                  basename_template = paste0("part-{i}-", gi, ".arrow"))
  }
}

if (sys.nframe() == 0) {
  option_list = list(
    make_option(c("-c", "--cmap"), type="character",
                default=NULL,
                help="UHGP arrow file",
                metavar="character"),
    make_option(c("-g", "--gffs"), type="character",
                default=NULL,
                help="directory for gffs",
                metavar="character"),
    make_option(c("-n", "--firstn"), type="integer",
                default=13L,
                help="how many characters to take for the partition string",
                metavar="character"),
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  main(opt)
}