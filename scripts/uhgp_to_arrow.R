#!/usr/bin/env Rscript

library(tidyverse)
library(arrow)
library(optparse)
library(ape)

main <- function(opt) {
  tsvd <- open_tsv_dataset(opt$input,
                           schema=schema(field("clusterID", string()),
                                         field("proteinID", string())))
  tsva <- tsvd %>% mutate(proteinID = gsub("GUT_GENOME", "", proteinID)) %>%
    mutate(proteinID = gsub("_", "", proteinID)) %>%
    mutate(clusterID = gsub("GUT_GENOME", "", clusterID)) %>%
    mutate(clusterID = gsub("_", "", clusterID)) %>%
    mutate(proteinID = cast(proteinID, uint64())) %>%
    mutate(clusterID = cast(clusterID, uint64())) %>%
    mutate(clusterG = as.integer(clusterID / 1e9)) %>%
    #mutate(genome = floor(clusterID / 100000)) %>%
    #mutate(genome = cast(genome, uint32())) %>%
    #select(clusterID, proteinID, genome, clusterG) %>%
    select(clusterID, proteinID, clusterG) %>%
    group_by(clusterG)
  write_dataset(tsva, opt$output, format="arrow")
}

if (sys.nframe() == 0) {
  option_list = list(
    make_option(c("-i", "--input"), type="character",
                default=NULL,
                help="input UHGP cluster file",
                metavar="character"),
    make_option(c("-o", "--output"), type="character",
                default=NULL,
                help="directory to output arrow file",
                metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  main(opt)
}
