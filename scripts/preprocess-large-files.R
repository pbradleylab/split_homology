#!/usr/bin/env Rscript
library(tidyverse)
library(arrow)

args <- commandArgs(trailingOnly = TRUE)

eggnog_full <- read_tsv(args[1])
eggnog_part <- read_tsv(args[2])
output_dir <- args[3]
#eggnog_full <- read_tsv("HumanUPR_0.67_src_20000_full70_eggnog.tsv")
#eggnog_part <- read_tsv("HumanUPR_0.67_src_20000_part70_eggnog.tsv")

max_reasonable_pident <- mean(eggnog_full$pident) + sd(eggnog_full$pident)*3
write(x = max_reasonable_pident, file=file.path(output_dir, "max-reasonable-pident.txt"))

eggnog_full_contams <- eggnog_full %>%
  filter(pident > max_reasonable_pident)
eggnog_part_contams <- eggnog_part %>%
  filter(pident > max_reasonable_pident)
write_tsv(eggnog_full_contams, file.path(output_dir, "full_contaminants.tsv"))
write_tsv(eggnog_part_contams, file.path(output_dir, "part_contaminants.tsv"))

eggnog_full_nocontams <- eggnog_full %>%
  filter(pident <= max_reasonable_pident) %>%
  separate_wider_delim(hum_protein, delim="|", names=c("sp", "Entry", "Entry Name"), cols_remove=FALSE)
eggnog_part_nocontams <- eggnog_part %>%
  filter(pident <= max_reasonable_pident) %>%
  separate_wider_delim(hum_protein, delim="|", names=c("sp", "Entry", "Entry Name"), cols_remove=FALSE)
eggnog_part_nocontams_prep <- eggnog_part_nocontams %>% select(-best_indiv_humcover, -best_indiv_humcover_right, -tot_length_right, -tot_humcover_right, -length_hum_right)
write_tsv(eggnog_full_nocontams, file.path(output_dir, "full_no_contaminants.tsv"))
write_ipc_file(eggnog_full_nocontams, file.path(output_dir, "full_no_contaminants.ipc"))

write_tsv(eggnog_part_nocontams_prep, file.path(output_dir, "part_no_contaminants.tsv"))
write_ipc_file(eggnog_part_nocontams_prep, file.path(output_dir, "part_no_contaminants.ipc"))

eggnog_fn <- eggnog_full_nocontams %>% group_by(hum_protein, Lineage) %>% nest()
eggnog_fc <- eggnog_fn %>% ungroup() %>% group_by(hum_protein) %>% count(name = "nFull")
write_tsv(eggnog_fc, file.path(output_dir, "eggnog_fc.tsv"))

# next use best_hits_only on full_no_contaminants.tsv to simplify
# and rerun humcover on part_no_contaminants.ipc
