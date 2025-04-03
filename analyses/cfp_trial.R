library(arrow)
library(tidyverse)
library(ape)
select <- dplyr::select
count <- dplyr::count

hc3_path <- "../data/processed/humcover3/"
prep_path <- "../data/processed/preprocess/"
compare_path <- "./param_compare"
dir.create(compare_path, showWarnings=FALSE)

params_to_test <- list.files(hc3_path)

uniprot_annotations <- read_tsv("../data/manual/uniprotkb_AND_model_organism_9606_2024_10_24.tsv.gz")
uniprot_annotations <- uniprot_annotations %>%
  mutate(mito = grepl("[Mm]itochondr", `Subcellular location [CC]`)) %>%
  filter(Reviewed == 'reviewed')

param_results <- map(params_to_test, \(.x) {
	eggnog_part <- read_ipc_file(file.path(hc3_path, .x, "part_humcover3.ipc"))
	eggnog_fc <- read_tsv(file.path(prep_path, .x, "eggnog_fc.tsv"))
	eggnog_pn <- eggnog_part %>% group_by(hum_protein, Lineage) %>% nest()
	eggnog_pc <- eggnog_pn %>% ungroup() %>% group_by(hum_protein) %>% count(name = "nPart")
	eggnog_all <- full_join(eggnog_fc, eggnog_pc) %>%
	  replace_na(list(nFull = 0, nPart = 0)) %>% 
	  separate(hum_protein, sep="\\|", into=c("sp", "Entry", "Entry Name"), remove=FALSE)
	eggnog_upr <- left_join(select(eggnog_all, -`Entry Name`), uniprot_annotations, by=c("Entry"))
	eggnog_missing <- setdiff(eggnog_all$Entry, uniprot_annotations$Entry)
	  print("Missing:")
	  print(eggnog_missing)
	eggnog_upr_full <- left_join(uniprot_annotations, select(eggnog_all, -`Entry Name`)) %>%
	  replace_na(list(nFull = 0, nPart = 0))
	any_split <- eggnog_upr_full %>% filter(nPart > 0) %>% select(Entry, `Entry Name`, `Protein names`, nFull, nPart) %>% mutate(`Protein names`=gsub("^([^\\[\\(]*) [\\[\\(].*", "\\1", `Protein names`)) %>% arrange(-nPart)
	write_csv(any_split,
	  file.path(compare_path, paste0(.x, "_any_split.csv")))
	any_full <- eggnog_upr_full %>% filter(nFull > 0) %>% select(Entry, `Entry Name`, `Protein names`, nFull, nPart) %>% mutate(`Protein names`=gsub("^([^\\[\\(]*) [\\[\\(].*", "\\1", `Protein names`)) %>% arrange(-nFull)
	write_csv(any_full,
	  file.path(compare_path, paste0(.x, "_any_full.csv")))
	return(list(split=any_split, full=any_full, overall=eggnog_upr_full))
})

names(param_results) <- params_to_test
all_split <- Reduce(bind_rows, imap(param_results, ~ mutate(.x$split, type="split", paramset=.y)))
all_full <- Reduce(bind_rows, imap(param_results, ~ mutate(.x$full, type="full", paramset=.y)))
overall_results <- bind_rows(all_split, all_full)
write_csv(overall_results, file.path(compare_path, "overall_results.csv"))

split_compare_f <- overall_results %>% filter(type=="split") %>% arrange(-nPart) %>% group_by(paramset) %>% select(paramset, Entry, `Entry Name`, nPart) %>% arrange(paramset) %>%  mutate(paramset=gsub("HumanUPR_(.*)_src_20000_(.*)", "p\\1_\\2", paramset)) %>% pivot_wider(names_from=paramset, values_from=nPart, values_fill=0)
full_compare_f <- overall_results %>% filter(type=="full") %>% arrange(-nFull) %>% group_by(paramset) %>% select(paramset, Entry, `Entry Name`, nFull) %>% arrange(paramset) %>%  mutate(paramset=gsub("HumanUPR_(.*)_src_20000_(.*)", "p\\1_\\2", paramset)) %>% pivot_wider(names_from=paramset, values_from=nFull, values_fill=0)
write_tsv(split_compare_f, "split_compare_f.tsv")
write_tsv(full_compare_f, "full_compare_f.tsv")

any_split_greater <- overall_results %>% filter(type=="split") %>% arrange(-nPart) %>% filter(nPart > nFull) %>% select(Entry) %>% distinct() %>% deframe()
any_full_greater <- overall_results %>% filter(type=="full") %>% arrange(-nFull) %>% filter(nFull > nPart) %>% select(Entry) %>% distinct() %>% deframe()
split_compare_2f <- overall_results %>% filter(type=="split", Entry %in% any_split_greater) %>% arrange(-nPart) %>% group_by(paramset) %>% select(paramset, Entry, `Entry Name`, nPart) %>% arrange(paramset) %>%  mutate(paramset=gsub("HumanUPR_(.*)_src_20000_(.*)", "p\\1_\\2", paramset)) %>% pivot_wider(names_from=paramset, values_from=nPart, values_fill=0)
full_compare_2f <- overall_results %>% filter(type=="full", Entry %in% any_full_greater) %>% arrange(-nFull) %>% group_by(paramset) %>% select(paramset, Entry, `Entry Name`, nFull) %>% arrange(paramset) %>%  mutate(paramset=gsub("HumanUPR_(.*)_src_20000_(.*)", "p\\1_\\2", paramset)) %>% pivot_wider(names_from=paramset, values_from=nFull, values_fill=0)

split_sensitivity_tbl <- split_compare_2f %>% pivot_longer(p0.5_60:p0.75_80) %>% mutate(found=value > 0) %>% group_by(Entry, `Entry Name`) %>% left_join(., summarize(., n_params_detected = sum(found), median_found = median(value))) %>% select(-found) %>% pivot_wider() %>% arrange(-n_params_detected, -median_found) %>% mutate(found_in_ours = p0.67_70>0)
write_tsv(split_sensitivity_tbl, "split_sensitivity_tbl.tsv")

full_sensitivity_tbl <- full_compare_2f %>% pivot_longer(p0.5_60:p0.75_80) %>% mutate(found=value > 0) %>% group_by(Entry, `Entry Name`) %>% left_join(., summarize(., n_params_detected = sum(found), median_found = median(value))) %>% select(-found) %>% pivot_wider() %>% arrange(-n_params_detected, -median_found) %>% mutate(found_in_ours = p0.67_70>0)
write_tsv(full_sensitivity_tbl, "full_sensitivity_tbl.tsv")
n_full_total <- full_sensitivity_tbl %>% nrow
mini_full_sens_tbl <- full_sensitivity_tbl %>% group_by(n_params_detected, found_in_ours) %>% count() %>% arrange(-n_params_detected, -found_in_ours) %>% mutate(tot=n_full_total) %>% mutate(frac=n/tot)
write_tsv(mini_full_sens_tbl, "mini_full_sens_tbl.tsv")

cum_frac_full <- full_sensitivity_tbl %>% filter(found_in_ours) %>% group_by(n_params_detected) %>% count() %>% arrange(-n_params_detected) %>% ungroup() %>% mutate(run=cumsum(n)) %>% mutate(frac=n/2569) %>% mutate(runfrac = cumsum(frac))
write_tsv(cum_frac_full, "cum_frac_full.tsv")

split_cmp <- overall_results %>% select(-type) %>% distinct() %>% mutate(paramset=gsub("HumanUPR_([^_]+)_src_20000_(.*)", "B\\1_H\\2", paramset)) %>% rename(split=nPart,full=nFull) %>% pivot_wider(values_from=c(full,split), names_from=paramset, values_fill=0, names_vary="slowest") %>% rowwise() %>% mutate(full_median=median(c_across(starts_with("full"))), split_median=median(c_across(starts_with("split")))) %>% relocate(Entry, `Entry Name`, `Protein names`, full_median, split_median) %>% arrange(-split_median, -full_median)
write_csv(split_cmp, "sensitivity_analysis_table.csv")

full_comparison <- overall_results %>% filter(nFull > nPart) %>% select(-type) %>% distinct() %>% group_by(Entry, `Entry Name`) %>% summarize(n_paramsets=n(), median_nFull = median(nFull), median_nPart =median(nPart), found="HumanUPR_0.67_src_20000_70" %in% paramset) %>% ungroup()
full_comparison_2 <- left_join(full_comparison, full_compare_f, by=c("Entry", "Entry Name"))%>% arrange(desc(found), desc(n_paramsets), desc(median_nFull))
write_csv(full_comparison_2, "SuppTable2.csv")
full_summary <- full_comparison %>% group_by(found, n_paramsets) %>% count
write_csv(full_summary, "full_summary.csv")

split_comparison <- overall_results %>% filter(nFull < nPart) %>% select(-type) %>% distinct() %>% group_by(Entry, `Entry Name`) %>% summarize(n_paramsets=n(), median_nFull = median(nFull), median_nPart =median(nPart), found="HumanUPR_0.67_src_20000_70" %in% paramset) %>% ungroup()
split_comparison_2 <- left_join(split_comparison, split_compare_f, by=c("Entry", "Entry Name")) %>% arrange(desc(found), desc(n_paramsets), desc(median_nPart))
write_csv(split_comparison, "SuppTable1.csv")
split_summary <- split_comparison %>% group_by(found, n_paramsets) %>% count
write_csv(split_summary, "split_summary.csv")
