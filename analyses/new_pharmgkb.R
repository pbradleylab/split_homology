# New PharmGKB analysis

library(tidyverse)
library(seqinr)
library(hgnc)
map <- purrr::map

#str_split_1(str, ", ") %>% str_split_fixed(., ":", 2) %>% as_tibble(colnames=c("key","value"))

hgnc_tsv <- read_tsv("../data/manual/hgnc_complete_set_2024-10-01.txt", guess_max=1e7, show_col_types=FALSE)
pharmgkb_genes <- read_tsv("../data/raw/pharmgkb/genes/genes.tsv")
pharmgkb_genes_crossref <- pharmgkb_genes %>%
  mutate(crossref = map(`Cross-references`, ~ {
    y <- str_split_1(.x, ", ") %>%
      str_split_fixed(":", 2)
    colnames(y) <- c("key", "value")
    as_tibble(y)
  }))
crossref_dbs <- pharmgkb_genes_crossref %>% unnest(cols=c(crossref)) %>% select(key) %>% distinct() %>% deframe()

crossref_hgnc <- pharmgkb_genes_crossref %>% inner_join(., hgnc_tsv, by=c("HGNC ID"="hgnc_id"))

# Deprecated - some will be NA if you do it this way because of old/outdated UniProt IDs
# pgkb_upr <- pharmgkb_genes_crossref %>% unnest(cols=c(crossref)) %>% filter(key=="UniProtKB") %>% inner_join(., eggnog_upr_full, by=c("value"="Entry"))

# This works better
pgkb_upr <- crossref_hgnc %>% filter(!is.na(uniprot_ids)) %>% separate_longer_delim(uniprot_ids, delim="|") %>% inner_join(., eggnog_upr_full, by=c("uniprot_ids"="Entry"))

all_pway_files <- list.files("../data/raw/pharmgkb/pathways-tsv/", "P.*Pathway.*tsv", full.names = TRUE)
concat_pways <- map(all_pway_files, ~ read_tsv(.x, col_types="ccccccccccc")) %>% bind_rows()

pharmgkb_chems <- read_tsv("../data/raw/pharmgkb/chemicals/chemicals.tsv")
pharmgkb_drugs <- pharmgkb_chems %>%
  filter(!grepl("Biological Intermediate", Type)) %>%
  mutate(chem_class = map(Type, ~ str_split_1(.x, ", "))) %>%
  unnest(chem_class) %>% 
  filter(chem_class %in% c("Drug", "Drug Class", "Prodrug", "Metabolite")) %>% # Note: metabolite here is drug metabolites, not biol intermediates
  nest(chem_class=chem_class)

pway_controllers <- concat_pways %>%
  filter(From != To) %>%
  filter(From %in% pharmgkb_drugs$Name | To %in% pharmgkb_drugs$Name) %>%
  filter(!is.na(Controller)) %>%
  filter(`Reaction Type` != "Transport") %>%
  mutate(Controller = gsub(", ", ",", Controller)) %>%
  separate_longer_delim(Controller, delim=",")

length(pway_controllers$Controller %>% unique) # note, only ~190

# mostly not current gene names, or classes of genes ("GSH")
missing_controllers <- setdiff(pway_controllers$Controller, pgkb_upr$Symbol)

pgkb_part_full <- inner_join(pway_controllers, pgkb_upr, by=c("Controller"="Symbol")) %>% select(From, To, Controller, Drugs, uniprot_ids, nPart, nFull)
write_tsv(pgkb_part_full, "pgkb_part_full.tsv")

#split_proteins <- pgkb_part_full %>% filter(nPart > nFull) %>% group_by(Controller, uniprot_ids, nPart, nFull) %>% nest() %>% arrange(-nPart) %>% mutate(From = map_chr(data, ~ paste0(unique(.x$From), collapse=";")))
split_drugs <- pgkb_part_full %>% filter(nPart > nFull) %>% select(Drugs) %>% filter(!is.na(Drugs)) %>% distinct() %>% deframe()

#full_proteins <- pgkb_part_full %>% filter(nPart < nFull) %>% group_by(Controller, uniprot_ids, nPart, nFull) %>% nest() %>% arrange(-nPart) %>% mutate(From = map_chr(data, ~ paste0(unique(.x$From), collapse=";")))
full_drugs <- pgkb_part_full %>% filter(nFull > nPart) %>% select(Drugs) %>% filter(!is.na(Drugs)) %>% distinct() %>% deframe() %>% map(., ~ str_split_1(.,", ")) %>% Reduce(c, .) %>% unique

full_metabs <- pgkb_part_full %>% filter(nFull > nPart) %>% select(From) %>% unlist() %>% unique() %>% map(., ~ str_split_1(.,", ")) %>% Reduce(c, .) %>% unique
all_metabs <- pgkb_part_full %>% select(From) %>% unlist() %>% unique() %>% gsub(",choline", ", choline",.) %>% map(., ~ str_split_1(.,", ")) %>% Reduce(c, .) %>% unique

pgkb_split <- pgkb_part_full %>% mutate(From=gsub(",choline", ", choline", From)) %>% mutate(From=map(From, ~ str_split_1(., ", "))) %>% unnest(From) %>% distinct()
total_nfullpart <- pgkb_split %>% select(From, uniprot_ids, nPart, nFull) %>% distinct() %>% group_by(From) %>% summarize(sum_nFull = sum(nFull), sum_nPart = sum(nPart), max_nFull = max(nFull), max_nPart = max(nPart))
full_proteins <- pgkb_split %>% filter(nPart < nFull) %>% group_by(Controller, uniprot_ids, nPart, nFull) %>% nest() %>% arrange(-nPart) %>% mutate(From = map_chr(data, ~ paste0(unique(.x$From), collapse=";")))
split_proteins <- pgkb_split %>% filter(nPart > nFull) %>% group_by(Controller, uniprot_ids, nPart, nFull) %>% nest() %>% arrange(-nPart) %>% mutate(From = map_chr(data, ~ paste0(unique(.x$From), collapse=";")))
all_proteins <- pgkb_split %>% group_by(Controller, uniprot_ids, nPart, nFull) %>% nest() %>% arrange(-nPart) %>% mutate(From = map_chr(data, ~ paste0(unique(.x$From), collapse=";")))


write_csv(split_proteins %>% select(-data), "Table3.csv")
write_csv(full_proteins %>% select(-data) %>% arrange(-nFull), "SuppTable3.csv")

pgkb_drugs <- pgkb_split %>% mutate(Drugs = map2_chr(Drugs, From, ~ { if (is.na(.x)) .y else .x }))



# What classes?

# fmn metabolism GO:0046444, fad metabolism GO:0046443, nicotinamide GO:0046496

# nucleobase processes
all_nucl_go <- GOBPOFFSPRING["GO:0006139"] %>% as.list() %>% .[[1]]
all_acyl_coa <- GOBPOFFSPRING["GO:0006637"] %>% as.list() %>% .[[1]]
all_coa_metab <- GOBPOFFSPRING["GO:0015936"] %>% as.list() %>% .[[1]]
all_fmn_metab <- c(GOBPOFFSPRING["GO:0046444"] %>% as.list() %>% .[[1]], "GO:0046444")
all_fad_metab <- c(GOBPOFFSPRING["GO:0046443"] %>% as.list() %>% .[[1]], "GO:0046443")
all_pyrid_metab <- c(GOBPOFFSPRING["GO:0019362"] %>% as.list() %>% .[[1]], "GO:0019362")

other_nucl <- Reduce(union, list(
  all_acyl_coa,
  all_coa_metab,
  all_fmn_metab,
  all_fad_metab,
  all_pyrid_metab,
  c("GO:0006637", "GO:0006763")
))

xeno_nearly_all_entries <- lapply(xeno_all, \(x) x$Entry)
existing_xeno <- Reduce(union, xeno_nearly_all_entries)

# What else could go here? Identifying a few other categories...
nucl_non_coa <- setdiff(union("GO:0006139", all_nucl_go),
                       other_nucl)
other_nucl_all <- upr_bp_unnest %>% filter(bp %in% other_nucl) %>% left_join(eggnog_upr_full) %>% select(-bp) %>% distinct()
nucl_all <- upr_bp_unnest %>% filter(bp %in% nucl_non_coa) %>% left_join(eggnog_upr_full) %>% select(-bp) %>% distinct() %>% filter(!(Entry %in% other_nucl_all$Entry))
all_redox_go <- c(GOMFOFFSPRING["GO:0016491"] %>% as.list() %>% .[[1]], "GO:0016491")
redox_all <- upr_mf_unnest %>% filter(mf %in% all_redox_go) %>% filter(`Entry Name` %in% mostly_full) %>% left_join(eggnog_upr_full) %>% select(-mf) %>% distinct()

# Label only genes that don't already have a class assigned
xeno_all_entries <- xeno_nearly_all_entries
xeno_all_entries[["redox"]] <- setdiff(redox_all$Entry, existing_xeno)
xeno_all_entries[["nucl"]] <- setdiff(nucl_all$Entry, union(existing_xeno, redox_all$Entry))
xeno_tbl <- xeno_all_entries %>% enframe(name="xeno_class", value="uniprot_ids") %>% unnest(cols=c(uniprot_ids))

full_xeno_cat <- left_join(full_proteins, xeno_tbl) %>%
  mutate(xeno_class = replace_na(xeno_class, "other"))

full_xeno_count <- full_xeno_cat %>%
  group_by(xeno_class) %>%
  count() %>%
  arrange(-n)

full_xeno_count %>% ungroup() %>% summarize(total=sum(n))
full_xeno_count %>% filter(!(xeno_class %in% c("nucl","redox","other"))) %>% ungroup() %>% summarize(known_total=sum(n))

split_xeno_cat <- left_join(split_proteins, xeno_tbl) %>%
  mutate(xeno_class = replace_na(xeno_class, "other"))

split_xeno_count <- split_xeno_cat %>%
  group_by(xeno_class) %>%
  count() %>%
  arrange(-n)

all_xeno_cat <- left_join(all_proteins, xeno_tbl) %>%
  mutate(xeno_class = replace_na(xeno_class, "other"))

all_xeno_count <- all_xeno_cat %>%
  group_by(xeno_class) %>%
  count() %>%
  arrange(-n)


####

pgkb_drugs <- pgkb_split %>% mutate(Drugs = map2_chr(Drugs, From, ~ { if (is.na(.x)) .y else .x })) %>% left_join(., select(all_xeno_cat, uniprot_ids, xeno_class))

pgkb_SH_drugs <- eggnog_part %>%
  filter(Entry %in% pgkb_drugs$uniprot_ids) %>%
  left_join(., pgkb_drugs, by=c("Entry"="uniprot_ids"),
            relationship = 'many-to-many') %>%
  select(-sp) %>%
  relocate(Drugs, .before=Entry) %>%
  arrange(Drugs, Lineage, g, featnum, strand, contig) %>%
  select(Drugs, From, To, Lineage, g, xeno_class, `Entry Name`, Entry, bac_protein, pident, length)
write_csv(pgkb_SH_drugs, "pgkb_SH_drugs.csv")


pgkb_FH_drugs <- full_best %>%
  filter(Entry %in% pgkb_drugs$uniprot_ids) %>%
  collect() %>%
  left_join(., pgkb_drugs, by=c("Entry"="uniprot_ids"),
            relationship = 'many-to-many') %>%
  select(-sp) %>%
  relocate(Drugs, .before=Entry) %>%
  arrange(Drugs, Lineage, nGenomes) %>%
  select(Drugs, From, To, Lineage, nGenomes, xeno_class, `Entry Name`, Entry, bac_protein, pident, length)

write_csv(pgkb_FH_drugs, "pgkb_FH_drugs.csv")
