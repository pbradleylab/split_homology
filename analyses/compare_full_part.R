library(arrow)
library(tidyverse)
library(topGO)
library(readxl)
library(ape)
library(ggtree)
library(picante)
library(phytools)
library(ggtreeExtra)
library(ggnewscale)
select <- dplyr::select
count <- dplyr::count

genomes_md <- read_tsv("../data/raw/genomes-all_metadata.tsv")
genome_tree <- read.tree("../data/raw/bac120_iqtree.nwk")

param_list <- list(c(0.5,60), c(0.67,70), c(0.75,80))
#param_list <- list(c(0.67,70))

for (params in param_list) {
  print(paste0("bac_length: ", params[1], "; hum_length: ", params[2]))
  bac_length = params[1]
  hum_length = params[2]
  
  eggnog_part <- read_ipc_file(
    file.path("../data/processed/humcover3/",
              paste0("HumanUPR_", bac_length, "_src_20000_", hum_length),
              "part_humcover3.ipc")
  )
  eggnog_fc <- read_tsv(
    file.path("../data/processed/preprocess/",
              paste0("HumanUPR_", bac_length, "_src_20000_", hum_length),
              "eggnog_fc.tsv")
  )
  full_best <- open_dataset(
      file.path("../data/processed/humcover3/",
                paste0("HumanUPR_", bac_length, "_src_20000_", hum_length),
                "full_humcover3.ipc"),
      format="ipc"
  )
  output_dir <- file.path("./output/", paste0("B", bac_length), paste0("H", hum_length))
  dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
  
  eggnog_pn <- eggnog_part %>% group_by(hum_protein, Lineage) %>% nest()
  eggnog_pc <- eggnog_pn %>% ungroup() %>% group_by(hum_protein) %>% count(name = "nPart")
  
  
  uniprot_annotations <- read_tsv("../data/manual/uniprotkb_AND_model_organism_9606_2024_10_24.tsv.gz")
  uniprot_annotations <- uniprot_annotations %>%
    mutate(mito = grepl("[Mm]itochondr", `Subcellular location [CC]`)) %>%
    filter(Reviewed == 'reviewed')
  
  eggnog_all <- full_join(eggnog_fc, eggnog_pc) %>%
    replace_na(list(nFull = 0, nPart = 0)) %>% 
    separate(hum_protein, sep="\\|", into=c("sp", "Entry", "Entry Name"), remove=FALSE)
  
  eggnog_upr <- left_join(select(eggnog_all, -`Entry Name`), uniprot_annotations, by=c("Entry"))
  eggnog_missing <- setdiff(eggnog_all$Entry, uniprot_annotations$Entry)
  print("Missing:")
  print(eggnog_missing)
  
  eggnog_upr_full_na <- full_join(uniprot_annotations, eggnog_all) %>%
    replace_na(list(nFull = 0, nPart = 0))
  eggnog_upr_full <- eggnog_upr_full_na %>% filter(!is.na(mito))
  
  table1 <-   eggnog_upr_full %>% filter(nPart > nFull) %>% select(Entry, `Entry Name`, `Protein names`, nFull, nPart) %>% mutate(`Protein names`=gsub("^([^\\[\\(]*) [\\[\\(].*", "\\1", `Protein names`)) %>% arrange(-nPart)
  write_csv(table1, file.path(output_dir, "Table1.csv"))
  
  # Mitochondrial localization
  
  full_contingency <- eggnog_upr_full %>%
    mutate(full_only = coalesce(((nFull > 0) & (nPart == 0)), FALSE)) %>%
    group_by(full_only, mito) %>%
    count() %>%
    pivot_wider(names_from="mito", values_from="n") %>%
    ungroup() %>%
    select(-full_only) %>%
    data.matrix
  fisher.test(full_contingency)
  # p < 2.2e-16, OR = 4.45
  
  part_contingency <- eggnog_upr_full %>%
    mutate(part_only = coalesce(((nPart > 0) & (nFull == 0)), FALSE)) %>%
    group_by(part_only, mito) %>%
    count() %>%
    pivot_wider(names_from="mito", values_from="n") %>%
    ungroup() %>%
    select(-part_only) %>%
    data.matrix
  fisher.test(part_contingency)
  # p = 0.29, OR = 1.97
  
  full_only_fracs <- eggnog_upr_full %>%
    filter(nPart == 0) %>%
    arrange(-nFull) %>%
    mutate(num_prots = cumsum(row_number())) %>%
    mutate(frac = cummean(mito)) %>% 
    select(nFull, frac, num_prots)%>% 
    distinct() %>%
    mutate(which="full") %>%
    dplyr::rename(n = nFull)
  
  full_only_genes <- eggnog_upr_full %>%
    filter(nPart == 0, nFull > 0) %>%
    arrange(-nFull) %>%
    select(`Entry Name`) %>%
    deframe
  
  part_only_genes <- eggnog_upr_full %>%
    filter(nFull == 0, nPart > 0) %>%
    arrange(-nFull) %>%
    select(`Entry Name`) %>%
    deframe
  
  has_part_genes <- eggnog_upr_full %>%
    filter(nPart > 0) %>%
    arrange(-nFull) %>%
    select(`Entry Name`) %>%
    deframe
  
  part_only_fracs <- eggnog_upr_full %>%
    filter(nFull == 0) %>%
    arrange(-nPart) %>%
    mutate(num_prots = cumsum(row_number())) %>%
    mutate(frac = cummean(mito)) %>% 
    select(nPart, frac, num_prots) %>% 
    distinct() %>%
    mutate(which="part") %>%
    dplyr::rename(n = nPart)
  
  all_fracs <- bind_rows(full_only_fracs, part_only_fracs)
  
  # Generate mitochondrial enrichment supplementary figure
  pdf(file.path(output_dir, "supplementary-mitochondrial-enr.pdf"))
  ggplot(all_fracs, aes(x=n, y=frac, color=which)) +
    geom_line(alpha=1, lwd=0.5) +
    geom_abline(slope=0, intercept = mean(eggnog_upr_full$mito), lty=2) +
    geom_point(aes(size=log10(num_prots))) +
    scale_size_continuous(range=c(0.3,3)) +
    theme_minimal()
  dev.off()
  
  # n=16 only have splits, not full-length, while 23 have more splits than full-length
  print(eggnog_upr_full %>% filter(nPart > 0, nFull == 0) %>% arrange(-nPart))
  print(eggnog_upr_full %>% filter(nPart > nFull) %>% arrange(-nPart))
  
  # Enrichment analysis
  select <- dplyr::select
  go_term_tbl <- enframe(unlist(Term(GOTERM)))
  
  map <- purrr::map
  upr_mf_unnest <- uniprot_annotations %>%
    mutate(mf = str_extract_all(`Gene Ontology (molecular function)`, "\\[(GO:.......)\\]")) %>%
    mutate(mf = map(mf, ~ gsub("\\[(.*)\\]", "\\1", .))) %>%
    select(`Entry Name`, mf) %>%
    unnest(mf) %>%
    left_join(., select(eggnog_upr_full, `Entry Name`, nFull, nPart)) %>%
    mutate(full_only = `Entry Name` %in% full_only_genes) %>%
    mutate(part_only = `Entry Name` %in% part_only_genes)
  
  upr_bp <- uniprot_annotations %>%
    mutate(bp = str_extract_all(`Gene Ontology (biological process)`, "\\[(GO:.......)\\]")) %>%
    mutate(bp = map(bp, ~ gsub("\\[(.*)\\]", "\\1", .)))
  upr_bp_unnest <- upr_bp %>%
    select(`Entry Name`, bp) %>%
    unnest(bp) %>%
    left_join(., select(eggnog_upr_full, `Entry Name`, nFull, nPart)) %>%
    mutate(full_only = `Entry Name` %in% full_only_genes) %>%
    mutate(part_only = `Entry Name` %in% part_only_genes)
  upr_bp_GO_mapping <- upr_bp_unnest %>%
    group_by(`Entry Name`) %>%
    nest() %>%
    mutate(data = map(data, ~ c(.$bp))) %>%
    deframe
  
  
  
  upr_mf_GO_mapping <- upr_mf_unnest %>%
    group_by(`Entry Name`) %>%
    nest() %>%
    mutate(data = map(data, ~ c(.$mf))) %>%
    deframe
  
  # Wrapper to run TopGO based on the information we typically actually have
  get_go_enrichment <- function(sig_genes,
                                all_genes,
                                desc="enrichments",
                                enr_sig=0.25,
                                write_to=NULL,
                                mapping=upr_bp_GO_mapping,
                                go_tbl=go_term_tbl,
                                min_nodesize=5,
                                adj_method='BH',
                                onto="BP") {
    sig_fac <- factor(as.integer(all_genes %in% sig_genes))
    names(sig_fac) <- all_genes
    topgo_data <- new("topGOdata",
                      description = "full_only enrichments",
                      ontology = onto,
                      annotationFun = annFUN.gene2GO,
                      gene2GO = mapping,
                      allGenes = sig_fac,
                      nodeSize = min_nodesize
    )
    topgo_results <- runTest(topgo_data, algorithm="weight01", statistic="fisher")
    topgo_pv_corr <- enframe(p.adjust(score(topgo_results), adj_method), value = "corr_pval")
    topgo_sig <- filter(topgo_pv_corr, corr_pval <= 0.05)$name
    topgo_tbl <- left_join(topgo_pv_corr, go_tbl, by="name")
    if (!(is.null(write_to))) {
      write_csv(topgo_tbl %>% filter(corr_pval <= enr_sig), file = write_to)
    }
    return(list(data=topgo_data, results=topgo_results, tbl=topgo_tbl))
  }
  
  # actually conduct enrichments
  mostly_full <- eggnog_upr_full %>% filter(nPart < nFull) %>% select(`Entry Name`) %>% deframe
  mostly_split <- eggnog_upr_full %>% filter(nPart > nFull) %>% select(`Entry Name`) %>% deframe
  all_genes <- eggnog_upr_full$`Entry Name` %>% unique
  all_ortho_genes <- eggnog_upr_full %>% filter(nPart > 0 | nFull > 0) %>% select(`Entry Name`) %>% distinct() %>% deframe
  
  mostly_full_no_mito <- eggnog_upr_full %>% filter(!mito, nPart < nFull) %>% select(`Entry Name`) %>% deframe
  all_genes_no_mito <- eggnog_upr_full %>% filter(!mito) %>% select(`Entry Name`) %>% deframe %>% unique
  
  # make sure none are exactly equal
  equal.nonzero <- eggnog_upr_full %>% filter(nPart > 0) %>% filter(nPart == nFull) %>% select(`Entry Name`) %>% deframe
  print(length(equal.nonzero) == 0)
  
  go_gaf <- read_tsv("../data/manual/goa_human.gaf.gz", comment="!", col_names=FALSE)
  noIEA_annots <- left_join(go_gaf, uniprot_annotations, by=c("X2"="Entry")) %>% filter(!grepl("NOT", X4), X7 != "IEA") %>% select(`Entry Name`, X5) %>% inner_join(., upr_bp_unnest,  by=c("Entry Name", "X5"="bp")) %>% distinct() %>% group_by(`Entry Name`) %>% nest %>%  mutate(data = map(data, ~ c(.$X5))) %>%deframe
  vstringent_annots <- left_join(go_gaf, uniprot_annotations, by=c("X2"="Entry")) %>% filter(!grepl("NOT", X4), (X7 %in% c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP"))) %>% select(`Entry Name`, X5) %>% inner_join(., upr_bp_unnest,  by=c("Entry Name", "X5"="bp")) %>% distinct() %>% group_by(`Entry Name`) %>% nest %>%  mutate(data = map(data, ~ c(.$X5))) %>%deframe
  
  mostly_full_vs_everything <- get_go_enrichment(
    mostly_full,
    all_genes,
    desc="Proteins with mostly full-length homologs vs. UniProt",
    write_to=file.path(output_dir, "mostly_full_vs_everything_GO.csv"))
  mostly_full_vs_everything_noIEA <- get_go_enrichment(
    mostly_full,
    intersect(all_genes, names(noIEA_annots)),
    mapping=noIEA_annots,
    desc="Proteins with mostly full-length homologs vs. UniProt, no IEA",
    write_to=file.path(output_dir, "mostly_full_vs_everything_GO_noIEA.csv"))
  mostly_full_vs_everything_vs <- get_go_enrichment(
    mostly_full,
    intersect(all_genes, names(vstringent_annots)),
    mapping=vstringent_annots,
    desc="Proteins with mostly full-length homologs vs. UniProt, very stringent",
    write_to=file.path(output_dir, "mostly_full_vs_everything_GO_vs.csv"))
  mostly_full_vs_everything_no_mito <- get_go_enrichment(
    mostly_full_no_mito,
    all_genes_no_mito,
    desc="Proteins with mostly full-length homologs vs. UniProt",
    write_to=file.path(output_dir, "mostly_full_vs_everything_GO.csv"))
  mostly_split_vs_everything <- get_go_enrichment(
    mostly_split,
    all_genes,
    desc="Proteins with mostly split homologs vs. UniProt",
    write_to=file.path(output_dir, "mostly_split_vs_everything_GO.csv"))
  mostly_split_vs_everything_noIEA <- get_go_enrichment(
    mostly_split,
    intersect(all_genes, names(noIEA_annots)),
    mapping=noIEA_annots,
    desc="Proteins with mostly split homologs vs. UniProt, no IEA",
    write_to=file.path(output_dir, "mostly_split_vs_everything_GO_noIEA.csv"))
  mostly_split_vs_everything_vs <- get_go_enrichment(
    mostly_split,
    intersect(all_genes, names(vstringent_annots)),
    mapping=vstringent_annots,
    desc="Proteins with mostly split homologs vs. UniProt, very stringent",
    write_to=file.path(output_dir, "mostly_split_vs_everything_GO_vs.csv"))
  
  mfe_sig <- mostly_full_vs_everything$tbl %>% arrange(corr_pval) %>% filter(corr_pval <= 0.05)
  mfe_nm_sig <- mostly_full_vs_everything_no_mito$tbl %>% arrange(corr_pval) %>% filter(corr_pval <= 0.05)
  msplit_sig <- mostly_split_vs_everything$tbl %>% arrange(corr_pval) %>% filter(corr_pval <= 0.05)
  
  mfe_noIEA_sig <- mostly_full_vs_everything_noIEA$tbl %>% arrange(corr_pval) %>% filter(corr_pval <= 0.05)
  mfe_vs_sig <- mostly_full_vs_everything_vs$tbl %>% arrange(corr_pval) %>% filter(corr_pval <= 0.05)
  msplit_noIEA_sig <- mostly_split_vs_everything_noIEA$tbl %>% arrange(corr_pval) %>% filter(corr_pval <= 0.05)
  msplit_vs_sig <- mostly_split_vs_everything_vs$tbl %>% arrange(corr_pval) %>% filter(corr_pval <= 0.05)
  
  
  supptbl_full_sig <- bind_rows(
    imap(list(all=mostly_full_vs_everything,
              noIEA=mostly_full_vs_everything_noIEA,
              onlyExp=mostly_full_vs_everything_vs),
       \(x, n) mutate(x$tbl, annotations=n))) %>%
    group_by(name, value) %>%
    filter(any(corr_pval <= 0.05)) %>%
    pivot_wider(values_from=corr_pval, names_from=annotations) %>%
    arrange(all)
  
  supptbl_split_sig <- bind_rows(
    imap(list(all=mostly_split_vs_everything,
              noIEA=mostly_split_vs_everything_noIEA,
              onlyExp=mostly_split_vs_everything_vs),
         \(x, n) mutate(x$tbl, annotations=n))) %>%
    group_by(name, value) %>%
    filter(any(corr_pval <= 0.05)) %>%
    pivot_wider(values_from=corr_pval, names_from=annotations) %>%
    arrange(all)
  
  write_tsv(mfe_sig,
            file.path(output_dir, "tbl-full-vs-all-enr.tsv"))
  write_tsv(mfe_nm_sig,
            file.path(output_dir, "tbl-full-vs-all-no-mito-enr.tsv"))
  write_tsv(msplit_sig,
            file.path(output_dir, "tbl-split-vs-all-enr.tsv"))
  write_tsv(supptbl_full_sig,
            file.path(output_dir, "SuppTable3.tsv"))
  write_tsv(supptbl_split_sig,
            file.path(output_dir, "SuppTable5.tsv"))
  
  split_go_results <- upr_bp_unnest %>% filter(bp %in% msplit_sig$name, `Entry Name` %in% mostly_split)
  
  # categories + tree
  
  gt_rooted <- midpoint_root(genome_tree)
  
  chart_cols <- c(salmon = "#f9977b",      # 1
                  oceanblue = "#4d779e",   # 2
                  green = "#007600",       # 3
                  orange2 = "#ee9a00",     # 4
                  tomato = "#ff6347",      # 1
                  steelblue2 = "#4f94cd",  # 2
                  lime = "#28e158",        # 3
                  peach = "#e3bb65",       # 4
                  skyblue = "#18c9ef",     # 5
                  violetred = "#e889b4",   # 6
                  tan = "#a8552c",         # 7
                  puce = "#9e466e",        # 6
                  umber = "#6e2d0d",       # 7
                  darkgray = "#444444",    # 8
                  pewter = "#a4a4a4",      # 8
                  lightgray = "#cccccc",   # 8
                  black = "#000000")       # 8
  
  
  get_pd_row <- function(h_protein, full=full_best, md=genomes_md, phy=gt_rooted, rowname=NULL) {
    these_proteins <- full %>% filter(hum_protein %in% h_protein) %>% select(bac_protein) %>% collect() %>% deframe()
    these_genomes <- gsub("(GUT_GENOME......)_.....", "\\1", these_proteins) %>% unique()
    these_species <- md %>% filter(Genome %in% these_genomes) %>% select(Species_rep) %>% distinct() %>% deframe()
    TipLabels <- phy$tip.label
    PresAbs <- as.numeric(TipLabels %in% these_species)
    matrix(nr=1, nc=length(PresAbs), data = PresAbs, dimnames=list(rowname, TipLabels))
  }
  
  eggnog_upr_mostly_full <- eggnog_upr_full %>% filter(nFull > nPart)
  # nucleobase processes
  all_xeno_go <- c(GOBPOFFSPRING["GO:0006805"] %>% as.list() %>% .[[1]], "GO:0006805")
  all_nucl_go <- GOBPOFFSPRING["GO:0006139"] %>% as.list() %>% .[[1]]
  all_acyl_coa <- GOBPOFFSPRING["GO:0006637"] %>% as.list() %>% .[[1]]
  all_coa_metab <- GOBPOFFSPRING["GO:0015936"] %>% as.list() %>% .[[1]]
  nucl_non_coa <- setdiff(c("GO:0006139", all_nucl_go), union(all_acyl_coa, union(all_coa_metab, c("GO:0006637", "GO:0006763"))))
  nucl_all <- upr_bp_unnest %>% filter(bp %in% nucl_non_coa) %>% left_join(eggnog_upr_full) %>% select(-bp) %>% distinct()
  nucl_mf_all <- nucl_all %>% filter(`Entry Name` %in% mostly_full)
  
  all_akr_go <- c(GOMFOFFSPRING["GO:0004033"] %>% as.list() %>% .[[1]], "GO:0004033")
  akr_all <- upr_mf_unnest %>% filter(mf %in% all_akr_go) %>% filter(`Entry Name` %in% mostly_full) %>% left_join(eggnog_upr_full) %>% select(-mf) %>% distinct()
  akr_mf_all <- akr_all %>%  filter(`Entry Name` %in% mostly_full)
  
  xeno_mf <- upr_bp_unnest %>%
    filter(bp %in% all_xeno_go) %>%
    filter(`Entry Name` %in% mostly_full) %>%
    left_join(eggnog_upr_full) 
  xeno_mp <- upr_bp_unnest %>%
    filter(bp %in% all_xeno_go) %>%
    filter(`Entry Name` %in% mostly_split) %>%
    left_join(eggnog_upr_full) 
  
  xeno_all_AK <- eggnog_upr_full %>% filter(grepl("Aldo/keto reductase family", `Protein families`))
  xeno_all_UDP <- eggnog_upr_full %>% filter(grepl("UDP-glycosyltransferase family", `Protein families`))
  xeno_all_GST <- eggnog_upr_full %>% filter(grepl("GST superfamily", `Protein families`)) 
  #xeno_all_aryl <- eggnog_upr_full %>% filter(grepl("Arylamine N-acetyltransferase family|\'GDXG\' lipolytic enzyme family", `Protein families`))
  xeno_all_aryl <- eggnog_upr_full %>% filter(grepl("Arylamine N-acetyltransferase family", `Protein families`))
  xeno_all_gdxg <- eggnog_upr_full %>% filter(grepl("\'GDXG\' lipolytic enzyme family", `Protein families`))
  
  xeno_all_cyto <- eggnog_upr_full %>% filter(grepl("Cytochrome P450 family", `Protein families`))# %>% filter(hum_protein %in% xeno$hum_protein)
  #xeno_all_ester <- eggnog_upr_full %>% filter(grepl("[Ee]sterase", `Protein names`))
  xeno_all_ester <- eggnog_upr_full %>% filter(grepl("Type-B carboxylesterase/lipase family", `Protein families`))
  xeno_all_flavin <- eggnog_upr_full %>% filter(grepl("Flavin monoamine oxidase family|FMO family", `Protein families`))
  xeno_all_sdr <- eggnog_upr_full %>% filter(grepl("Short-chain dehydrogenases/reductases \\(SDR\\)", `Protein families`))
  xeno_all_quin <- eggnog_upr_full %>% filter(grepl("Quinone oxidoreductase subfamily", `Protein families`))
  xeno_all_aldh <- eggnog_upr_full %>% filter(grepl("Aldehyde dehydrogenase family", `Protein families`))
  xeno_all <- list(aldo=xeno_all_AK,
                   udp=xeno_all_UDP,
                   gst=xeno_all_GST,
                   ester=xeno_all_ester,
                   aryl=xeno_all_aryl,
                   gdxg=xeno_all_gdxg,
                   cyto=xeno_all_cyto,
                   flavin=xeno_all_flavin,
                   aldh = xeno_all_aldh,
                   quin = xeno_all_quin,
                   sdr = xeno_all_sdr
  )
  xeno_all_hp <- lapply(xeno_all, \(x) x$hum_protein)
  xeno_mf_all_hp <- lapply(xeno_all_hp, \(x) intersect(x, eggnog_upr_mostly_full$hum_protein))
  xeno_mf_tbl <- xeno_mf_all_hp %>% enframe(name="xeno_class", value="hum_protein") %>% unnest(cols=c(hum_protein))
  xeno_st <- Reduce(rbind, map2(xeno_mf_all_hp, names(xeno_all_hp), ~ get_pd_row(.x, rowname=.y)))
  xeno_st_tbl <- xeno_st %>% t %>% as_tibble(rownames="tip")
  xeno_fac_tbl <- xeno_st_tbl %>% mutate(across(aldo:sdr, \(x) factor(x))) %>%
    mutate(across(aldo:sdr, \(x) { levels(x) <- c("N", "Y"); x}))
  xeno_na_tbl <- xeno_fac_tbl %>% mutate(across(aldo:sdr, ~ na_if(.x, "N")))
  xeno_st_pd <- picante::pd(xeno_st, gt_rooted) %>% arrange(-PD) %>% filter(PD > 0)
  
  
  # n_bacprot here is how many distinct uhgp-90 bacterial proteins detected in each family, n_prot is how many distinct human proteins had homologs in that family
  family_xenos <- lapply(xeno_mf_all_hp[rownames(xeno_st_pd)],  # this re-orders by descending PD
                         function(x) {
                           full_best %>%
                             filter(hum_protein %in% x) %>%
                             collect() %>% 
                             separate_wider_delim(Lineage, delim=";", names=c("domain","phylum","class","order","family","genus","species")) %>%
                             select(hum_protein, domain:species, bac_protein) %>%
                             group_by(family) %>%
                             summarize(n_sp = length(unique(species)), n_bacprot = length(unique(bac_protein)), n_prot = length(unique(hum_protein))) %>%
                             arrange(-n_sp)
                         })
  n_found_per_fam <- map2(family_xenos, names(family_xenos), ~ mutate(.x, xeno_class=.y)) %>% Reduce(bind_rows, .)
  top_3_per_fam <- n_found_per_fam %>% ungroup() %>% group_by(xeno_class) %>% slice_max(n_sp, n=3)
  top_3_per_fam_bacprot <- n_found_per_fam %>% ungroup() %>% group_by(xeno_class) %>% slice_max(n_bacprot, n=3)
  all_top_3_fams <- union(top_3_per_fam$family, top_3_per_fam_bacprot$family)
  
  species_md <- genomes_md %>% select(Species_rep, Lineage) %>% distinct() %>% separate_wider_delim(Lineage, delim=";", names=c("domain","phylum","class","order","family","genus","species"))
  
  
  table2a <- n_found_per_fam %>% select(-n_bacprot, -n_prot) %>% filter(family %in% all_top_3_fams) %>%
    pivot_wider(names_from=xeno_class, values_from=n_sp, values_fill = 0, names_prefix = "sp_") %>%
    inner_join(., select(species_md, phylum:family) %>% distinct(), by="family") %>%
    select(-class, -order) %>% 
    relocate(phylum, .before="family")
  table2b <- n_found_per_fam %>% select(-n_sp, -n_prot) %>%
    filter(family %in% all_top_3_fams) %>%
    pivot_wider(names_from=xeno_class, values_from=n_bacprot, values_fill = 0, names_prefix = "bp_") %>%
    inner_join(., select(species_md, phylum:family) %>% distinct(), by="family") %>%
    select(-class, -order) %>% 
    relocate(phylum, .before="family")
  table2 <- left_join(table2a, table2b, by=c("phylum", "family")) %>%
    mutate(phylum = gsub("p__", "", phylum)) %>%
    mutate(family = gsub("f__", "", family))
  
  highlight_tbl <- table2 %>% select(-phylum, -family) %>% mutate(across(everything(), ~ rank(-.x) <= 3)) %>% mutate(i = row_number()) %>% pivot_longer(-i) %>% filter(value) %>% mutate(j = map_int(name, ~ which(colnames(table2)==.x))) %>% arrange(j, i) %>% select(-value, -name)
  highlight_tbla <- table2a %>% select(-phylum, -family) %>% mutate(across(everything(), ~ rank(-.x) <= 3)) %>% mutate(i = row_number()) %>% pivot_longer(-i) %>% filter(value) %>% mutate(j = map_int(name, ~ which(colnames(table2a)==.x))) %>% arrange(j, i) %>% select(-value, -name)
  highlight_tblb <- table2b %>% select(-phylum, -family) %>% mutate(across(everything(), ~ rank(-.x) <= 3)) %>% mutate(i = row_number()) %>% pivot_longer(-i) %>% filter(value) %>% mutate(j = map_int(name, ~ which(colnames(table2b)==.x))) %>% arrange(j, i) %>% select(-value, -name)
  
  pretty_print_table2 <- function(tbl) {
    highlight_tbl <- tbl %>% select(-phylum, -family) %>%
      mutate(across(everything(), ~ rank(-.x) <= 3)) %>%
      mutate(i = row_number()) %>%
      pivot_longer(-i) %>%
      filter(value) %>%
      mutate(j = map_int(name, ~ which(colnames(tbl)==.x))) %>%
      arrange(j, i) %>%
      select(-value, -name)
    tbl <- tbl %>% mutate(phylum = gsub("p__", "", phylum)) %>%
      mutate(family = gsub("f__", "", family))
    flext2 <- flextable(tbl) %>%
      font(part = "all", fontname = "Aptos Narrow") %>%
      set_header_labels(values=setNames(gsub("[^_]+_", "",
                                             colnames(tbl)[3:ncol(tbl)]),
                                        colnames(tbl)[3:ncol(tbl)]))
    for (r in 1:nrow(highlight_tbl)) {
      flext2 <- flext2 %>% bg(i=highlight_tbl[r,]$i, j=highlight_tbl[r,]$j, "#D0D0D0")
    }
    # flext2 <- flext2 %>%
    #   add_header_row(values=c("", "Species", "Bacterial proteins"), colwidths=c(2,10,10))
    flext2 %>% padding(padding.left = 0, padding.right = 0, padding.top = 2, padding.bottom = 2) %>% line_spacing(space=0.9) %>% autofit(add_w = 0, add_h=0)
  }
  
  library(flexlsx)
  library(flextable)
  wb <- openxlsx2::wb_workbook()$add_worksheet("xenobiotics_sp")$add_worksheet("xenobiotics_prot")
  wb <- wb_add_flextable(wb, "xenobiotics_sp", pretty_print_table2(table2a), dims = "A1")
  wb <- wb_add_flextable(wb, "xenobiotics_prot", pretty_print_table2(table2b), dims = "A1")
  wb$save(file.path(output_dir, "Table2_Redo.xlsx"))
  save_as_docx(pretty_print_table2(table2a) %>% fit_to_width(max_width=7.25, unit="in"), path = file.path(output_dir, "Table2.docx"))
  save_as_docx(pretty_print_table2(table2b) %>% fit_to_width(max_width=7.25, unit="in"), path = file.path(output_dir, "Table3.docx"))
  
  
  write_csv(table2, file.path(output_dir, "Table2.csv"))
  #save_as_docx(flext2, path="Table2.docx")
  
  
  xeno_n_human <- xeno_mf_tbl %>% count(xeno_class)
  
  
  
  #other_clades_to_label <- c("Mycobacteriaceae", "Paenibacillaceae", "Lachnospiraceae", "Enterobacteriaceae", "Clostridiaceae", "Bacillaceae_G", "Bacteroidaceae", "Ruminococcaceae", "Burkholderiaceae")
  other_clades_to_label <- c("f__Lachnospiraceae",
                             "f__Bacteroidaceae",
                             "f__Enterobacteriaceae", 
                             "f__Lactobacillaceae", 
                             "f__Ruminococcaceae", 
                             "f__Acutalibacteraceae", 
                             "f__Oscillospiraceae", 
                             "f__Burkholderiaceae", 
                             "f__Paenibacillaceae", 
                             "f__Mycobacteriaceae")
  tips_per_clade <- map(other_clades_to_label, ~ filter(species_md, family==paste0(.x))$Species_rep)
  names(tips_per_clade) <- other_clades_to_label
  tpc_sort <- tips_per_clade[order(map_int(tips_per_clade, length), decreasing = TRUE)]
  mrca_data <- enframe(map_int(tpc_sort, ~ape::getMRCA(gt_rooted, .x)), name="annotation", value="id") %>% mutate(row_num = row_number())
  
  
  ggtr <- ggtree(gt_rooted, layout='fan', open.angle=15)
  ggtr_cladelab <- ggtr + geom_cladelab(data=mrca_data, mapping=aes(label=row_num, node=id), offset=0, align=TRUE, barsize=2, offset.text=.175, barcolor="#888888", textcolor="#222222")
  
  #ggtr$data
  
  top_phyla <- species_md %>% group_by(phylum) %>% count() %>% arrange(-n) %>% .[1:10,] %>% .$phylum
  phyla_plot_data <- species_md %>% rename(Species_rep="tip") %>% select(tip, phylum) %>% mutate(phylum = map_chr(phylum, ~ ifelse(.x %in% top_phyla, .x, "other")))
  ggtr_phylum <- ggtr_cladelab + geom_fruit(data=phyla_plot_data, geom=geom_tile, mapping=aes(fill=phylum, y=tip, x=0), offset=0.2) + scale_fill_brewer(palette="Set3")
  colors_in_order <- c("green","oceanblue","puce","tan","salmon","orange2","darkgray","tomato","violetred","skyblue", "umber")
  xeno_cols <- chart_cols[colors_in_order]
  print(xeno_st_pd %>% as_tibble(rownames="xeno_class") %>% mutate(color = colors_in_order))
  
  ggtr_f <- ggtr_phylum; for (.r in 1:nrow(xeno_st_pd)) {
    r <- rownames(xeno_st_pd)[.r]
    .offset <- ifelse(.r > 1, -0.125, .25)
    ggtr_f <- ggtr_f + 
      new_scale_fill() +
      geom_fruit(show.legend=FALSE, data=xeno_na_tbl, geom=geom_tile, axis.params=list(color="black",axis="y"), mapping=aes(fill=.data[[r]], y=tip, x=-0.25), pwidth=0.25,width=0.25, offset=.offset) +
      scale_fill_manual(values=c(N="#FFFFFF00", Y=xeno_cols[[.r]]), na.value=NA)
  }
  
  pdf(width=10,height=10,file.path(output_dir, "xenobiotic_types.pdf"))
  print(ggtr_f)
  dev.off()
  
  
  family_xenos <- lapply(xeno_mf_all_hp,
                         function(x) {
                           full_best %>%
                             filter(hum_protein %in% x) %>%
                             collect() %>% 
                             select(hum_protein, Lineage) %>%
                             separate_wider_delim(Lineage, delim=";", names=c("domain","phylum","class","order","family","genus","species")) %>%
                             select(hum_protein, domain:species) %>%
                             distinct() %>%
                             group_by(family) %>%
                             summarize(n_sp = length(unique(species)), n_prot = length(hum_protein)) %>%
                             arrange(-n_sp)
                         })
  n_found_per_fam <- map2(family_xenos, names(family_xenos), ~ mutate(.x, xeno_class=.y)) %>% Reduce(bind_rows, .)
  top_3_per_fam <- n_found_per_fam %>% ungroup() %>% group_by(xeno_class) %>% slice_max(n_sp, n=3)
  top_3_per_fam_prot <- n_found_per_fam %>% ungroup() %>% group_by(xeno_class) %>% slice_max(n_prot, n=3)
  
  sorted_line_plot <- function(tbl) {
    tbl %>%
      arrange(hum_protein, sstart) %>%
      group_by(hum_protein) %>%
      mutate(r = row_number(), p=percent_rank(r)) %>%
      select(hum_protein, g, sstart, send, r, p, phylum, length_hum) %>%
      pivot_longer(values_to="location", names_to="which_end", sstart:send) %>%
      ungroup() %>%
      group_by(hum_protein) %>%
      mutate(xnorm=location/(length_hum)) %>%
      ggplot(aes(x=xnorm, y=p, group=r)) +
      geom_line() +
      geom_point(color="#000000", size=0.025, alpha=0.25) +
      facet_wrap(~hum_protein) +
      scale_color_brewer(palette="Set3")+
      theme_pander() +
      labs(x = "normalized position on human protein", y = "genomes") +
      scale_x_continuous(breaks = c(0, 1), minor_breaks = c(0,0.25,0.5,0.75,1)) + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank()) + 
      scale_linewidth(range=c(0.25,2))
  }
  
  sorted_mega_line_plot <- function(tbl) {
    tbl %>%
      arrange(hum_protein, Lineage, g, sstart) %>%
      group_by(hum_protein, g) %>%
      nest() %>%
      ungroup() %>%
      group_by(hum_protein) %>%
      mutate(r = row_number(), p=percent_rank(r)) %>%
      unnest() %>%
      ungroup() %>%
      group_by(hum_protein) %>%
      mutate(p_lwd = 1/(length(unique(g)))) %>%
      select(hum_protein, g, sstart, send, r, p, p_lwd, phylum, length_hum) %>%
      pivot_longer(values_to="location", names_to="which_end", sstart:send) %>%
      ungroup() %>%
      group_by(hum_protein) %>%
      mutate(xnorm=location/(length_hum)) %>%
      ggplot(aes(x=xnorm, y=p, group=r, color=phylum, lwd=p_lwd)) +
      geom_line() +
      geom_point(color="#000000", aes(size=p_lwd, alpha=p_lwd)) +
      facet_wrap(~hum_protein) +
      scale_color_brewer(palette="Set3")+
      theme_pander() +
      labs(x = "normalized position on human protein", y = "genomes") +
      scale_x_continuous(breaks = c(0, 1), minor_breaks = c(0,0.25,0.5,0.75,1)) + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank()) + 
      scale_linewidth(range=c(0.25,2)) +
      scale_size_continuous(range =  c(0.01,1.25)) + scale_alpha_continuous(range=c(0.1,1), transform="log2")
  }
  
  genes_mostly_split <- eggnog_part %>%
    separate_wider_delim(hum_protein,names= c("sp.a","ID","hum_protein"),delim = "|") %>%
    filter(hum_protein %in% mostly_split)
  
  genes_mostly_split_display <- genes_mostly_split %>%
    mutate(hum_protein = gsub("_HUMAN", "", hum_protein))
  
  # get correct ordering
  split_genes_ordering <- genes_mostly_split_display %>%
    #  select(hum_protein, Lineage) %>%
    select(hum_protein, g) %>%
    distinct() %>%
    group_by(hum_protein) %>%
    count() %>%
    arrange(-n)
  
  genes_mostly_split_display <- genes_mostly_split_display %>%
    mutate(hum_protein = factor(hum_protein, levels=split_genes_ordering$hum_protein)) %>%
    separate_wider_delim(Lineage, delim=";", names=c("domain","phylum","class","order","family","genus","species"), cols_remove = FALSE)
  
  
  split_top_phyla <- species_md %>% group_by(phylum) %>% count() %>% arrange(-n) %>% .[1:10,] %>% .$phylum
  split_display_data <- genes_mostly_split_display %>% mutate(phylum = map_chr(phylum, ~ ifelse(.x %in% top_phyla, .x, "other")))
  
  library(ggthemes)
  
  pdf(width=6.5, height=6.5, file.path(output_dir, "Figure2.pdf"))
  sorted_mega_line_plot(split_display_data)
  dev.off()
  
  sorted_line_plot(split_display_data)
  
  ### Now show distribution of these genes
  
  eggnog_upr_full %>% filter(`Entry Name` %in% mostly_split)
  if ((params[1] == 0.67) && (params[2] == 70)) { 
  
  
  genome_map <- genomes_md %>% select(Species_rep, Genome) %>% distinct()
  split_presabs_tbl_1 <- eggnog_part %>%
    filter(hum_protein %in% filter(eggnog_upr_full, nPart > nFull)$hum_protein) %>%
    group_by(hum_protein, Genome) %>%
    count() %>%
    filter(n>=2) %>%
    ungroup() %>%
    left_join(., genome_map, by="Genome") %>%
    group_by(hum_protein, Species_rep) %>%
    count() %>%
    mutate(present = ifelse(n>0, 1,0)) %>%
    select(hum_protein, Species_rep, present)
  
  split_presabs_tbl <- split_presabs_tbl_1 %>%
    pivot_wider(names_from=Species_rep, values_from=present, values_fill = 0) %>%
    right_join(select(eggnog_upr_full, hum_protein, `Entry Name`) , .) %>%
    select(-hum_protein)
  split_presabs_mtx <- split_presabs_tbl[-1] %>% data.matrix(); rownames(split_presabs_mtx) <- split_presabs_tbl[[1]]
  split_presabs_tree_tbl <- t(split_presabs_mtx) %>% as_tibble(rownames="tip")
  
  total_n_per_fam <- species_md %>% filter(Species_rep %in% gt_rooted$tip.label) %>% group_by(family) %>% count(name = "total_n")
  presence_by_family_b <- split_presabs_tree_tbl %>%
    group_by(tip) %>%
    summarize(n = sum(across(OPLA_HUMAN:THNS1_HUMAN))) %>%
    left_join(., species_md, by=c("tip"="Species_rep")) %>%
    ungroup() %>%
    group_by(family) %>%
    left_join(., total_n_per_fam) %>%
    summarize(total_det = sum(n),
              nspecies_det = sum(n>0),
              total_nspecies = unique(total_n),
              frac_det = nspecies_det / total_nspecies) %>%
    arrange(-nspecies_det)
  print(presence_by_family_b)
  
  other_clades_to_label_b <- c("Lachnospiraceae", "Ruminococcaceae", "Lactobacillaceae", "Oscillospiraceae", "Acutalibacteraceae", "Eggerthellaceae", "CAG-272", "CAG-727", "Bacillaceae_A", "Paenibacillaceae")
  tips_per_clade_b <- map(other_clades_to_label_b, ~ filter(species_md, family==paste0("f__", .x))$Species_rep)
  names(tips_per_clade_b) <- other_clades_to_label_b
  tpc_sort_b <- tips_per_clade_b[order(map_int(tips_per_clade_b, length), decreasing = TRUE)]
  mrca_data_b <- enframe(map_int(tpc_sort_b, ~ape::getMRCA(gt_rooted, .x)), name="annotation", value="id") %>% mutate(row_num = row_number())
  ggtr_cladelab_b <- ggtr + geom_cladelab(data=mrca_data_b, mapping=aes(label=row_num, node=id), offset=0, align=TRUE, barsize=2, offset.text=.175, barcolor="#888888", textcolor="#222222")
  ggtr_phylum_b <- ggtr_cladelab_b + geom_fruit(data=phyla_plot_data, geom=geom_tile, mapping=aes(fill=phylum, y=tip, x=0), offset=0.2) + scale_fill_brewer(palette="Set3")
  
  # note: a few are archaeal!!!
  n_bacteria <- split_presabs_mtx[, intersect(colnames(split_presabs_mtx), gt_rooted$tip.label)] %>% rowSums
  present_in_bac <- names(which(n_bacteria >= 1))
  split_presabs_pd <- picante::pd(split_presabs_mtx[present_in_bac, ], gt_rooted) %>% arrange(-PD)
  
  
  ggtr_split_f <- ggtr_phylum_b; for (.r in 1:15) { #nrow(split_presabs_pd)) { Do top 15 instead
    r <- rownames(split_presabs_pd)[.r]
    .offset <- ifelse(.r > 1, -.05, .2)
    ggtr_split_f <- ggtr_split_f + 
      new_scale_fill() +
      geom_fruit(show.legend=FALSE, data=split_presabs_tree_tbl[, c("tip", r)], geom=geom_tile, axis.params=list(color="black",axis="n"), mapping=aes(fill=.data[[r]], y=tip, x=-0.1), pwidth=0.1,width=0.1, offset=.offset) +
      scale_fill_steps(low="#FFFFFF", high="#000000", na.value=NA)
  }
  
  split_subset <- c("XDH_HUMAN", "DPYD_HUMAN")
  ggtr_split_f2 <- ggtr_phylum_b; for (.r in 1:2) { #nrow(split_presabs_pd)) { Do top 15 instead
    r <- split_subset[.r]
    .offset <- ifelse(.r > 1, -.05, .2)
    ggtr_split_f2 <- ggtr_split_f2 + 
      new_scale_fill() +
      geom_fruit(show.legend=FALSE, data=split_presabs_tree_tbl[, c("tip", r)], geom=geom_tile, axis.params=list(color="black",axis="n"), mapping=aes(fill=.data[[r]], y=tip, x=-0.1), pwidth=0.1,width=0.1, offset=.offset) +
      scale_fill_steps(low="#FFFFFF", high="#000000", na.value=NA)
  }
  
  pdf(width=10,height=10,file.path(output_dir, "xenobiotic_split_types.pdf"))
  print(ggtr_split_f2)
  dev.off()
  
  
  
  # Find DesE homologues
   desE <- full_best %>% filter(bac_protein=="GUT_GENOME228173_01934") %>% collect()
   print(desE)
   write_csv(desE, "SuppTable7.csv")
   desE_best <- desE %>% select(`Entry Name`, hum_protein) %>% distinct()
   print(desE_best)
  }
}

