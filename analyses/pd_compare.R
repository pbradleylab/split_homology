# pd comparison
param_list <- list(c(0.5,60), c(0.67,70), c(0.75,80))
pd_results <- map(param_list, ~ {
  b <- .x[1]
  h <- .x[2]
  pd_tbl <- read_tsv(paste0("output/B", b, "/H", h, "/pd.tsv")) %>%
    mutate(param = paste0("B",b*100,"H",h))
}) %>%
  bind_rows
pd_wide <- pd_results %>% pivot_wider(names_from=param, values_from=c("PD","SR"))
write_csv(pd_wide, "SuppTable8.csv")

# go comparison
go_full_results <- map(param_list, ~ {
  b <- .x[1]
  h <- .x[2]
  pd_tbl <- read_tsv(paste0("output/B", b, "/H", h, "/tbl-full-vs-all-enr.tsv")) %>%
    mutate(param = paste0("B",b*100,"H",h))
}) %>%
  bind_rows
go_full_wide <- go_full_results %>% pivot_wider(names_from=param, values_from=corr_pval) %>% arrange(B67H70)
write_csv(go_full_wide, "SuppTable4.csv")

go_split_results <- map(param_list, ~ {
  b <- .x[1]
  h <- .x[2]
  pd_tbl <- read_tsv(paste0("output/B", b, "/H", h, "/tbl-split-vs-all-enr.tsv")) %>%
    mutate(param = paste0("B",b*100,"H",h))
}) %>%
  bind_rows
go_split_wide <- go_split_results %>% pivot_wider(names_from=param, values_from=corr_pval) %>% arrange(B67H70)
write_csv(go_split_wide, "SuppTable6.csv")

# convert formats for publication
st1 <- read_tsv("output/B0.67/H70/SuppTable1.tsv")
write_csv(st1, "SuppTable1.csv")
st3 <- read_tsv("output/B0.67/H70/SuppTable3.tsv")
write_csv(st3, "SuppTable3.csv")
st5 <- read_tsv("output/B0.67/H70/SuppTable5.tsv")
write_csv(st3, "SuppTable5.csv")