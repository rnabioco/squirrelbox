require(tidyverse)
require(feather)

# prep samples
pheno_path <- c("/Users/rf/sandy/brain_phenodata.csv", "/Users/rf/sandy/medulla_phenodata.csv", "/Users/rf/sandy/hypothalamus_phenodata.csv")
pheno <- map(pheno_path, function(x) {
  read.csv(x) %>% select(ids, population)
})
dict <- do.call(rbind, pheno) %>%
  mutate(ids = str_remove(ids, "[A-Za-z]")) %>%
  distinct()

# expression
loc <- "/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130"
f <- read_tsv(paste0(loc, "/Forebrain.tsv.gz")) %>%
  select(-contains("_"), -source, -cluster, gene_id) %>%
  mutate(region = "Forebrain") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>%
  left_join(dict %>% mutate(ids = str_c("B", ids, sep = "")), by = c("sample" = "ids"))
m <- read_tsv(paste0(loc, "/Medulla.tsv.gz")) %>%
  select(-contains("_"), -source, -cluster, gene_id) %>%
  mutate(region = "Medulla") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>%
  left_join(dict %>% mutate(ids = str_c("M", ids, sep = "")), by = c("sample" = "ids"))
h <- read_tsv(paste0(loc, "/Hypothalamus.tsv.gz")) %>%
  select(-contains("_"), -source, -cluster, gene_id) %>%
  mutate(region = "Hypothalamus") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>%
  left_join(dict %>% mutate(ids = str_c("H", ids, sep = "")), by = c("sample" = "ids"))

combined <- do.call(rbind, list(f, m, h)) %>% rename(state = "population")
write_feather(combined, "combined2.feather") # <- combined2.feather

# gene_info
h <- read_tsv(paste0(loc, "/Hypothalamus.tsv.gz")) %>%
  select(gene_id, unique_gene_symbol, gene_symbol, clean_gene_symbol, original_gene_name, source)
f <- read_tsv(paste0(loc, "/Forebrain.tsv.gz")) %>%
  select(gene_id, unique_gene_symbol, gene_symbol, clean_gene_symbol, original_gene_name, source)
m <- read_tsv(paste0(loc, "/Medulla.tsv.gz")) %>%
  select(gene_id, unique_gene_symbol, gene_symbol, clean_gene_symbol, original_gene_name, source)

combined3 <- do.call(rbind, list(f, h, m)) %>% distinct()
write_feather(combined3, "combined3.feather") # <- combined3.feather

# cluster_info
mod_f <- read_tsv(paste0(loc, "/Forebrain.tsv.gz")) %>%
  select(gene = unique_gene_symbol, cluster) %>%
  mutate(cluster = ifelse(is.na(cluster), "insig", cluster))
mod_h <- read_tsv(paste0(loc, "/Hypothalamus.tsv.gz")) %>%
  select(gene = unique_gene_symbol, cluster) %>%
  mutate(cluster = ifelse(is.na(cluster), "insig", cluster))
mod_m <- read_tsv(paste0(loc, "/Medulla.tsv.gz")) %>%
  select(gene = unique_gene_symbol, cluster) %>%
  mutate(cluster = ifelse(is.na(cluster), "insig", cluster))
mod <- mod_f %>%
  full_join(mod_h, by = "gene", suffix = c("_f", "_h")) %>%
  full_join(mod_m, by = "gene", suffix = c("", "_m")) %>%
  rename(cluster_m = cluster) %>%
  mutate(cluster_f = ifelse(is.na(cluster_f), "filtered", cluster_f)) %>%
  mutate(cluster_h = ifelse(is.na(cluster_h), "filtered", cluster_h)) %>%
  mutate(cluster_m = ifelse(is.na(cluster_m), "filtered", cluster_m))
write_feather(mod, "clusters.feather")

# orf
bed <- read_tsv("/Users/rf/Downloads/final_tx_annotations_20200201.tsv.gz",
  col_names = T,
  col_types = list(blockSizes = "c")
) %>%
  select(1:6, 10:14, 19)
colnames(bed) <- c(
  "chrom", "start", "end", "transcript_id", "score",
  "strand", "exons", "e_l", "site", "gene_id", "unique_gene_symbol", "majiq"
)

get_length <- function(e_l) {
  vec1 <- as.numeric(head(unlist(str_split(e_l, ",")), -1))
  sum(vec1)
}

get_rc <- function(strand, seq) {
  if (strand == "-") {
    seq <- chartr("ATGC", "TACG", seq)
    seq <- stringi::stri_reverse(seq)
  }
  seq
}

get_translation <- function(seq, fromnum) {
  as.character(translate(DNAString(str_sub(seq, fromnum)), if.fuzzy.codon = "solve"))
}

findorf <- function(chrom, start, end, strand, e_l, site) {
  s <- start
  temp_tbl <- data.frame(
    "start" = as.numeric(head(unlist(str_split(site, ",")), -1)),
    "length" = as.numeric(head(unlist(str_split(e_l, ",")), -1))
  ) %>%
    mutate(start = start + s) %>%
    mutate(end = start + length)
  tbl_gr <- GRanges(
    seqnames = Rle(chrom),
    ranges = IRanges(temp_tbl$start + 1, end = temp_tbl$end)
  )
  seq <- paste(getSeq(BSgenome.SQ1, tbl_gr, as.character = TRUE), collapse = "")
  if (strand == "-") {
    seq <- get_rc("-", seq)
  }
  aaseq1 <- get_translation(seq, 1)
  aaseq2 <- get_translation(seq, 2)
  aaseq3 <- get_translation(seq, 3)
  aaall <- c(unlist(str_split(aaseq1, "\\*")), unlist(str_split(aaseq2, "\\*")), unlist(str_split(aaseq3, "\\*")))
  aaall <- na.omit(str_extract(aaall, "M.*"))
  aalengths <- str_length(aaall)
  aalengthsmax <- which.max(aalengths)
  aaall[aalengthsmax]
}
library(pbmcapply)
library("BSgenome")
library("BSgenome.SQ1")
a <- bed %>% mutate(orf = pbmcapply::pbmcmapply(findorf, chrom, start, end, strand, e_l, site))
a2 <- a %>%
  mutate(orf = ifelse(sapply(a$orf, length) == 1, orf, "")) %>%
  mutate(orf = as.character(unlist(orf))) %>%
  mutate(len = str_length(orf)) %>%
  mutate(transcript = mapply(get_length, e_l))
saveRDS(a2, "orf_frombed.rds")

a2 <- readRDS("orf_frombed.rds")
a2 <- a2 %>%
  select(chrom, start, end, site, orf, len, transcript) %>%
  left_join(bed)

dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Hypothalamus.tsv.gz")
d <- dataref %>%
  select(gene_id, contains("_padj")) %>%
  select_all(.funs = funs(paste0("hy_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a2, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Forebrain.tsv.gz")
d <- dataref %>%
  select(gene_id, contains("_padj")) %>%
  select_all(.funs = funs(paste0("fb_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Medulla.tsv.gz")
d <- dataref %>%
  select(gene_id, contains("_padj")) %>%
  select_all(.funs = funs(paste0("med_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Adrenal.tsv.gz")
d <- dataref %>%
  select(gene_id, contains("_padj")) %>%
  select_all(.funs = funs(paste0("adr_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Liver.tsv.gz")
d <- dataref %>%
  select(gene_id, contains("_padj")) %>%
  select_all(.funs = funs(paste0("liv_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Kidney.tsv.gz")
d <- dataref %>%
  select(gene_id, contains("_padj")) %>%
  select_all(.funs = funs(paste0("kid_", .)))
colnames(d)[1] <- "gene_id"
padj_orf <- left_join(a3, d, by = c("gene_id" = "gene_id"))
feather::write_feather(padj_orf, "squirrelBox/data/padj_orf.feather")

# domain annotation
a3 <- read_csv("/Users/rf/newselect201905/sq_gen/sq_gen/squirrelBox/padj_orf.csv")
a4 <- a3 %>%
  filter(str_detect(unique_gene_symbol, "_like|_contain|^G[0-9]+")) %>%
  filter(len > 0) %>%
  group_by(gene_id) %>%
  arrange(desc(len)) %>%
  dplyr::slice(1)
writefasta <- function(full_list, name) {
  for (n in c(1:nrow(full_list))) {
    gene <- as.character(full_list[n, 10])
    print(gene)
    seq <- as.character(full_list[n, 11])
    print(seq)
    write(paste0(">", gene), file = name, append = TRUE)
    write(seq, file = name, append = TRUE)
  }
}
writefasta(a4, "novels2020.fa")
# hmmsearch --tblout novelse2_2020 --cpu 6 -E 1e-2 /Users/rf/Downloads/Pfam-A.hmm /Users/rf/newselect201905/sq_gen/sq_gen/novels2020.fa
domains <- read_table2("/Users/rf/novelse2_2020", comment = "#", col_names = FALSE) %>%
  mutate(gene_id = X1, domain = X3) %>%
  select(gene_id, domain)
write_csv(domains, "novel_domains.csv")

# magik
mag <- read_tsv("/Users/rf/Downloads/clustered_retained_inton_lsvs/MAJIQ_dpsi_summary_sig_squirrelBox.tsv.gz")

# gene lists from transition pairs
loc <- "/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130"
state_order <- c(
  "SA",
  "IBA",
  "Ent",
  "LT",
  "Ar",
  "SpD"
)
h <- read_tsv(paste0(loc, "/Hypothalamus.tsv.gz")) %>%
  select(unique_gene_symbol, contains("ashr"), contains("wald_padj"))
h2 <- h %>%
  pivot_longer(-unique_gene_symbol,
    names_to = c("state", "type"),
    names_pattern = "([a-zA-Z]+_vs_[a-zA-Z]+)_(.*)",
    values_to = "val"
  ) %>%
  pivot_wider(names_from = "type", values_from = "val") %>%
  rename(gene = "unique_gene_symbol") %>%
  mutate_at(vars(contains("_")), as.numeric)
h3 <- split(h2, h2$state)

f <- read_tsv(paste0(loc, "/Forebrain.tsv.gz")) %>%
  select(unique_gene_symbol, contains("ashr"), contains("wald_padj"))
f2 <- f %>%
  pivot_longer(-unique_gene_symbol,
    names_to = c("state", "type"),
    names_pattern = "([a-zA-Z]+_vs_[a-zA-Z]+)_(.*)",
    values_to = "val"
  ) %>%
  pivot_wider(names_from = "type", values_from = "val") %>%
  rename(gene = "unique_gene_symbol") %>%
  mutate_at(vars(contains("_")), as.numeric)
f3 <- split(f2, f2$state)

m <- read_tsv(paste0(loc, "/Medulla.tsv.gz")) %>%
  select(unique_gene_symbol, contains("ashr"), contains("wald_padj"))
m2 <- m %>%
  pivot_longer(-unique_gene_symbol,
    names_to = c("state", "type"),
    names_pattern = "([a-zA-Z]+_vs_[a-zA-Z]+)_(.*)",
    values_to = "val"
  ) %>%
  pivot_wider(names_from = "type", values_from = "val") %>%
  rename(gene = "unique_gene_symbol") %>%
  mutate_at(vars(contains("_")), as.numeric)
m3 <- split(m2, m2$state)

write_sig <- function(df, p = 0.001, region) {
  state <- df$state[1]
  up <- df %>%
    filter(wald_padj <= p) %>%
    filter(log2FC_ashr > 0)
  down <- df %>%
    filter(wald_padj <= p) %>%
    filter(log2FC_ashr < 0)
  write_csv(
    up,
    paste0(region, "_", state, "_up.csv")
  )
  write_csv(
    down,
    paste0(region, "_", state, "_down.csv")
  )
}

sapply(h3, function(x) write_sig(x, region = "hy"))
sapply(f3, function(x) write_sig(x, region = "fb"))
sapply(m3, function(x) write_sig(x, region = "med"))
