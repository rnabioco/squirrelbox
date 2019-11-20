# prep samples
combined2 <- read_csv("squirrelBox/combined2.csv.gz")
dict <- combined2 %>% distinct(sample, state)

a <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Adrenal.tsv.gz") %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Adrenal") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
f <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Forebrain.tsv.gz") %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Forebrain") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
k <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Kidney.tsv.gz") %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Kidney") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
l <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Liver.tsv.gz") %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Liver") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
m <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Medulla.tsv.gz") %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Medulla") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
h <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Hypothalamus.tsv.gz") %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Hypothalamus") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)

combined <- do.call(rbind, list(a,f,k,l,m,h))
write_feather(combined, "combined2.feather")

a <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Adrenal.tsv.gz") %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
h <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Hypothalamus.tsv.gz") %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
f <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Forebrain.tsv.gz") %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
l <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Liver.tsv.gz") %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
m <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Medulla.tsv.gz") %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
k <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Kidney.tsv.gz") %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)

combined3 <- do.call(rbind, list(a,f,k,l,m,h)) %>% distinct()
write_feather(combined3, "combined3.feather")

bed <- read_tsv("squirrelBox/final_annot_20191112.bed12", col_names = F, col_types = "cnncncnnnnncccc") %>%
  select(1:6,13)
colnames(bed) <- c("chrom", "start", "end", "unique_gene_symbol", "score", "strand", "gene_id")
# gene_list <- bed %>% select(gene_id, chrom) %>% distinct(gene_id, .keep_all = TRUE)
# gene_list_orf <- gene_list %>% mutate(orf = pbmcapply::pbmcmapply(findorf, gene_id)) 
# saveRDS(gene_list_orf, "gene_id_orf11.rds")
gene_list_orf <- readRDS("gene_id_orf11.rds")
get_exons <- function(gene){
  #(bed %>% filter(gene_id == gene) %>% arrange(desc(exon_number)))$exon_number[1] + 1
  bed %>% filter(gene_id == gene) %>% nrow()
}

get_length <- function(gene){
  sum((bed %>% filter(gene_id == gene) %>% bed_merge() %>% mutate(length = end - start))$length)
}
gene_list_orf2 <- gene_list_orf %>%
  mutate(orf = ifelse(sapply(gene_list_orf$orf, length) == 1, orf, "")) %>%
  mutate(orf = unlist(orf)) %>% mutate(len = str_length(orf)) %>%
  mutate(exons = mapply(get_exons, gene_id), transcript = mapply(get_length, gene_id))

a2 <- gene_list_orf2
dataref <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Hypothalamus.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("hy_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a2, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Forebrain.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("fore_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESeq2_salmon_rlog_20191113/Medulla.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("med_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
d4 <- bed %>% select(gene_id, unique_gene_symbol) %>% distinct()
a4 <- a3 %>% left_join(d4, by = "gene_id")
write_csv(a4, "padj_orf11.csv")