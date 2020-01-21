# prep samples
require(tidyverse)
require(feather)
combined2 <- read_feather("squirrelBox/combined2.feather")
dict <- combined2 %>% distinct(sample, state)
dict2 <- read_tsv("/Users/rf/Downloads/SequenceFileInfoKidney.txt") %>% 
  filter(SampleID != "Sample", !is.na(SampleID)) %>% 
  select(sample = SampleID, state = State)
dict3 <- read_tsv("/Users/rf/Downloads/AllKidneySamples.txt") %>% 
  select(sample = Sample, state = group)
dict <- rbind(dict, dict2) %>% rbind(dict3) %>% 
  distinct()
# pheno_path <- c("/Users/rf/sandy/brain_phenodata.csv","/Users/rf/sandy/medulla_phenodata.csv","/Users/rf/sandy/hypothalamus_phenodata.csv")
# pheno <- map(pheno_path, function(x) {read.csv(x) %>% select(ids, population)})
# dict <- do.call(rbind, pheno) %>% mutate(ids = str_remove(ids, "[A-Za-z]")) %>% distinct()
loc <- "/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200117"

a <- read_tsv(paste0(loc, "/Adrenal.tsv.gz")) %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Adrenal") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
f <- read_tsv(paste0(loc, "/Forebrain.tsv.gz")) %>%
  select(-contains("_"), -source, -cluster, gene_id) %>%
  mutate(region = "Forebrain") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
k <- read_tsv(paste0(loc, "/Kidney.tsv.gz")) %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Kidney") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
l <- read_tsv(paste0(loc, "/Liver.tsv.gz")) %>%
  select(-contains("_"), -source, gene_id) %>%
  mutate(region = "Liver") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
m <- read_tsv(paste0(loc, "/Medulla.tsv.gz")) %>%
  select(-contains("_"), -source, -cluster, gene_id) %>%
  mutate(region = "Medulla") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)
h <- read_tsv(paste0(loc, "/Hypothalamus.tsv.gz")) %>%
  select(-contains("_"), -source, -cluster, gene_id) %>%
  mutate(region = "Hypothalamus") %>%
  pivot_longer(-c(gene_id, region), names_to = "sample", values_to = "log2_counts") %>% 
  left_join(dict)

combined <- do.call(rbind, list(a,f,k,l,m,h))
write_feather(combined, "combined2.feather") # <- combined2.feather

a <- read_tsv(paste0(loc, "/Adrenal.tsv.gz")) %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
h <- read_tsv(paste0(loc, "/Hypothalamus.tsv.gz")) %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
f <- read_tsv(paste0(loc, "/Forebrain.tsv.gz")) %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
l <- read_tsv(paste0(loc, "/Liver.tsv.gz")) %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
m <- read_tsv(paste0(loc, "/Medulla.tsv.gz")) %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)
k <- read_tsv(paste0(loc, "/Kidney.tsv.gz")) %>%
  select(gene_id,unique_gene_symbol,gene_symbol,clean_gene_symbol,original_gene_name,source)

combined3 <- do.call(rbind, list(a,f,k,l,m,h)) %>% distinct()
write_feather(combined3, "combined3.feather") # <- combined3.feather

mod_f <- read_tsv(paste0(loc, "/Forebrain.tsv.gz")) %>%
  select(gene = unique_gene_symbol, cluster) %>%
  mutate(cluster = ifelse(is.na(cluster), "insig", cluster))
mod_h <- read_tsv(paste0(loc, "/Hypothalamus.tsv.gz")) %>%
  select(gene = unique_gene_symbol, cluster) %>%
  mutate(cluster = ifelse(is.na(cluster), "insig", cluster))
mod_m <- read_tsv(paste0(loc, "/Medulla.tsv.gz")) %>%
  select(gene = unique_gene_symbol, cluster) %>%
  mutate(cluster = ifelse(is.na(cluster), "insig", cluster))
mod <- mod_f %>% full_join(mod_h, by = "gene", suffix = c("_f", "_h")) %>%
  full_join(mod_m, by = "gene", suffix = c("", "_m")) %>% 
  rename(cluster_m = cluster) %>% 
  mutate(cluster_f = ifelse(is.na(cluster_f), "filtered", cluster_f)) %>% 
  mutate(cluster_h = ifelse(is.na(cluster_h), "filtered", cluster_h)) %>% 
  mutate(cluster_m = ifelse(is.na(cluster_m), "filtered", cluster_m))
write_feather(mod, "clusters.feather")
# bed <- read_tsv("squirrelBox/final_annot_20191112.bed12", col_names = F, col_types = "cnncncnnnnncccc") %>%
#   select(1:6,13)
# colnames(bed) <- c("chrom", "start", "end", "unique_gene_symbol", "score", "strand", "gene_id")
bed <- read_tsv("squirrelBox/final_annot_20191112.bed12", col_names = F, col_types = list(X11 = "c", X12 = "c")) %>%
  select(1:6,10:13)
colnames(bed) <- c("chrom", "start", "end", "unique_gene_symbol", "score", "strand", "exons", "e_l", "site", "gene_id")

get_length <- function(e_l){
  vec1 <- as.numeric(head(unlist(str_split(e_l, ",")),-1))
  sum(vec1)
}

get_rc <- function(strand, seq){
  if(strand == "-"){
    seq <- chartr("ATGC","TACG", seq)
    seq <- stringi::stri_reverse(seq)
  }
  seq
}

get_translation <- function(seq, fromnum){
  as.character(translate(DNAString(str_sub(seq,fromnum)), if.fuzzy.codon = "solve"))
}

findorf <- function(chrom, start, end, strand, e_l, site){
  s <- start
  temp_tbl <- data.frame("start" = as.numeric(head(unlist(str_split(site, ",")), -1)),
                         "length" = as.numeric(head(unlist(str_split(e_l, ",")), -1))) %>%
    mutate(start = start + s) %>%
    mutate(end = start + length)
  tbl_gr <- GRanges(seqnames = Rle(chrom),
                    ranges = IRanges(temp_tbl$start + 1, end = temp_tbl$end))
  seq <- paste(getSeq(BSgenome.SQ1, tbl_gr, as.character = TRUE), collapse = '')
  if(strand == "-"){
    seq <- get_rc("-", seq)
  }
  aaseq1 <- get_translation(seq,1)
  aaseq2 <- get_translation(seq,2)
  aaseq3 <- get_translation(seq,3)
  aaall <- c(unlist(str_split(aaseq1, "\\*")), unlist(str_split(aaseq2, "\\*")), unlist(str_split(aaseq3, "\\*")))
  aaall <- na.omit(str_extract(aaall, "M.*"))
  aalengths <- str_length(aaall)
  aalengthsmax <- which.max(aalengths)
  aaall[aalengthsmax]
}

library("BSgenome")
library("BSgenome.SQ1")
a <- bed %>% mutate(orf = pbmcapply::pbmcmapply(findorf, chrom, start, end, strand, e_l, site))
a2 <- a %>% mutate(orf = ifelse(sapply(a$orf, length) == 1, orf, "")) %>% mutate(orf = as.character(unlist(orf))) %>% mutate(len = str_length(orf)) %>% mutate(transcript = mapply(get_length, e_l))
saveRDS(a2, "orf_frombed.rds")

# gene_list <- bed %>% select(gene_id, chrom) %>% distinct(gene_id, .keep_all = TRUE)
# gene_list_orf <- gene_list %>% mutate(orf = pbmcapply::pbmcmapply(findorf, gene_id)) 
# saveRDS(gene_list_orf, "gene_id_orf11.rds")
# gene_list_orf <- readRDS("gene_id_orf11.rds")
# get_exons <- function(gene){
#   #(bed %>% filter(gene_id == gene) %>% arrange(desc(exon_number)))$exon_number[1] + 1
#   bed %>% filter(gene_id == gene) %>% nrow()
# }
# 
# get_length <- function(gene){
#   sum((bed %>% filter(gene_id == gene) %>% bed_merge() %>% mutate(length = end - start))$length)
# }
# gene_list_orf2 <- gene_list_orf %>%
#   mutate(orf = ifelse(sapply(gene_list_orf$orf, length) == 1, orf, "")) %>%
#   mutate(orf = unlist(orf)) %>% mutate(len = str_length(orf)) %>%
#   mutate(exons = mapply(get_exons, gene_id), transcript = mapply(get_length, gene_id))
# 
# a2 <- gene_list_orf2
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
#d4 <- bed %>% select(gene_id, unique_gene_symbol) %>% distinct()
#a4 <- a3 %>% left_join(d4, by = "gene_id")
write_csv(a3, "padj_orf.csv") # <- padj_orf.csv

find_elbow <- function(x, y, after) {
  d1 <- diff(y) / diff(x) # first derivative
  d2 <- diff(d1) / diff(x[-1]) # second derivative
  maxd <- max(d2[after:length(d2)])
  mind <- min(d2[after:length(d2)])
  print(maxd)
  print(mind)
  indices <- which(d2 == maxd)  
  return(indices)
}

r1 <- read_csv("/Users/rf/newselect201911/hy/hy_mean_5000.csv", col_names = F)
r2 <- read_csv("/Users/rf/newselect201911/med/med_mean_5000.csv", col_names = F)
r3 <- read_csv("/Users/rf/newselect201911/fore/fore_mean_5000.csv", col_names = F)
r <- full_join(r1, r2, by = "X1") %>% full_join(r3, by = "X1")
r[is.na(r)] <- 0
colnames(r) <- c("unique_gene_symbol", "hy_imp", "med_imp", "fore_imp")
r <- r %>% mutate(hy_rank = rank(-hy_imp), med_rank = rank(-med_imp), fore_rank = rank(-fore_imp))
write_csv(r, "rf_imp.csv")

qsave(combined2, "combined2.qs")
combined2 <- qread("combined2.qs")