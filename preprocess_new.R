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
loc <- "/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130"

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
bed <- read_tsv("/Users/rf/Downloads/final_tx_annotations_20200201.tsv.gz", 
                col_names = T,
                col_types = list(blockSizes= "c")) %>%
  select(1:6,10:14,19)
colnames(bed) <- c("chrom", "start", "end", "transcript_id", "score", 
                   "strand", "exons", "e_l", "site", "gene_id", "unique_gene_symbol", "majiq")

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
library(pbmcapply)
library("BSgenome")
library("BSgenome.SQ1")
a <- bed %>% mutate(orf = pbmcapply::pbmcmapply(findorf, chrom, start, end, strand, e_l, site))
a2 <- a %>% mutate(orf = ifelse(sapply(a$orf, length) == 1, orf, "")) %>% mutate(orf = as.character(unlist(orf))) %>% mutate(len = str_length(orf)) %>% mutate(transcript = mapply(get_length, e_l))
saveRDS(a2, "orf_frombed.rds")

a2 <- readRDS("orf_frombed.rds")
a2 <- a2 %>% select(chrom, start, end, site, orf, len, transcript) %>% left_join(bed)

dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Hypothalamus.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("hy_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a2, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Forebrain.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("fore_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Medulla.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("med_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Adrenal.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("adr_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Liver.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("liv_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
dataref <- read_tsv("/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130/Kidney.tsv.gz")
d <- dataref %>% select(gene_id, contains("_padj")) %>% select_all(.funs = funs(paste0("kid_", .)))
colnames(d)[1] <- "gene_id"
a3 <- left_join(a3, d, by = c("gene_id" = "gene_id"))
write_csv(a3, "padj_orf.csv")

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


a3 <- read_csv("/Users/rf/newselect201905/sq_gen/sq_gen/squirrelBox/padj_orf.csv")
a4 <- a3 %>% filter(str_detect(unique_gene_symbol, "_like|_contain|^G[0-9]+")) %>%
  filter(len > 0) %>%
  group_by(gene_id) %>%
  arrange(desc(len)) %>% dplyr::slice(1)
writefasta <- function(full_list, name){
  for(n in c(1:nrow(full_list))){
    gene <- as.character(full_list[n,10])
    print(gene)
    seq <- as.character(full_list[n,11])
    print(seq)
    write(paste0(">",gene), file = name, append=TRUE)
    write(seq, file = name, append=TRUE)
  }
}

writefasta(a4, "novels2020.fa")

# hmmsearch --tblout novelse2_2020 --cpu 6 -E 1e-2 /Users/rf/Downloads/Pfam-A.hmm /Users/rf/newselect201905/sq_gen/sq_gen/novels2020.fa
domains <- read_table2("/Users/rf/novelse2_2020", comment = "#", col_names = FALSE) %>% mutate(gene_id = X1, domain = X3) %>% select(gene_id, domain)
write_csv(domains, "novel_domains.csv")

# venn diagram
library(eulerr)
loc <- "/Users/rf/sandy/newtab/DESesq2_salmon_rlog_20200130"
a <- read_tsv(paste0(loc, "/Adrenal.tsv.gz"))
f <- read_tsv(paste0(loc, "/Forebrain.tsv.gz"))
k <- read_tsv(paste0(loc, "/Kidney.tsv.gz"))
l <- read_tsv(paste0(loc, "/Liver.tsv.gz"))
m <- read_tsv(paste0(loc, "/Medulla.tsv.gz"))
h <- read_tsv(paste0(loc, "/Hypothalamus.tsv.gz"))

data1 <- rbind(data.frame(gene = f$gene_id, test = "fore"),
              data.frame(gene = h$gene_id, test = "hy"),
              data.frame(gene = m$gene_id, test = "med"))
data1 <- data1 %>% group_by(gene) %>% 
  summarize(test = bazar::concat(test, sep = "&"), n = n())

counts1 <- table(data1$test)
ps1 <- as.vector(counts1)
names(ps1) <- names(counts1)
vd1 <- euler(ps1)
plot(vd1, labels = list(cex = .5), quantities = list(cex = .5), adjust_labels = FALSE)
pdf(file = "./3brainvd_small.pdf", width = 3, height = 3)
plot(vd1, labels = list(cex = .5), quantities = list(cex = .5), adjust_labels = FALSE)
dev.off()
a <- read_tsv(paste0(loc, "/Adrenal_significant_001.tsv.gz"))
f <- read_tsv(paste0(loc, "/Forebrain_significant_001.tsv.gz"))
k <- read_tsv(paste0(loc, "/Kidney_significant_001.tsv.gz"))
l <- read_tsv(paste0(loc, "/Liver_significant_001.tsv.gz"))
m <- read_tsv(paste0(loc, "/Medulla_significant_001.tsv.gz"))
h <- read_tsv(paste0(loc, "/Hypothalamus_significant_001.tsv.gz"))
# data2 <- rbind(data.frame(gene = f %>% filter(LRT_padj <= 0.001) %>% pull(gene_id), test = "fore"),
#                data.frame(gene = h %>% filter(LRT_padj <= 0.001) %>% pull(gene_id), test = "hy"),
#                data.frame(gene = m %>% filter(LRT_padj <= 0.001) %>% pull(gene_id), test = "med"))
data2 <- rbind(data.frame(gene = f$gene_id, test = "fore"),
               data.frame(gene = h$gene_id, test = "hy"),
               data.frame(gene = m$gene_id, test = "med"))
data2 <- data2 %>% group_by(gene) %>% 
  summarize(test = bazar::concat(test, sep = "&"), n = n())

counts2 <- table(data2$test)
ps2 <- as.vector(counts2)
names(ps2) <- names(counts2)
vd2 <- euler(ps2)
plot(vd2, labels = list(cex = .5), quantities = list(cex = .5), adjust_labels = FALSE)

pdf(file = "./sig3DEbrainvd_small.pdf", width = 3, height = 3)
plot(vd2, labels = list(cex = .5), quantities = list(cex = .5), adjust_labels = FALSE)
dev.off()

# magik
mag <- read_tsv('/Users/rf/Downloads/clustered_retained_inton_lsvs/MAJIQ_dpsi_summary_sig_squirrelBox.tsv.gz')

#
a <-dataref %>% filter(IBA_vs_LT_wald_padj <= 0.001, IBA_vs_LT_log2FC_ashr > 0)
b <-dataref %>% filter(IBA_vs_LT_wald_padj <= 0.001, IBA_vs_LT_log2FC_ashr < 0)
write_csv(a %>% select(unique_gene_symbol), "med_IBAoverLT.csv")
write_csv(b %>% select(unique_gene_symbol), "med_IBAunderLT.csv")

intersect(read_csv("med_IBAoverLT.csv") %>% pull(1), 
  read_csv("hy_IBAoverLT.csv") %>% pull(1)) %>% intersect(read_csv("fore_IBAoverLT.csv") %>% pull(1)) -> IBAup

intersect(read_csv("med_IBAunderLT.csv") %>% pull(1), 
          read_csv("hy_IBAunderLT.csv") %>% pull(1)) %>% intersect(read_csv("fore_IBAunderLT.csv") %>% pull(1)) -> IBAdown

c(read_csv("med_IBAoverLT.csv") %>% pull(1), 
  read_csv("hy_IBAoverLT.csv") %>% pull(1),
  read_csv("fore_IBAoverLT.csv") %>% pull(1)) %>% unique() -> IBAup
c(read_csv("med_IBAunderLT.csv") %>% pull(1), 
  read_csv("hy_IBAunderLT.csv") %>% pull(1),
  read_csv("fore_IBAunderLT.csv") %>% pull(1)) %>% unique() -> IBAdown
write_lines(IBAdown, "IBAdown.txt")
write_lines(IBAup, "IBAup.txt")

# finding transition pairs
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
               values_to = "val") %>% 
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
               values_to = "val") %>% 
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
               values_to = "val") %>% 
  pivot_wider(names_from = "type", values_from = "val") %>% 
  rename(gene = "unique_gene_symbol") %>% 
  mutate_at(vars(contains("_")), as.numeric)
m3 <- split(m2, m2$state)

write_sig <- function(df, p = 0.001, region) {
  state <- df$state[1]
  up <- df %>% filter(wald_padj <= p) %>% 
    filter(log2FC_ashr > 0)
  down <- df %>% filter(wald_padj <= p) %>% 
    filter(log2FC_ashr < 0)
  write_csv(up, 
            paste0(region, "_", state, "_up.csv"))
  write_csv(down, 
            paste0(region, "_", state, "_down.csv"))
}

sapply(h3, function(x) write_sig(x, region = "hy"))
sapply(f3, function(x) write_sig(x, region = "fore"))
sapply(m3, function(x) write_sig(x, region = "med"))

