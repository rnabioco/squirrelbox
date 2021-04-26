library(shiny)
library(dplyr)
library(tibble)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(feather)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap) # bioconductor
library(DT)
library(igraph)
library(valr)
library(transite) # bioconductor
library(plotly)
library(crosstalk)
library(shinyjs)
library(shinythemes)
library(shinycustomloader)
library(shinyjqui)
library(shinyWidgets)
library(rintrojs)
library(shinyBS)
library(bsplus)
library(BioCircos)

### source R
source("config.R")
source(paste0(rpath, "/ggvenn.R"))

if (file.exists(paste0(annotpath, "/genes.csv"))) {
  genes_df <- read_csv(paste0(annotpath, "/genes.csv"), col_names = FALSE, comment = "#")
  if (ncol(genes_df) > 1) {
    genes_df[[2]] <- ifelse(is.na(genes_df[[2]]), genes_df[[1]], genes_df[[2]])
    genes_tips <- genes_df %>%
      pull(2) %>%
      str_c(collapse = "\', \'")
    genes_tips <- paste0("\'", genes_tips, "\'")
  } else {
    genes_tips <- genes_df[[1]]
  }
}

if (file.exists(paste0(annotpath, "/alt.csv"))) {
  maj_df <- read_csv(paste0(annotpath, "/alt.csv"), col_names = FALSE)
  if (ncol(maj_df) > 1) {
    maj_df[[2]] <- ifelse(is.na(maj_df[[2]]), maj_df[[1]], maj_df[[2]])
    maj_tips <- maj_df %>%
      pull(2) %>%
      str_c(collapse = "\', \'")
    maj_tips <- paste0("\'", maj_tips, "\'")
  } else {
    maj_tips <- maj_df[[1]]
  }
}

# read database
if (file.exists(paste0(datapath, "/combined2.feather"))) {
  combined2 <- read_feather(paste0(datapath, "/combined2.feather"))
  combined2_akl <- read_feather(paste0(datapath, "/combined2_akl.feather"))
  combined2 <- bind_rows(combined2, combined2_akl)
  combined3 <- read_feather(paste0(datapath, "/full_combined3.feather"))
} else if (file.exists(paste0(datapath, "/combined2.csv"))) {
  combined2 <- fread(paste0(datapath, "/combined2.csv"), nThread = ncore)
  combined3 <- fread(paste0(datapath, "/combined3.csv"), nThread = ncore)
}

# use corrected
combined2_cor <- read_feather(paste0(datapath, "/full_combined2_cor.feather"))
combined2 <- bind_rows(combined2_cor, combined2 %>% filter(!(region %in% (combined2_cor$region %>% unique()))))

# fix EAr
combined2 <- combined2 %>% mutate(state = ifelse(state == "EAr", "Ar", as.character(state))) %>% filter(state != "ST")

# groseq 
combined_gro <- read_feather(paste0(datapath, "/combined_groseq.feather")) %>%
  mutate(region = "Liver_GRO-seq") %>% distinct()
combined2 <- bind_rows(combined2, combined_gro)
combined3_gro <- read_feather(paste0(datapath, "/combined3_groseq.feather"))
combined3 <- bind_rows(combined3, combined3_gro)

combined_gropro <- read_feather(paste0(datapath, "/combined_groseq_promoter.feather")) %>%
  mutate(region = "Liver_GRO-seq_promoter") %>% distinct()
combined2 <- bind_rows(combined2, combined_gropro)

# hiding data
combined2 <- combined2 %>% filter(!(region %in% hide_region))
hide_index <- !(region_order %in% hide_region)
region_order <- region_order[hide_index]
region_main2 <- region_main2[hide_index[-c(1:length(region_main))]]
region_short <- region_short[hide_index]
region_one <- region_one[hide_index]

# name fix
combined2 <- combined2 %>% mutate(region = ifelse(str_detect(region, "GRO-seq"), region, str_c(region, "_RNA-seq")))
region_order <- ifelse(str_detect(region_order, "GRO-seq"), region_order, str_c(region_order, "_RNA-seq"))
region_main <- ifelse(str_detect(region_main, "GRO-seq"), region_main, str_c(region_main, "_RNA-seq"))
region_main2 <- ifelse(str_detect(region_main2, "GRO-seq"), region_main2, str_c(region_main2, "_RNA-seq"))

# read annotation file to find ucsc track
bed <- suppressWarnings(read_tsv(paste0(annotpath, "/final_tx_annotations_20200201.tsv.gz"),
  col_names = c(
    "chrom",
    "start",
    "end",
    "transcript_id",
    "score",
    "strand",
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    "gene_id",
    "unique_gene_symbol",
    NA,
    "clean_gene_symbol",
    NA,
    NA,
    "majiq_directed"
  ), skip = 1
)) %>%
  select(-contains("X")) %>%
  mutate(majiq_directed = factor(ifelse(is.na(majiq_directed), 0, 1)))

# read modules/clusters
mod <- read_feather(paste0(datapath, "/full_clusters.feather")) %>% select(gene, any_of(str_c("cluster_", region_short)))
mod <- mod[, c("gene", intersect(str_c("cluster_", region_short), colnames(mod)))] # edit for full

eigen <- suppressWarnings(read_tsv(paste0(datapath, "/cluster_patterns_matrices/reference_patterns.tsv"))) %>%
  rename(state = X1) %>%
  mutate(state = factor(state, levels = state_order))

eigen_gg <- list()
for (clu in colnames(eigen[, -1])) {
  df_plot <- data.frame(state = eigen$state, value = eigen[[clu]])
  eigen_gg[[clu]] <- ggplot(df_plot, aes(state, value, group = 1)) +
    geom_line() +
    geom_point(aes(color = state)) +
    scale_color_manual(values = state_cols) +
    ylab("expr") +
    theme(legend.position = "none")
}
eigen_gg[["filtered"]] <- ggplot(df_plot, aes(state, value, group = 1)) +
  ylab("expr") +
  theme(legend.position = "none")
eigen_gg[["insig"]] <- ggplot(df_plot, aes(state, value, group = 1)) +
  ylab("expr") +
  theme(legend.position = "none")
eigen_gg[["Unassigned"]] <- ggplot(df_plot, aes(state, value, group = 1)) +
  ylab("expr") +
  theme(legend.position = "none")

# query function
comb_fil_factor <- function(combined2, combined3, inid) {
  t1 <- Sys.time()

  combined3 <- combined3 %>% filter(unique_gene_symbol %in% inid)
  if (nrow(combined3) == 0) {
    combined3 <- combined3 %>% filter(gene_id %in% inid)
  }

  combined2 <- combined2 %>%
    filter(gene_id %in% (combined3$gene_id %>% unique())) %>%
    mutate(sample = (str_remove(sample, "[A-Z]+")))
  combined <- combined3 %>% inner_join(combined2, by = "gene_id")
  res <- combined %>% mutate(
    state = factor(state,
      levels = state_order
    ),
    region = factor(region,
      levels = region_order
    )
  )

  if (verbose_bench) {
    print(paste0("comb_fil_factor step: ", Sys.time() - t1))
  }
  res
}

# lists of genes
autocomplete_list <- str_sort(c(
  combined3$unique_gene_symbol,
  combined3$gene_id
) %>% unique(), decreasing = F)

find_spelling <- function(entry, dict) {
  temp <- intersect(
    c(entry, str_to_upper(entry), str_to_title(entry)),
    dict
  )
  if (length(temp) == 0) {
    return("")
  } else {
    temp[1]
  }
}

namedvec <- combined3$clean_gene_symbol
names(namedvec) <- combined3$unique_gene_symbol

unique_to_clean <- function(genevec, namedvec, na_omit = T) {
  temp <- namedvec[genevec]
  if (na_omit) {
    temp %>%
      na.omit() %>%
      unique()
  } else {
    temp
  }
}

# read go terms and TFs
gmt_lookup <- read_csv(paste0(annotpath, "/","gmt_id.csv"))
gmt_to_list <- function(path,
                        cutoff = 0,
                        sep = "\thttp://www\\..*?.org/gsea/msigdb/cards/.*?\t",
                        rm = "REACTOME_",
                        per = TRUE) {
  df <- readr::read_csv(path,
    col_names = F, col_types = cols()
  )
  df <- tidyr::separate(df,
    X1,
    sep = sep,
    into = c("path", "genes")
  ) %>% dplyr::mutate(path = stringr::str_remove(path, rm))
  if (per) {
    df <- df %>% mutate(
      genes = stringr::str_split(
        genes,
        "\t"
      )
    )
    return(tidyr::unnest_legacy(df, genes))
  } else {
    l <- str_split(df$genes, "\t")
    names(l) <- df$path
    return(l)
  }
}

gmt <- gmt_to_list(paste0(annotpath, "/", gmt_file), rm = gmt_short)

construct_gmtlist <- function(gmt_file, genes, n) {
  temp <- list()
  temp[["n"]] <- n

  temp[["Biological Process"]] <- gmt_to_list(paste0(annotpath, "/", str_replace(gmt_file, "all", "bp")), rm = gmt_short, per = FALSE)
  temp[["Biological Process"]] <- sapply(temp[["Biological Process"]], function(x) {
    intersect(x, genes)
  })
  temp[["Biological Process"]] <- temp[["Biological Process"]][sapply(temp[["Biological Process"]], length) >= 5]

  temp[["Cellular Component"]] <- gmt_to_list(paste0(annotpath, "/", str_replace(gmt_file, "all", "cc")), rm = gmt_short, per = FALSE)
  temp[["Cellular Component"]] <- sapply(temp[["Cellular Component"]], function(x) {
    intersect(x, genes)
  })
  temp[["Cellular Component"]] <- temp[["Cellular Component"]][sapply(temp[["Cellular Component"]], length) >= 5]

  temp[["Molecular Function"]] <- gmt_to_list(paste0(annotpath, "/", str_replace(gmt_file, "all", "mf")), rm = gmt_short, per = FALSE)
  temp[["Molecular Function"]] <- sapply(temp[["Molecular Function"]], function(x) {
    intersect(x, genes)
  })
  temp[["Molecular Function"]] <- temp[["Molecular Function"]][sapply(temp[["Molecular Function"]], length) >= 5]

  temp
}

if (file.exists(paste0(annotpath, "/", gmt_file, "_br.rds"))) {
  gmtlist_br <- readRDS(paste0(annotpath, "/", gmt_file, "_br.rds"))
} else {
  gmtlist_br <- construct_gmtlist(
    gmt_file,
    str_to_upper(combined3$clean_gene_symbol %>% unique()),
    length(combined3$clean_gene_symbol %>% unique())
  )
  saveRDS(gmtlist_br, paste0(annotpath, "/", gmt_file, "_br.rds"))
}

if (file.exists(paste0(annotpath, "/", gmt_file, ".rds"))) {
  gmtlist_sq <- readRDS(paste0(annotpath, "/", gmt_file, ".rds"))
} else {
  gmtlist_sq <- construct_gmtlist(
    gmt_file,
    str_to_upper(bed$clean_gene_symbol %>% unique()),
    length(bed$clean_gene_symbol %>% unique())
  )
  saveRDS(gmtlist_sq, paste0(annotpath, "/", gmt_file, ".rds"))
}

if (file.exists(paste0(annotpath, "/", gmt_file, "_human.rds"))) {
  gmtlist_human <- readRDS(paste0(annotpath, "/", gmt_file, "_human.rds"))
} else {
  gmtlist_human <- construct_gmtlist(gmt_file, gmt$genes %>% unique(), 60662)
  saveRDS(gmtlist_human, paste0(annotpath, "/", gmt_file, "_human.rds"))
}

domains <- read_csv(paste0(datapath, "/novel_domains.csv"), col_types = "cc")

br_expr <- combined2 %>%
  filter(region %in% region_main) %>%
  pull(gene_id) %>%
  unique()

# load orf predictions
gro_padj <- read_feather(paste0(datapath, "/padj_groseq.feather"))
gropro_padj <- read_feather(paste0(datapath, "/padj_groseq_promoter.feather"))
orfs <- read_feather(paste0(datapath, "/full_padj_orf_cor.feather")) 
colnames(orfs) <- str_replace(colnames(orfs), "EAr", "Ar") # fix
orfs <- orfs %>%
  select(gene_id,
    orf_len = len,
    exons,
    rna_len = transcript,
    orf,
    unique_gene_symbol,
    everything()
  ) %>%
  left_join(gro_padj, by = "unique_gene_symbol") %>% 
  left_join(gropro_padj, by = "unique_gene_symbol") %>% 
  mutate(novel = factor(ifelse(str_detect(gene_id, "^G"), 1, 0))) %>%
  mutate(min_padj = select(., contains("LRT_padj")) %>% 
           reduce(pmin, na.rm = TRUE)) %>% 
  mutate(domains = factor(ifelse(gene_id %in% domains$gene_id, 1, 0))) %>%
  mutate(
    br_expr = factor(ifelse(gene_id %in% br_expr, 1, 0))
  ) %>%
  mutate(majiq_directed = factor(ifelse(is.na(majiq), 0, 1)))

# chromosome sizes
sq1 <- as.list(readRDS(paste0(annotpath, "/sq_chr.rds")))[1:chrlimit]
sq_g <- data.frame(
  chrom = names(sq1),
  size = unlist(sq1)
)
bed_fc <- readRDS(paste0(datapath, "/fc.rds"))
bed_fc_liv <- readRDS(paste0(datapath, "/fc_liv.rds")) 
bed_fc <- bind_rows(bed_fc, bed_fc_liv) %>% 
  bed_sort(by_chrom = T)

# editing sites
edits <- read_tsv(paste0(datapath, "/brain_editing_site_proportions.bed"))
edits$max <- apply(edits[, c(7:96)], 1, function(x) {
  max(x, na.rm = T)
})
edits <- edits %>%
  bed_intersect(bed, suffix = c("", "_gene")) %>%
  select(c("chrom", "start", "end", "site", "max"), unique_gene_symbol = unique_gene_symbol_gene) %>%
  distinct()
edits <- edits %>% bed_slop(sq_g, both = 2000000)
edited_genes <- edits$unique_gene_symbol %>% unique()
orfs <- orfs %>% mutate(edited = ifelse(unique_gene_symbol %in% edited_genes, TRUE, FALSE))

# load short orf predictions
sorf <- read_csv(paste0(datapath, "/MiPepid_pred.csv"))
sorf_blast <- read_csv(paste0(datapath, "/SmProt_blast.csv"))
fulltbl <- combined3 %>%
  select(-c(gene_symbol, clean_gene_symbol, original_gene_name)) %>%
  left_join(orfs %>% select(any_of(orf_cols_join), any_of(str_c(region_short, "_LRT_padj")), min_padj), by = "gene_id") %>%
  left_join(mod, by = c("unique_gene_symbol" = "gene")) %>%
  mutate(source = factor(source)) %>%
  mutate_at(vars(contains("cluster")), factor) %>%
  mutate_at(vars(contains("chrom")), factor)
fulltbl_collapse <- fulltbl %>%
  group_by(gene_id) %>%
  arrange(desc(orf_len), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  ungroup()
fulltbl <- fulltbl %>%
  mutate(micropeptide_pred = ifelse(transcript_id %in% sorf$transcript_DNA_sequence_ID, TRUE, FALSE)) %>%
  left_join(sorf_blast %>% select(transcript_id, micropeptide_homology = stitle), by = "transcript_id") %>%
  select(any_of(orf_cols_join), everything()) %>%
  distinct()
fulltbl_sorf <- fulltbl %>% filter(micropeptide_pred)
fulltbl_collapse <- fulltbl_collapse %>%
  mutate(micropeptide_pred = ifelse(gene_id %in% fulltbl_sorf$gene_id, TRUE, FALSE)) %>%
  left_join(fulltbl_sorf %>% select(gene_id, micropeptide_homology), by = "gene_id") %>%
  select(any_of(orf_cols_join), everything()) %>%
  distinct()

length_detected_genes <- orfs %>%
  filter(br_expr == 1) %>%
  pull(gene_id) %>%
  unique() %>%
  length()
orfs <- orfs %>% mutate(micropeptide_pred = ifelse(gene_id %in% fulltbl_sorf$gene_id, TRUE, FALSE))

# padj functions
find_padj <- function(region, state, tbl) {
  region
  temp <- str_c(tbl[str_remove(tbl, "_.+") == region_short[region_order == region] &
    str_detect(tbl, state)],
  collapse = "<br>"
  )
  if (length(temp) == 0) {
    temp <- "NA"
  }
  temp
}

sig_cols <- colnames(orfs)[colnames(orfs) %>% str_detect("vs")]
sig_sym <- data.frame(
  call = rep(1, length(sig_cols)), # place holder only
  row.names = sig_cols
)

calls_sig <- function(padj, sig_sym, pval) {
  t1 <- Sys.time()
  temp <- cbind(padj, sig_sym)
  temp <- temp %>%
    rownames_to_column("comp") %>%
    mutate(call1 = ifelse(padj <= pval, 1, 0))

  if (verbose_bench) {
    print(paste0("calls_sig: ", Sys.time() - t1))
  }
  temp
}

find_groups_igraph <- function(df) {
  df <- df[, !(names(df) == "region"), drop = F]
  g <- graph_from_data_frame(df, directed = FALSE)
  cg <- max_cliques(g)
  res <- lapply(cg, names)
  res
}

sort_groups <- function(groups, states, state_order) {
  t1 <- Sys.time()
  all_groups <- states
  leftout <- as.list(setdiff(all_groups, unlist(groups)))
  full <- c(groups, leftout)
  if ((length(full)) == length(states)) {
    full4 <- as.matrix(full) %>% t()
  } else {
    full4 <- sapply(full, function(x) x[order(match(x, state_order))])
    if (class(full4) != "matrix") {
      full4 <- sapply(full4, "length<-", max(lengths(full4)))
    }
  }
  full2 <- full4 %>%
    t() %>%
    as.data.frame()
  colnames(full2) <- str_c("V", 1:ncol(full2))
  full2 <- full2 %>%
    mutate_all(factor, levels = state_order)
  full2 <- full2[order(full2$V1, method = "radix"), , drop = F]

  full3 <- full2 %>%
    mutate(letter = letters[1:n()]) %>%
    gather(-letter, value = "state", key = "NA", na.rm = T)

  full3 <- full3[order(full3$letter, method = "radix"), , drop = F] %>%
    group_by(state) %>%
    summarize(letter = str_c(letter, collapse = ""))

  if (verbose_bench) {
    print(paste0("sort_groups: ", Sys.time() - t1))
  }
  full3
}

groups_to_letters_igraph <- function(df) {
  reg <- df$region %>% unique()
  g <- lapply(reg, function(x) {
    g <- df %>% filter(region == x)
    states <- unique(c(g$state1, g$state2))
    g <- g %>% filter(call1 == 0)
    g2 <- find_groups_igraph(g) #<<-?
    g3 <- sort_groups(g2, states, state_order)
    g3$region <- x
    return(g3)
  })
  do.call(rbind, g)
}

fisher <- function(genevec, gmtlist, length_detected_genes = NA, top = Inf) {
  gmtlist_vec <- unlist(gmtlist) %>% unique()
  genevec <- intersect(genevec, gmtlist_vec)
  sampleSize <- length(genevec)
  if (is.na(length_detected_genes)) {
    length_detected_genes <- gmtlist_vec %>% length()
  }
  res <- sapply(gmtlist, function(x) {
    if (length(x) == 0) {
      return(c(stringofhits = "", pval = 1))
    }
    hitInSample <- length(intersect(genevec, x))
    hitInPop <- length(x)
    failInPop <- length_detected_genes - hitInPop
    stringofhits <- intersect(genevec, x) %>% str_c(collapse = ",")
    if (length(stringofhits) == 0) {
      stringofhits <- ""
    }
    pval <- phyper(
      hitInSample - 1,
      hitInPop,
      failInPop,
      sampleSize,
      lower.tail = FALSE
    )
    return(c(hits = stringofhits, pval = pval))
  }, simplify = FALSE)

  res <- data.frame(res) %>% data.table::transpose()
  names(res) <- c("hits", "pval")
  res$pathway <- names(gmtlist)
  res$padj <- p.adjust(as.numeric(res$pval), method = "fdr")
  res %>%
    mutate(pval = as.numeric(pval)) %>%
    mutate(minuslog10 = -log10(padj)) %>%
    mutate(len = unlist(map(str_split(hits, ","), length))) %>%
    mutate(len = ifelse(hits == "", 0, len)) %>%
    mutate(go_len = unlist(lapply(gmtlist, length))) %>%
    arrange(desc(minuslog10), desc(len)) %>%
    select(pathway, pval, padj, minuslog10, pval, hits, len, go_len)
}

maj <- read_tsv(paste0(datapath, "/MAJIQ_dpsi_summary_sig_squirrelBox.tsv.gz")) %>%
  mutate(
    region = factor(region),
    comp = factor(comp)
  ) %>%
  rename(comp_pair = "comp") %>%
  distinct()

# seqs for kmer
seqs <- read_feather(paste0(datapath, "/utrs_sq.feather")) %>%
  filter(gene_id %in% combined3$gene_id)

if (file.exists(paste0(datapath, "/seqs_precal_noG.rds"))) {
  seqs_precal <- readRDS(paste0(datapath, "/seqs_precal_noG.rds"))
} else {
  seqs_precal <- list()
  seqs_precal[["5mers_utr3"]] <- generate_kmers(seqs %>% filter(
    str_length(utr3) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr3),
  k = 5
  )
  seqs_precal[["6mers_utr3"]] <- generate_kmers(seqs %>% filter(
    str_length(utr3) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr3),
  k = 6
  )
  seqs_precal[["7mers_utr3"]] <- generate_kmers(seqs %>% filter(
    str_length(utr3) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr3),
  k = 7
  )
  seqs_precal[["5mers_utr5"]] <- generate_kmers(seqs %>% filter(
    str_length(utr5) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr5),
  k = 5
  )
  seqs_precal[["6mers_utr5"]] <- generate_kmers(seqs %>% filter(
    str_length(utr5) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr5),
  k = 6
  )
  seqs_precal[["7mers_utr5"]] <- generate_kmers(seqs %>% filter(
    str_length(utr5) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr5),
  k = 7
  )
  saveRDS(seqs_precal, paste0(datapath, "/seqs_precal_noG.rds"))
}

comp_kmer <- function(df = seqs,
                      gene_vec,
                      col = "utr3",
                      bac = seqs_precal[["5mers_utr3"]],
                      k = 5,
                      cutoff = 200,
                      recal_bac = FALSE) {
  enq <- df %>%
    filter(str_to_upper(unique_gene_symbol) %in% str_to_upper(gene_vec)) %>%
    pull(col) %>%
    na.omit()
  if (length(enq) == 0) {
    return(NA)
  }
  enq <- enq[str_length(enq) >= cutoff]
  enq_res <- generate_kmers(enq, k)

  if (recal_bac) {
    bac <- df %>%
      filter(!(str_to_upper(unique_gene_symbol) %in% str_to_upper(gene_vec))) %>%
      pull(col) %>%
      na.omit()
    bac <- bac[str_length(bac) >= cutoff]
    bac <- generate_kmers(bac, k)
  }

  res <- compute_kmer_enrichment(enq_res,
    bac,
    permutation = FALSE,
    chisq_p_value_threshold = 0,
    p_adjust_method = "fdr"
  )
  res$kmer <- str_replace_all(names(enq_res), "T", "U")
  res[order(res$adj_p_value, method = "radix"), , drop = F]
}

fivemers <- read_csv(paste0(annotpath, "/RBP_5mer.csv"))
sevenmers <- read_csv(paste0(annotpath, "/mir_7mer.csv"))
sixmers <- read_csv(paste0(annotpath, "/RBP_6mer.csv"))

# load curated gene lists
lists_vec <- list.files(listpath)
gene_list <- sapply(lists_vec, function(x) {
  read_csv(paste0(listpath, "/", x)) %>% pull(1)
}, simplify = FALSE)
names(gene_list) <- names(gene_list) %>% str_remove("\\..+")

# brief preview of each file
data_files <- list.files(path = datapath, full.names = T)
generate_prev <- function(x, ncol = 5, nrow = 3) {
  if (str_detect(x, ".bed") | str_detect(x, ".tsv")) {
    temp <- read_tsv(x, n_max = (nrow + 1))
    return(temp[1:min(nrow(temp), nrow), 1:min(ncol(temp), ncol), drop = FALSE])
  } else if (str_detect(x, ".feather")) {
    temp <- read_feather(x)
    return(temp[1:min(nrow(temp), nrow), 1:min(ncol(temp), ncol), drop = FALSE])
  } else if (str_detect(x, ".csv")) {
    temp <- read_csv(x, n_max = (nrow + 1))
    return(temp[1:min(nrow(temp), nrow), 1:min(ncol(temp), ncol), drop = FALSE])
  } else if (str_detect(x, ".rds")) {
    temp <- readRDS(x)
    if (!is.data.frame(temp)) {
      temp <- temp[[1]]
    }
    if (is.vector(temp)) {
      temp <- as.data.frame(temp)
    }
    return(temp[1:min(nrow(temp), nrow), 1:min(ncol(temp), ncol), drop = FALSE])
  } 
}
data_prev <- sapply(data_files, generate_prev, simplify = F)
names(data_prev) <- str_remove(names(data_prev), paste0(datapath, "/"))
                               
layout_ggplotly <- function(gg, x = -0.06, y = -0.12){
  gg[['x']][['layout']][['annotations']][[1]][['y']] <- y
  gg[['x']][['layout']][['annotations']][[2]][['x']] <- x
  gg
}  

proxy_height <- paste0(plot_height * 100 / 2, "px")

# check for qc report
if (!(file.exists(paste0("www/", qc_report)))) {
  qc_report <- NULL
}

# revigo
go <- read_tsv(paste0(datapath, "/ClusterGeneEnrichmentsLiver.txt")) %>% 
  mutate(cutsig = ifelse(Benjamini <= 0.001, 1, 0)) %>% 
  mutate(minuslogp = -log10(Benjamini))
bpall <- read_csv(paste0(datapath, "/REVIGO_BP_09.csv") )
bp <- read_csv(paste0(datapath, "/REVIGO_BP_09_clustered.csv") )