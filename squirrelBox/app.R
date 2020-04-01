library(dplyr)
library(tibble)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(feather)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(cowplot)
library(ComplexHeatmap)
library(igraph)
library(valr)
library(transite)
library(plotly)
library(crosstalk)
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinythemes)
library(shinycustomloader)
library(shinyjqui)
library(shinyWidgets)
library(rintrojs)


shinyOptions(cache = diskCache("./app-cache", max_size = 100 * 1024^2))
options(readr.num_columns = 0)
options(stringsAsFactors = FALSE)
options(spinner.type = 6)
theme_set(theme_cowplot())
# options(shiny.reactlog = TRUE)

### folders
datapath <- "data"
annotpath <- "annot"
listpath <- "data/lists"

### general data settings
versionN <- 0.97
geoN <- "G1234"
pageN <- 10
warningN <- 100
plot_width <- 8
plot_height <- 6
set_shinytheme <- "paper"
track_name <- "hub_1512849_KG_HiC"
track_url <- "https://squirrelhub.s3-us-west-1.amazonaws.com/hub/hub.txt"
gmt_file <- "c5.bp.v7.0.symbols.gmt" # add MF back in
gmt_short <- "GO_"
sig_cut <- 0.001
ncore <- parallel::detectCores() - 1

### choose and order columns
table_cols <- c(
  "gene_id", # comment out to remove from table outputs
  "unique_gene_symbol",
  "clean_gene_symbol",
  "original_gene_name",
  "source"
)
columns_tips <- c(
  "unique to each gene_id, from ensemble/NCBI or newly assigned",
  "gene clusters assigned by expression pattern for tissue",
  "gene clusters assigned by expression pattern for tissue",
  "gene clusters assigned by expression pattern for tissue",
  "adjusted pval for tissue, type 0...0.001 to filter for sig",
  "adjusted pval for tissue, type 0...0.001 to filter for sig",
  "adjusted pval for tissue, type 0...0.001 to filter for sig",
  "adjusted pval for tissue, type 0...0.001 to filter for sig",
  "adjusted pval for tissue, type 0...0.001 to filter for sig",
  "adjusted pval for tissue, type 0...0.001 to filter for sig",
  "newly annotated genes start with G",
  "reference denotes carry-over from ensemble/NCBI",
  "aa length of longest orf",
  "number of exons",
  "length of longest transcript",
  "newly annotated genes, gene_id start with G",
  "minimum adjusted pval across tissues",
  "whether domains were found for orf",
  "expressed above cutoff in at least one brain tissue",
  "expressed above cutoff in at least one non-brain tissue",
  "if collapsed, only longest transcript_id listed",
  "whether symbol assignment was assisted by majiq splic junctions"
) %>% str_c(collapse = "\', \'")
columns_tips <- paste0("\'", columns_tips, "\'")

orf_cols_join <- c(
  "gene_id",
  "orf_len",
  "exons",
  "rna_len",
  "novel",
  "min_sig",
  "domains",
  "br_expr",
  "nonbr_expr",
  "transcript_id",
  "majiq_directed"
)

orf_cols <- c(
  "gene_id",
  "transcript_id",
  "rna_len",
  "orf_len",
  "exons"
)

maj_cols <- c(
  "region",
  "comp_pair",
  "significant_90",
  "significant_99",
  "LSV_ID",
  "A5SS",
  "A3SS",
  "ES"
)

### sample settings, define state colors and order, region order
state_cols <- c(
  SA = rgb(255, 0, 0, maxColorValue = 255),
  IBA = rgb(67, 205, 128, maxColorValue = 255),
  Ent = rgb(155, 48, 255, maxColorValue = 255),
  LT = rgb(25, 25, 112, maxColorValue = 255),
  EAr = rgb(0, 0, 255, maxColorValue = 255),
  Ar = rgb(0, 0, 255, maxColorValue = 255),
  LAr = rgb(0, 102, 102, maxColorValue = 255),
  SpD = rgb(255, 165, 0, maxColorValue = 255)
)
state_order <- c(
  "SA",
  "IBA",
  "Ent",
  "LT",
  "EAr",
  "Ar",
  "LAr",
  "SpD"
)
region_order <- c(
  "Forebrain",
  "Hypothalamus",
  "Medulla",
  "Adrenal",
  "Kidney",
  "Liver"
)
region_main <- c(
  "Forebrain",
  "Hypothalamus",
  "Medulla"
)
region_main2 <- c(
  "Adrenal",
  "Kidney",
  "Liver"
)
region_short <- c(
  "fore",
  "hy",
  "med",
  "adr",
  "kid",
  "liv"
)
region_one <- c(
  "f",
  "h",
  "m",
  "a",
  "k",
  "l"
)

# read database
if (file.exists(paste0(datapath, "/combined2.feather"))) {
  combined2 <- read_feather(paste0(datapath, "/combined2.feather"))
  combined3 <- read_feather(paste0(datapath, "/combined3.feather"))
} else if (file.exists(paste0(datapath, "/combined2.csv"))) {
  combined2 <- fread(paste0(datapath, "/combined2.csv"), nThread = ncore)
  combined3 <- fread(paste0(datapath, "/combined3.csv"), nThread = ncore)
}

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
mod <- read_feather(paste0(datapath, "/clusters.feather"))
mod <- mod[, c("gene", intersect(str_c("cluster_", region_one), colnames(mod)))]

eigen <- read_tsv(paste0(datapath, "/cluster_patterns_matrices/reference_patterns.tsv")) %>%
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
  combined3 <- combined3 %>%
    filter(gene_id %in% inid | unique_gene_symbol %in% inid)
  combined2 <- combined2 %>%
    filter(gene_id %in% c(inid, combined3$gene_id %>%
      unique())) %>%
    mutate(sample = (str_remove(sample, "[A-Z]+")))
  combined <- combined3 %>% inner_join(combined2, by = "gene_id")
  combined %>% mutate(
    state = factor(state,
      levels = state_order
    ),
    region = factor(region,
      levels = region_order
    )
  )
}

# lists of genes
autocomplete_list <- str_sort(c(
  combined3$unique_gene_symbol,
  combined3$gene_id
) %>% unique(), decreasing = T)

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

unique_to_clean <- function(genevec, namedvec) {
  namedvec[genevec] %>% na.omit()
}

# read go terms and TFs
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

if (file.exists(paste0(annotpath, "/", gmt_file, ".rds"))) {
  gmtlist <- readRDS(paste0(annotpath, "/", gmt_file, ".rds"))
} else {
  gmtlist <- gmt_to_list(paste0(annotpath, "/", gmt_file), rm = gmt_short, per = FALSE)
  gmtlist <- sapply(gmtlist, function(x) {
    intersect(x, str_to_upper(bed$clean_gene_symbol %>% unique()))
  })
  gmtlist <- gmtlist[sapply(gmtlist, length) >= 5]
  saveRDS(gmtlist, paste0(annotpath, "/", gmt_file, ".rds"))
}

domains <- read_csv(paste0(datapath, "/novel_domains.csv"), col_types = "cc")

br_expr <- combined2 %>%
  filter(region %in% region_main) %>%
  pull(gene_id) %>%
  unique()

nonbr_expr <- combined2 %>%
  filter(region %in% region_main2) %>%
  pull(gene_id) %>%
  unique()

# load orf predictions
orfs <- read_feather(paste0(datapath, "/padj_orf.feather")) %>%
  select(gene_id,
    orf_len = len,
    exons,
    rna_len = transcript,
    orf,
    unique_gene_symbol,
    everything()
  ) %>%
  mutate(novel = factor(ifelse(str_detect(gene_id, "^G"), 1, 0))) %>%
  mutate(min_sig = pmin(hy_LRT_padj, med_LRT_padj, fore_LRT_padj, na.rm = T)) %>%
  mutate(domains = factor(ifelse(gene_id %in% domains$gene_id, 1, 0))) %>%
  mutate(
    br_expr = factor(ifelse(gene_id %in% br_expr, 1, 0)),
    nonbr_expr = factor(ifelse(gene_id %in% nonbr_expr, 1, 0))
  ) %>%
  mutate(majiq_directed = factor(ifelse(is.na(majiq), 0, 1)))

fulltbl <- combined3 %>%
  select(-c(gene_symbol, clean_gene_symbol, original_gene_name)) %>%
  left_join(orfs %>% select(orf_cols_join, contains("LRT")), by = "gene_id") %>%
  left_join(mod, by = c("unique_gene_symbol" = "gene")) %>%
  mutate(source = factor(source)) %>%
  mutate_at(vars(contains("cluster")), factor) %>%
  distinct()

fulltbl_collapse <- fulltbl %>%
  group_by(gene_id) %>%
  arrange(desc(orf_len), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  ungroup()

length_detected_genes <- orfs %>%
  filter(br_expr == 1) %>%
  pull(gene_id) %>%
  unique() %>%
  length()

# padj functions
find_padj <- function(region, state, tbl) {
  temp <- str_c(tbl[str_sub(tbl, 1, 1) == str_to_lower(str_sub(region, 1, 1)) &
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

calls_sig <- function(padj, sig_sym) {
  temp <- cbind(padj, sig_sym)
  temp <- temp %>%
    rownames_to_column("comp") %>%
    mutate(call1 = ifelse(padj <= sig_cut, 1, 0))
  temp
}

find_groups_igraph <- function(df) {
  df <- df %>% select(-region)
  g <- graph_from_data_frame(df, directed = FALSE)
  cg <- max_cliques(g)
  lapply(cg, names)
}

sort_groups <- function(groups, states, state_order) {
  all_groups <- states
  leftout <- list(setdiff(all_groups, unlist(groups)))
  leftout <- as.list(unlist(leftout))
  full <- c(groups, leftout)
  full4 <- sapply(full, function(x) factor(x)[order(factor(x, levels = state_order))])
  if (class(full4) != "matrix") {
    full4 <- sapply(full4, "length<-", max(lengths(full4)))
  }
  full2 <- full4 %>%
    t() %>%
    as.data.frame()
  colnames(full2) <- str_c("V", 1:ncol(full2))
  full2 <- full2 %>%
    mutate_all(factor, levels = state_order) %>%
    arrange(V1)
  full3 <- full2 %>%
    mutate(letter = letters[1:n()]) %>%
    pivot_longer(-letter, names_to = "NA", values_to = "state") %>%
    filter(!(is.na(state))) %>%
    group_by(state) %>%
    summarize(letter = str_c(letter, collapse = "")) %>%
    ungroup()
  full3
}

groups_to_letters_igraph <- function(df) {
  reg <- df$region %>% unique()
  df2 <- df %>% filter(call1 == 0)
  g <- lapply(reg, function(x) {
    g <- df2 %>% filter(region == x)
    g2 <- find_groups_igraph(g)
    states <- unique(c(df2$state1, df2$state2))
    g3 <- sort_groups(g2, states, state_order)
    g3$region <- x
    return(g3)
  })
  do.call(rbind, g)
}

fisher <- function(genevec, gmtlist, length_detected_genes, top = Inf) {
  genevec <- intersect(genevec, unlist(gmtlist) %>% unique())
  sampleSize <- length(genevec)
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
    pval <- phyper(hitInSample - 1,
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
  seqs_precal[["5mers_utr3"]] <- generateKmers(seqs %>% filter(
    str_length(utr3) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr3),
  k = 5
  )
  seqs_precal[["6mers_utr3"]] <- generateKmers(seqs %>% filter(
    str_length(utr3) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr3),
  k = 6
  )
  seqs_precal[["7mers_utr3"]] <- generateKmers(seqs %>% filter(
    str_length(utr3) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr3),
  k = 7
  )
  seqs_precal[["5mers_utr5"]] <- generateKmers(seqs %>% filter(
    str_length(utr5) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr5),
  k = 5
  )
  seqs_precal[["6mers_utr5"]] <- generateKmers(seqs %>% filter(
    str_length(utr5) >= 200,
    str_length(cds) >= 200
  ) %>%
    pull(utr5),
  k = 6
  )
  seqs_precal[["7mers_utr5"]] <- generateKmers(seqs %>% filter(
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
  enq_res <- generateKmers(enq, k)

  if (recal_bac) {
    bac <- df %>%
      filter(!(str_to_upper(unique_gene_symbol) %in% str_to_upper(gene_vec))) %>%
      pull(col) %>%
      na.omit()
    bac <- bac[str_length(bac) >= cutoff]
    bac <- generateKmers(bac, k)
  }

  res <- computeKmerEnrichment(enq_res,
    bac,
    permutation = FALSE,
    chisq.p.value.threshold = 0,
    p.adjust.method = "fdr"
  )
  res$kmer <- str_replace_all(names(enq_res), "T", "U")
  res %>% arrange(adj.p.value)
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

# some other code for webpage functions
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'

# Define UI for application that draws the boxplot
ui <- fluidPage(
  title = "squirrelBox",
  theme = shinytheme(set_shinytheme),
  tags$style("
      .checkbox {
        line-height: 20px;
        margin-top: -15px;
        margin-bottom: -15px;
      }
      .header {
        position:fixed;
        z-index:100;
      }
  "),
  tags$head(tags$style(
    HTML(
      "
         #SIDE {
            background-color: #F8F8F8;
        }",
      ".btn-options {
         background-color: #F8F8F8;
        }"
    )
  )),
  useShinyjs(),
  introjsUI(),
  tags$style("
  .introjs-helperLayer {
  background: transparent;
}

.introjs-overlay {
  opacity: 0 !important;
  z-index: 99999999!important;
}

.introjs-helperLayer:before {
  opacity: 0;
  content: '';
  position: absolute;
  width: inherit;
  height: inherit;
  border-radius: 0.5em;
  border: .2em solid rgba(255, 217, 68, 0.8);
  box-shadow: 0 0 0 1000em rgba(0,0,0, .7);
  opacity: 1;
}

.introjs-helperLayer:after {
  content: '';
  left: 0;
  right: 0;
  top: 0;
  bottom: 0;
  position: absolute;
  z-index: 1000;
}
             "),
  tags$style(".mock {position:fixed;width:7%;margin-top: 60px;z-index:10;}"),
  tags$style("
             .download_this{
    margin-right: 1px;
}"),
  tags$head(tags$script(HTML(jscode))),
  tags$head(tags$style(
    type = "text/css",
    ".shiny-output-error { visibility: hidden; }",
    ".shiny-output-error:before { visibility: hidden; }",
    ".shiny-output-warning { visibility: hidden; }",
    ".shiny-output-warning:before { visibility: hidden; }"
  )),
  tags$head(tags$style(
    type = "text/css",
    "#ucscPlot img {max-width: 100%; width: 100%; height: auto}"
  )),
  titlePanel(div(
    class = "header", img(src = "logo.png", style = "width : 4%;"),
    "13-lined ground squirrel gene-level RNA-seq expression",
    style = "font-size:23px"
  )),
  fixedPanel(
    style = "z-index:100;",
    dropdownButton(
      circle = TRUE, icon = icon("gear"), status = "options",
      size = "sm", up = TRUE, inline = TRUE,
      tooltip = tooltipOptions(
        title = "interface options and tips",
        placement = "top"
      ), margin = "20px",
      br(),
      actionButton("tutorial", "show me around"),
      br(),
      div(id = "doLockdiv", checkboxInput("doLock", "lock main panel order", value = F, width = NULL)),
      div(id = "doBrdiv", checkboxInput("doBr", "plot brain data", value = T, width = NULL)),
      div(id = "doTisdiv", checkboxInput("doTis", "plot non-brain data", value = F, width = NULL))
    ),
    actionButton("back_to_top", label = "to_top"),
    bsButton("showpanel", "sidebar", type = "toggle", value = FALSE),
    right = 10,
    bottom = 10,
    bsTooltip("back_to_top", "reset to top of page", placement = "top"),
    bsTooltip("showpanel", "toggle sidebar display", placement = "top")
  ),
  sidebarLayout(
    sidebarPanel(
      id = "SIDE",
      style = "position:fixed;width:23%;margin-top: 60px;z-index:50;",
      width = 3,
      div(
        id = "sideall",
        div(
          style = "display: inline-block;vertical-align:top; width: 160px;",
          tagAppendAttributes(selectizeInput("geneID",
            label = NULL,
            selected = "",
            choices = ""
          ),
          `data-proxy-click` = "Find"
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 10px;",
          actionButton("Find", "Find", icon = icon("search"))
        ),
        bsTooltip("Find", "gene id/symbols accepted"),
        div(
          id = "sidediv",
          tabsetPanel(
            id = "side1",
            tabPanel(
              span(icon("link", class = NULL, lib = "font-awesome"), "links", title = "save files and external links"),
              fluidRow(
                column(
                  width = 1
                ),
                column(
                  width = 4,
                  uiOutput("tab"),
                  uiOutput("blastlink")
                ),
                column(
                  width = 1
                ),
                column(
                  width = 4,
                  uiOutput("tab2"),
                  uiOutput("tab3"),
                  uiOutput("tab4")
                )
              ),
              downloadButton("savePlot", label = "save plot", class = "download_this"),
              downloadButton(outputId = "saveTable", label = "save table", class = "download_this"),
              br(),
              bsTooltip("savePlot", "save plot as pdf"),
              bsTooltip("saveTable", "save filtered/result table as csv"),
            ),
            tabPanel(
              span("options", title = "options specific to each main tab"),
              br(.noWS = "outside"),
              # div(id = "doPlotlydiv", checkboxInput("doPlotly", "interactive plots", value = F, width = NULL)),
              # div(id = "doPadjdiv", checkboxInput("doPadj", "indicate sig", value = T, width = NULL)),
              # div(id = "doNamediv", checkboxInput("doName", "additional labels", value = F, width = NULL)),
              # div(id = "doBrdiv", checkboxInput("doBr", "plot brain data", value = T, width = NULL)),
              # div(id = "doTisdiv",checkboxInput("doTis", "plot non-brain data", value = F, width = NULL)),
              div(id = "doEigendiv", checkboxInput("doEigen", "plot model clusters", value = T, width = NULL)),
              checkboxInput("doUcsc", "pull track", value = T, width = NULL),
              checkboxInput("doMod", "find module", value = T, width = NULL),
              checkboxInput("doKegg", "GO terms", value = T, width = NULL),
              # div(id = "doNormdiv", checkboxInput("doNorm", "line plot norm to SA", value = F, width = NULL)),
              checkboxInput("doTooltips", "show hover tips", value = T, width = NULL),
              bsTooltip("doPlotlydiv", "toggles main and kmer plot",
                options = list(container = "body"), placement = "right"
              ),
              bsTooltip("doPadjdiv", str_c("label groups by p <= ", sig_cut),
                options = list(container = "body"), placement = "right"
              ),
              bsTooltip("doNamediv", "label points by sample", options = list(container = "body"), placement = "right"),
              bsTooltip("doName2div", "show toggleable legend", options = list(container = "body"), placement = "right"),
              bsTooltip("doEigendiv", "assigned per region"),
              bsTooltip("doNormdiv", "otherwise centered by mean", options = list(container = "body"), placement = "right"),
              bsTooltip("doLockdiv", "by default panels can be dragged and rearranged")
            )
          ),
          tabsetPanel(
            id = "side2",
            tabPanel(
              span(icon("file-alt", class = NULL, lib = "font-awesome"), "load", title = "load list of genes for analysis from file or interactive table"),
              div(id = "filediv", fileInput("file", label = NULL)),
              div(
                uiOutput("listn"),
                style = "display: inline-block;vertical-align:middle;"
              ),
              div(
                style = "display: inline-block;float:right;vertical-align:middle;",
                disabled(actionButton("Prev1", "Prev", icon = icon("angle-up"))),
                disabled(actionButton("Next1", "Next", icon = icon("angle-down")))
              ),
              bsTooltip("filediv", "expects gene symbols as first column, or comma separated"),
              bsTooltip("Prev1", "query previous gene on loaded list"),
              bsTooltip("Next1", "query next gene on loaded list"),
              DT::dataTableOutput("tbllist"),
              style = "height:300px; overflow-y: scroll;"
            ),
            tabPanel(
              value = "cart",
              introBox(
                span(icon("shopping-cart", class = NULL, lib = "font-awesome"), "cart", title = "cart list of genes to save and export"),
                data.step = 1,
                data.intro = "add genes to a cart list, which can be exported or moved to the loaded list"
              ),
              uiOutput("listn2"),
              actionButton("Add", "Add"),
              actionButton("Load", "Load"),
              downloadButton(
                outputId = "saveList",
                label = "export"
              ),
              bsTooltip("Add", "add current query gene to cart"),
              bsTooltip("Load", "send to loaded list in side panel"),
              DT::dataTableOutput("tbllist2"),
              style = "height:300px; overflow-y: scroll;"
            ),
            tabPanel(
              introBox(
                span(icon("history", class = NULL, lib = "font-awesome"), "history", title = "history list of query genes"),
                data.step = 2,
                data.intro = "list of genes searched",
                data.position = "right"
              ),
              DT::dataTableOutput("historyl"),
              style = "height:300px; overflow-y: scroll;"
            )
          )
        )
      )
    ),
    mainPanel(
      id = "MAIN",
      width = 9,
      style = "z-index:1;margin-top: 60px;",
      tabsetPanel(
        id = "tabMain",
        tabPanel(
          title = span(icon("pencil-ruler", class = NULL, lib = "font-awesome"),
            "main",
            title = "Plot expression box plot and other info of query gene"
          ),
          value = "plot",
          div(
            id = "sorted",
            DT::dataTableOutput("results"),
            div(
              div(
                style = "display: inline-block;vertical-align:top;",
                dropdownButton(
                  circle = FALSE, status = "options", icon = icon("gear"), width = "200px", size = "sm",
                  tooltip = tooltipOptions(title = "plotting options"), margin = "20px",
                  br(),
                  div(id = "doPlotlydiv", checkboxInput("doPlotly", "interactive plots", value = F, width = NULL)),
                  div(id = "doPadjdiv", checkboxInput("doPadj", "indicate sig", value = T, width = NULL)),
                  div(id = "doNamediv", checkboxInput("doName", "additional labels", value = F, width = NULL))
                )
              ),
              div(
                style = "display: inline-block;vertical-align:top;",
                uiOutput("boxPlotUI") %>% withLoader(proxy.height = paste0(plot_height * 100 / 2, "px"))
              )
            ),
            bsCollapse(
              id = "tabs", multiple = TRUE, open = NULL,
              bsCollapsePanel(
                uiOutput("EigenPlot") %>% withLoader(),
                title = "cluster_assignments",
                style = "danger"
              )
            ),
            bsCollapse(
              id = "tabs2", multiple = TRUE, open = NULL,
              bsCollapsePanel(DT::dataTableOutput("orfinfo") %>% withLoader(),
                title = "called_orfs",
                style = "primary"
              )
            ),
            bsCollapse(
              id = "tabs3", multiple = TRUE, open = "NULL",
              bsCollapsePanel(
                DT::dataTableOutput("majinfo") %>% withLoader(),
                title = "majiq_alternative_splicing",
                style = "warning"
              )
            ),
            bsCollapse(
              id = "tabs4", multiple = TRUE, open = "NULL",
              bsCollapsePanel(htmlOutput("ucscPlot") %>% withLoader(),
                title = "UCSC browser plot",
                style = "success"
              )
            ),
            bsCollapse(
              id = "tabs5", multiple = TRUE, open = "NULL",
              bsCollapsePanel(DT::dataTableOutput("gotab") %>% withLoader(),
                title = "go_terms/domains",
                style = "info"
              )
            )
          )
        ),
        tabPanel(
          title = span(icon("table", class = NULL, lib = "font-awesome"),
            "transcript_gene",
            title = "Table of expression and other info of all genes/transcripts"
          ),
          value = "table_data",
          div(
            id = "doCollapsediv",
            style = "display: inline-block;width: 160px;",
            checkboxInput("doCollapse",
              "longest transcript",
              value = T,
              width = NULL
            )
          ),
          bsTooltip("doCollapsediv", "only show longest orf transcript for each gene"),
          actionButton("loadtab", "load"),
          bsTooltip("loadtab", "send to loaded list in side panel"),
          DT::dataTableOutput("tbl")
        ),
        tabPanel(
          title = span("majiq_alt",
            title = "Table of majiq output for alternative splicing events"
          ),
          value = "table_AS",
          div(
            id = "doJoindiv",
            style = "display: inline-block;width: 160px;",
            checkboxInput("doJoin",
              "gene info",
              value = FALSE,
              width = NULL
            )
          ),
          bsTooltip("doJoindiv", "bring in gene info as last columns"),
          actionButton("loadtab2", "load"),
          bsTooltip("loadtab2", "send to loaded list in side panel"),
          DT::dataTableOutput("alt")
        ),
        tabPanel(
          title = span("line_plot",
            title = "Plot expression of loaded gene list"
          ),
          value = "line_plot",
          div(
            style = "display: inline-block;vertical-align:top;",
            dropdownButton(
              circle = FALSE, status = "options", icon = icon("gear"), width = "200px", size = "sm",
              tooltip = tooltipOptions(title = "plotting options"), margin = "20px",
              br(),
              div(id = "doName2div", checkboxInput("doName2", "additional labels", value = F, width = NULL)),
              div(id = "doNormdiv", checkboxInput("doNorm", "line plot norm to SA", value = F, width = NULL))
            )
          ),
          div(
            style = "display: inline-block;vertical-align:top;",
            plotlyOutput("linePlot") %>% withLoader()
          )
        ),
        tabPanel(
          title = span("heatmap",
            title = "Plot Z-Score of loaded gene list as heat map"
          ),
          value = "heat_plot",
          br(),
          fluidRow(
            column(
              width = 3,
              checkboxInput("doRowcluster",
                "cluster rows",
                value = T,
                width = NULL
              ),
              checkboxInput("doColumncluster",
                "cluster columns",
                value = F,
                width = NULL
              )
            ),
            column(
              width = 3,
              checkboxInput("doSplit",
                "split by region",
                value = TRUE,
                width = NULL
              ),
              checkboxInput("doPivot",
                "pivot plot",
                value = F,
                width = NULL
              )
            ),
            column(
              width = 3,
              checkboxInput("doLabelgene",
                "label genes",
                value = T,
                width = NULL
              ),
              checkboxInput("doAutoresize",
                "resize on saving",
                value = F,
                width = NULL
              )
            )
          ),
          uiOutput("heatPlotUI") %>% withLoader()
        ),
        tabPanel(
          title = span("GO_enrichment",
            title = "GO term enrichment for loaded gene list (slow)"
          ),
          value = "enrichment_plot",
          plotlyOutput("richPlot") %>% withLoader()
        ),
        tabPanel(
          title = span("kmer",
            title = "kmer enrichment analysis and annotation for loaded gene list (slow)"
          ),
          value = "kmer_analysis",
          tags$style(HTML(".radio-inline {margin-left: 5px;margin-right: 25px;}")),
          div(
            style = "display: inline-block;vertical-align:top;",
            radioButtons("utr", "UTR choice", c("5UTR", "3UTR", "RBP"), selected = "3UTR", inline = TRUE)
          ),
          div(
            id = "kmerdiv",
            style = "display: inline-block;vertical-align:top; width:175px;",
            radioButtons("km", "kmer length", c("5", "6", "7"), selected = "6", inline = TRUE)
          ),
          div(
            id = "kmlabdiv",
            style = "display: inline-block;vertical-align:top; width:250px;",
            tags$style(HTML(".radio-inline {margin-right: 10px;}")),
            radioButtons("kmlab", "annotate kmer", c("RBP/mir", "seq", "none"), selected = "RBP/mir", inline = TRUE)
          ),
          div(
            style = "display: inline-block;vertical-align:top;",
            selectInput("utrlen", NULL, choices = c(200, 500, 1000, "full length"), selected = "full length")
          ),
          div(
            id = "rbptermdiv",
            style = "display: inline-block;vertical-align:top;width:250px",
            textInput("rbpterm", "highlight annotation", value = "MEX3C")
          ),
          div(
            id = "doPlotly2div",
            checkboxInput("doPlotly2", "interactive plot", value = F, width = NULL),
            style = "width:200px",
          ),
          uiOutput("kmerPlotUI") %>% withLoader(loader = "pacman"),
          bsTooltip(
            "kmerdiv",
            "Annotations: 5mer - Ray2013 + Encode, 6mer - Transite R, 7mer TargetScan mir seed"
          ),
          bsTooltip(
            "kmlabdiv",
            "label points with annotations"
          ),
          bsTooltip(
            "rbptermdiv",
            "highlights annotation in black"
          ),
          bsTooltip("doPlotly2div", "toggles main and kmer plot",
            options = list(container = "body")
          ),
        ),
        tabPanel(
          title = span("genes_venn",
            title = "visualize gene overlap between regions by venn diagram, and retrieve lists"
          ),
          value = "venn",
          div(
            style = "display: inline-block;vertical-align:top; width: 160px;",
            selectizeInput("seta", "geneset_A",
              choices = c(
                "_none", "_load_list", "_cart_list",
                names(gene_list)
              ),
              selected = "fore_sig"
            )
          ),
          div(
            style = "display: inline-block;vertical-align:top; width: 160px;",
            selectizeInput("setb", "geneset_B",
              choices = c(
                "_none", "_load_list", "_cart_list",
                names(gene_list)
              ),
              selected = "hy_sig"
            )
          ),
          div(
            style = "display: inline-block;vertical-align:top; width: 160px;",
            selectizeInput("setc", "geneset_C",
              choices = c(
                "_none", "_load_list", "_cart_list",
                names(gene_list)
              ),
              selected = "med_sig"
            )
          ),
          div(
            id = "doUpperdiv",
            checkboxInput("doUpper",
              "ignore case",
              value = T,
              width = NULL
            ),
            style = "display: inline-block; width: 100px;"
          ),
          div(
            id = "loadalldiv",
            actionButton("Cart_all", "Cart_all"),
            style = "display: inline-block"
          ),
          bsTooltip("doUpperdiv", "coerce all gene symbols to upper case"),
          bsTooltip("loadalldiv", "add all genes from these sets into `cart` side panel"),
          plotlyOutput("vennPlot") %>% withLoader()
        ),
        tabPanel(
          span(icon("question", class = NULL, lib = "font-awesome"),
            "about",
            title = "View version and author info",
          ),
          value = "about",
          uiOutput("intro"),
          uiOutput("track"),
          uiOutput("rawdata"),
          uiOutput("GOversion"),
          uiOutput("version"),
          uiOutput("GitHub"),
          uiOutput("contact"),
          column(width = 4, DT::dataTableOutput("explain")),
          column(width = 1),
          column(width = 4, DT::dataTableOutput("explain2"))
        )
      )
    )
  )
)

# Define server logic required to draw the boxplot and render metadata table
server <- function(input, output, session) {
  rv <- reactiveValues()
  rv$run2 <- 0
  rv$mod <- 0
  rv$conn <- 0
  rv$init <- 0
  rv$actuallygo <- 0
  rv$old <- ""
  rv$blast <- ""
  rv$pval <- data.frame()
  rv$temp_orfs <- data.frame()
  rv$listn <- 1
  rv$listn2 <- 0
  rv$cart <- 0
  rv$xsel <- "NA"
  rv$richsel <- "NA"
  rv$line <- 0
  rv$line_refresh <- 0
  rv$listn2renew <- 0
  rv$plot_temp <- data.frame()
  rv$mod_df <- data.frame()
  rv$toolarge <- 0
  rv$go <- 0
  rv$tabinit_plot <- 0
  rv$tabinit_data <- 0
  rv$tabinit_enrich <- 0
  rv$tabinit_kmer <- 0
  rv$tabinit_venn <- 0

  # hide some checkboxes
  removeModal()
  hide("doKegg")
  hide("doMod")
  hide("doUcsc")
  hide("doEigendiv")
  hide("doTooltips")
  hide("utrlen")
  hide("utr")

  # empty history list to start
  historytab <- c()

  # init
  observeEvent(rv$init == 0, {
    if (rv$init == 0) {
      query <- parseQueryString(session$clientData$url_search)
      if (!is.null(query[["gene"]])) {
        updateSelectizeInput(session,
          inputId = "geneID",
          selected = query[["gene"]],
          choices = autocomplete_list,
          server = T
        )
      } else {
        updateSelectizeInput(session,
          inputId = "geneID",
          selected = "Mex3c",
          choices = autocomplete_list,
          server = T
        )
      }
      rv$init <- 1
      rv$run2 <- 1
    }
  })

  # sortable or not
  jqui_sortable(ui = "#sorted", operation = "enable", options = list(cancel = ".datatables"))
  jqui_sortable(ui = "#sidediv", operation = "enable", options = list(cancel = ".datatables"))
  jqui_sortable(ui = "#tabMain", operation = "enable")

  observe({
    if (input$doLock == TRUE) {
      jqui_sortable(ui = "#sorted", operation = "destroy")
      jqui_sortable(ui = "#sidediv", operation = "destroy")
      jqui_sortable(ui = "#tabMain", operation = "destroy")
    } else if (input$doLock != TRUE) {
      jqui_sortable(ui = "#sorted", operation = "enable", options = list(cancel = ".datatables"))
      jqui_sortable(ui = "#sidediv", operation = "enable", options = list(cancel = ".datatables"))
      jqui_sortable(ui = "#tabMain", operation = "enable")
    }
  })

  # jump to plot
  observeEvent(input$Find, {
    if ((input$geneID != "") & (input$geneID != rv$old)) {
      updateTabsetPanel(session,
        "tabMain",
        selected = "plot"
      )
      rv$actuallygo <- rv$actuallygo + 1
      if (rv$init == 1) {
        rv$init <<- 2
      } else {
        updateQueryString(paste0("?gene=", input$geneID), mode = "push")
      }
    }
  })

  # query
  inid <- eventReactive(rv$actuallygo,
    {
      rv$old <<- input$geneID
      shinyjs::runjs("window.scrollTo(0, 0)")
      input$geneID
    },
    ignoreNULL = T,
    ignoreInit = T
  )

  historytab <- c()
  historytablist <- c()
  carttablist <- c()

  # boxplot1
  boxPlot1 <- reactive({
    plot_temp <- rv$plot_temp
    if (nrow(plot_temp) == 0) {
      return(ggplot())
    }

    if (input$doTis & input$doBr) {
      mis <- setdiff(region_order, plot_temp$region %>% unique() %>% as.character())
    } else if (!(input$doTis) & input$doBr) {
      mis <- setdiff(region_main, plot_temp$region %>% unique() %>% as.character())
    } else {
      mis <- setdiff(region_main2, plot_temp$region %>% unique() %>% as.character())
    }

    if (length(mis) > 0) {
      for (element in mis) {
        l <- as.list(plot_temp[1, ])
        l$region <- element
        l$sample <- "0"
        l$log2_counts <- NA
        l$state <- state_order[1]
        plot_temp <- rbind(plot_temp, l)
      }
    }
    if (!input$doTis) {
      plot_temp <- plot_temp %>% filter(!(region %in% region_main2))
    }
    if (!input$doBr) {
      plot_temp <- plot_temp %>% filter(!(region %in% region_main))
    }
    if (nrow(rv$pval) != 0) {
      padj <- rv$pval %>%
        select(ends_with("_wald_padj")) %>%
        t()
      if (ncol(padj) > 1) {
        padj <- padj[, 1, drop = FALSE]
      }
      tbl <- str_c(rownames(padj), format(padj[, 1], digits = 2), sep = " : ")
      plot_temp <- plot_temp %>% mutate(text = mapply(
        find_padj,
        as.character(region),
        as.character(state),
        list(tbl)
      ))
    } else {
      plot_temp <- plot_temp %>% mutate(text = "NA")
    }

    set.seed(1)
    g <- ggplot(plot_temp, aes(state, log2_counts, text = text)) +
      ylab("rlog(counts)") +
      facet_wrap(~region, scales = "free") +
      theme(legend.position = "none")

    if (input$doName == T) {
      g <- g +
        geom_boxplot(aes(color = state), outlier.shape = NA) +
        scale_color_manual(values = state_cols) +
        geom_text(position = position_jitter(seed = 1), aes(label = sample))
    } else {
      g <- g +
        geom_boxplot(aes(fill = state), outlier.shape = NA) +
        scale_fill_manual(values = state_cols) +
        geom_point(position = position_jitter(seed = 1))
    }

    if (input$doPadj == T & nrow(rv$pval) != 0 & input$doPlotly == F) {
      temp2 <- calls_sig(padj, sig_sym)
      temp2 <- temp2 %>%
        replace_na(list(call1 = list(0))) %>%
        separate(comp, into = c("region", "state1", NA, "state2"), extra = "drop") %>%
        select(-padj, -call) %>%
        mutate(call1 = as.numeric(call1)) %>%
        mutate(region = as.character(region))
      temp2$region <- region_order[factor(temp2$region, level = region_short) %>% as.numeric()]
      temp3 <- groups_to_letters_igraph(temp2) %>%
        mutate(region = factor(region, level = region_order))

      if (!(is.na(plot_temp$log2_counts) %>% all())) {
        agg <- aggregate(log2_counts ~ state + region, plot_temp, max)
        agg_min <- aggregate(log2_counts ~ state + region, plot_temp, min)
        agg$min <- agg_min$log2_counts
        agg2 <- agg %>%
          mutate(region = factor(region, level = region_order)) %>%
          group_by(region) %>%
          mutate(
            maxy = max(log2_counts),
            miny = min(min),
            nudgey = (maxy - miny) * 0.1
          ) %>%
          ungroup()
        agg3 <- agg2 %>%
          left_join(temp3 %>% select(region, state, letter)) %>%
          replace_na(list(letter = list("")))

        g <- g +
          geom_text(data = agg3, aes(
            text = letter,
            label = letter,
            y = log2_counts + nudgey,
            x = state,
            group = NULL
          ))
      }
    }

    g
  })

  # boxplot size
  boxPlotr <- reactive({
    g <- boxPlot1()
    output$boxPlot <- renderPlot(g)
    if (input$doTis + input$doBr == 2) {
      plotOutput("boxPlot", width = plot_width * 100, height = plot_height * 100)
    } else {
      plotOutput("boxPlot", width = plot_width * 100, height = plot_height * 100 / 2)
    }
  })

  # boxplot-plotly
  boxPlotlyr <- reactive({
    g <- boxPlot1()
    output$boxPlot2 <- renderPlotly(ggplotly(g + facet_wrap(~region), tooltip = "text") %>%
      layout(hovermode = "closest"))
    if (input$doTis + input$doBr == 2) {
      plotlyOutput("boxPlot2", width = plot_width * 100, height = plot_height * 100)
    } else {
      plotlyOutput("boxPlot2", width = plot_width * 100, height = plot_height * 100 / 2)
    }
  })

  # actually draw boxplot
  output$boxPlotUI <- renderUI({
    if (input$doPlotly == FALSE) {
      boxPlotr()
    } else {
      boxPlotlyr()
    }
  })

  # actually draw model module/cluster
  output$EigenPlot <- renderUI({
    if (input$doEigen != T) {
      plotOutput("boxPlot3", height = 1)
    } else {
      plotOutput("boxPlot3", width = plot_width * 100, height = plot_height * 100 / 2)
    }
  })

  # boxplot - models
  output$boxPlot3 <- renderCachedPlot(
    {
      if (input$doEigen != T) {
        g <- ""
        g
      } else {
        if (nrow(rv$mod_df) == 0) {
          mods <- c("filtered", "filtered", "filtered")
          reg <- colnames(mod)[-1]
        } else {
          mods <- (rv$mod_df[1, ] %>% unlist())[-1]
          reg <- names(rv$mod_df)[-1]
        }
        cowplot::plot_grid(
          plotlist = map(mods, function(x) eigen_gg[[x]]),
          labels = str_c(str_remove(reg, "cluster_"), mods, sep = ": "),
          ncol = 3,
          label_x = .3, hjust = 0
        )
      }
    },
    cacheKeyExpr = {
      rv$mod_df %>% select(-1)
    },
    sizePolicy = sizeGrowthRatio(width = plot_width * 100, height = plot_height * 100 / 2, growthRate = 1.2)
  )

  # filter data
  outputtab <- reactive({
    inid <- inid()
    temp_orfs <- orfs %>% filter((gene_id == inid) | (unique_gene_symbol == inid))
    if (nrow(temp_orfs) > 0) {
      if (nrow(temp_orfs) == 0) {
        rv$blast <<- ""
        rv$pval <<- data.frame()
        rv$temp_orfs <<- data.frame()
      } else {
        rv$pval <<- temp_orfs
        rv$blast <<- temp_orfs$orf[1]
        rv$temp_orfs <<- temp_orfs
      }
    } else {
      rv$blast <<- ""
      rv$pval <<- data.frame()
      rv$temp_orfs <<- data.frame()
    }

    rv$plot_temp <<- comb_fil_factor(combined2, combined3, inid)

    filtered <- rv$plot_temp %>%
      select(table_cols) %>%
      unique()

    if (nrow(filtered) == 0) {
      rv$blast <<- ""
      rv$pval <<- data.frame()
      rv$temp_orfs <<- data.frame()
      return(NULL)
    }

    filtered2 <- bed %>%
      filter(gene_id == filtered$gene_id[1]) %>%
      select(c("chrom", "start", "end", "gene_id"))
    if (nrow(filtered2) > 1) {
      filtered2 <- cbind(bed_merge(filtered2), filtered2[1, "gene_id"])
    }

    if (nrow(filtered2 > 0)) {
      tempvec <- unique(c(filtered$unique_gene_symbol, historytab))
      if (length(tempvec) > 10) {
        tempvec <- tempvec[1:10]
      }
      historytab <<- tempvec
    }

    out <- cbind(filtered, filtered2)
    out <- out[, -which(duplicated(colnames(out)))]
    # clusters
    mod1 <- mod %>%
      filter(gene %in% out$unique_gene_symbol)
    if (length(mod1) == 0) {
      rv$mod_df <<- data.frame()
    } else {
      rv$mod_df <<- mod1
    }

    out
  })

  # display gene info
  output$results <- DT::renderDataTable({
    temp <- outputtab() %>% select(-clean_gene_symbol)
    DT::datatable(temp,
      class = "table-condensed",
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      options = list(
        searchable = FALSE,
        dom = "t",
        paging = FALSE,
        columnDefs = list(
          list(
            className = "dt-center",
            targets = 0:(ncol(temp) - 1)
          ),
          list(
            orderable = "false",
            targets = 0:(ncol(temp) - 1)
          )
        ),
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
          "}"
        )
      )
    )
  })

  # orf call table
  output$orfinfo <- DT::renderDataTable({
    outputtab()
    if (nrow(rv$temp_orfs) == 0) {
      temp <- data.frame(`no orf found` = "")
    } else {
      temp <- rv$temp_orfs %>% select(orf_cols)
    }
    DT::datatable(temp,
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      options = list(searchable = FALSE, dom = "t")
    )
  })

  # majik report table
  output$majinfo <- DT::renderDataTable({
    temp <- maj %>% filter(unique_gene_symbol == outputtab()$unique_gene_symbol[1])
    if (nrow(temp) == 0) {
      temp <- data.frame(`no alternative splicing` = "")
    } else {
      temp <- temp %>% select(maj_cols)
    }
    DT::datatable(temp,
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      filter = "top",
      options = list(
        searchable = FALSE, dom = "t",
        autoWidth = FALSE
      )
    )
  })

  # goterm table
  output$gotab <- DT::renderDataTable({
    if (input$doKegg != T) {
      return()
    }
    outputtab <- outputtab()
    temp <- gmt %>% filter(genes == str_to_upper(outputtab$clean_gene_symbol))
    if (nrow(temp) == 0) {
      temp <- domains %>% filter(gene_id %in% outputtab$gene_id)
    }
    DT::datatable(temp,
      escape = FALSE,
      selection = "multiple",
      rownames = FALSE,
      options = list(
        dom = "ft", searchHighlight = TRUE,
        autoWidth = F
      )
    )
  })

  # download ucscplot
  output$ucscPlot <- renderUI({
    if (input$doUcsc != T) {
      return()
    }
    outputtab <- outputtab()
    if (nrow(outputtab) > 1) {
      outputtab$end <- max(outputtab$end)
      outputtab$start <- min(outputtab$start)
      outputtab <- outputtab[1, ]
    }
    l <- as.integer((outputtab$end - outputtab$start) * 0.2)
    url <- str_c(
      "http://genome.ucsc.edu/cgi-bin/hgTracks?db=",
      track_name,
      "&hubUrl=",
      track_url,
      "&ignoreCookie=1&position=",
      outputtab$chrom,
      ":",
      outputtab$start - l,
      "-",
      outputtab$end + l
    )
    src <- str_c(
      "http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=",
      track_name,
      "&hubUrl=",
      track_url,
      "&ignoreCookie=1&position=",
      outputtab$chrom,
      ":",
      outputtab$start - l,
      "-",
      outputtab$end + l
    )
    HTML(str_c(
      '<a href ="',
      url,
      '">',
      '<img src="',
      src,
      '">',
      "<a/>"
    ))
  })

  # various links in sidebar
  output$tab <- renderUI({
    outputtab <- outputtab()
    if (nrow(outputtab) > 1) {
      outputtab <- outputtab[1, ]
    }
    url <- a("trackhub",
      href = str_c(
        "http://genome.ucsc.edu/cgi-bin/hgTracks?db=",
        track_name,
        "&hubUrl=",
        track_url,
        "&ignoreCookie=1&position=",
        outputtab$chrom,
        ":",
        outputtab$start,
        "-",
        outputtab$end
      )
    )
    tagList(url)
  })

  # link genbank
  output$tab2 <- renderUI({
    outputtab <- outputtab()
    if (nrow(outputtab) > 1) {
      outputtab <- outputtab[1, ]
    }
    clean <- a("genbank",
      href = str_c(
        "https://www.ncbi.nlm.nih.gov/gene/?term=",
        str_remove(outputtab$unique_gene_symbol, "_.+"),
        "[sym]+AND+human[ORGN]"
      )
    )
    tagList(clean)
  })

  # link hgnc
  output$tab3 <- renderUI({
    outputtab <- outputtab()
    if (nrow(outputtab) > 1) {
      outputtab <- outputtab[1, ]
    }
    clean <- a("hgnc",
      href = str_c(
        "https://www.genenames.org/data/gene-symbol-report/#!/symbol/",
        str_remove(outputtab$unique_gene_symbol, "_.+")
      )
    )
    tagList(clean)
  })

  # link genecard
  output$tab4 <- renderUI({
    outputtab <- outputtab()
    if (nrow(outputtab) > 1) {
      outputtab <- outputtab[1, ]
    }
    clean <- a("genecard",
      href = str_c(
        "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
        str_remove(outputtab$unique_gene_symbol, "_.+")
      )
    )
    tagList(clean)
  })

  # link blast
  output$blastlink <- renderUI({
    if (rv$blast != "" & !(is.na(rv$blast))) {
      outputtab <- outputtab()
      if (nrow(outputtab) > 1) {
        outputtab <- outputtab[1, ]
      }
      orf <- rv$blast
      url <- a("blast_orf",
        href = str_c(
          "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastp&DATABASE=nr&QUERY=",
          orf
        )
      )
      tagList(url)
    } else {
      return()
    }
  })

  # save plot as pdf
  savePlot <- downloadHandler(
    filename = function() {
      outputtab <- outputtab()
      if (nrow(outputtab) > 1) {
        outputtab <- outputtab[1, ]
      }
      sym <- tryCatch(outputtab$unique_gene_symbol, error = function(err) {
        return("wrong")
      })
      paste0(sym, ".pdf", sep = "")
    },
    content = function(file) {
      if (input$doTis + input$doBr == 2) {
        w <- plot_width
        h <- plot_height
      } else {
        w <- plot_width
        h <- plot_height / 2
      }
      ggplot2::ggsave(file, plot = boxPlot1(), device = "pdf", width = w, height = h)
    }
  )

  # line plot
  linetemp <- reactive({
    inid()
    rv$line_refresh
    if (length(historytablist) == 0) {
      return(data.frame(unique_gene_symbol = "load data first"))
    }
    plot_temp <- comb_fil_factor(combined2, combined3, historytablist) %>%
      group_by(region, state, unique_gene_symbol) %>%
      summarize(counts = mean(2^log2_counts)) %>%
      ungroup()
    if (input$doNorm == TRUE) {
      plot_temp <- plot_temp %>%
        group_by(unique_gene_symbol, region) %>%
        mutate(log2_counts = log2(counts / counts[1])) %>%
        ungroup()
    } else {
      plot_temp <- plot_temp %>%
        group_by(unique_gene_symbol, region) %>%
        mutate(log2_counts = log2(counts / mean(counts))) %>%
        ungroup()
    }
    mis <- setdiff(region_main, plot_temp$region %>% unique() %>% as.character())
    if (length(mis) > 0) {
      for (element in mis) {
        l <- as.list(plot_temp[1, ])
        l$region <- element
        l$log2_counts <- NA
        l$state <- state_order[1]
        plot_temp <- rbind(plot_temp, l)
      }
    }
    if (!input$doTis) {
      plot_temp <- plot_temp %>% filter(!(region %in% region_main2))
    }
    if (!input$doBr) {
      plot_temp <- plot_temp %>% filter(!(region %in% region_main))
    }
    plot_temp
  })

  # plotly interative parts of line plot
  observeEvent(linetemp(), {
    rv$line <- rv$line + 1
    d <<- SharedData$new(linetemp, ~unique_gene_symbol, as.character(rv$line))
    d$clearSelection()
  })

  # actually draw line plot
  linePlot1 <- reactive({
    set.seed(1)
    linetemp()
    if (length(historytablist) == 0) {
      return(ggplotly(ggplot() +
        ggtitle("no genes loaded")))
    }

    g <- ggplot(d, aes(state, log2_counts,
      group = unique_gene_symbol,
      text = unique_gene_symbol
    )) +
      ylab("log2fold") +
      facet_wrap(~region) +
      theme(legend.position = "none") +
      geom_point(aes(color = unique_gene_symbol)) +
      geom_line(aes(color = unique_gene_symbol))
    g
  })

  output$linePlot <- renderPlotly({
    g <- linePlot1()
    fac <- input$doTis + input$doBr
    if ((rv$toolarge == 0) | (rv$toolarge == 1 & rv$go == 2)) {
      if (input$doName2 == T) {
        ggplotly(g, tooltip = "text", height = plot_height * 100 * fac / 2, width = plot_width * 100) %>%
          layout(
            autosize = FALSE,
            showlegend = TRUE
          )
      } else {
        ggplotly(g, tooltip = "text", height = plot_height * 100 * fac / 2, width = plot_width * 100)
      }
    } else {
      if (rv$go <= 0) {
        showModal(modalWarn)
      }
      ggplotly(ggplot() +
        ggtitle("plotting cancelled"))
    }
  })

  savePlot5 <- downloadHandler(
    filename = "lineplot.pdf",
    content = function(file) {
      fac <- input$doTis + input$doBr
      ggplot2::ggsave(file, plot = linePlot1(), device = "pdf", width = plot_width, height = plot_height / 2 * fac)
    }
  )

  # heatmap
  heatPlot1 <- reactive({
    set.seed(1)
    temp <- linetemp()
    if (length(historytablist) == 0) {
      return(Heatmap(matrix(), column_title = "no genes loaded", show_heatmap_legend = FALSE))
    }
    temp2 <- temp %>%
      select(-log2_counts) %>%
      pivot_wider(names_from = state, values_from = counts) %>%
      unite(region, unique_gene_symbol, col = "id", sep = ":") %>%
      column_to_rownames("id")
    if (input$doPivot) {
      temp2 <- scale(t(temp2))
    } else {
      temp2 <- t(scale(t(temp2)))
    }

    if (input$doSplit) {
      if (input$doPivot) {
        Heatmap(temp2,
          cluster_rows = input$doRowcluster,
          cluster_columns = input$doColumncluster,
          column_split = str_remove(colnames(temp2), ":.+"),
          column_labels = str_remove(colnames(temp2), "^.+:"),
          show_column_names = input$doLabelgene,
          heatmap_legend_param = list(title = "Z-Score")
        )
      } else {
        Heatmap(temp2,
          cluster_rows = input$doRowcluster,
          cluster_columns = input$doColumncluster,
          row_split = str_remove(rownames(temp2), ":.+"),
          row_labels = str_remove(rownames(temp2), "^.+:"),
          show_row_names = input$doLabelgene,
          heatmap_legend_param = list(title = "Z-Score")
        )
      }
    } else {
      if (input$doPivot) {
        Heatmap(temp2,
          cluster_rows = input$doRowcluster,
          cluster_columns = input$doColumncluster,
          heatmap_legend_param = list(title = "Z-Score"),
          show_column_names = input$doLabelgene
        )
      } else {
        Heatmap(temp2,
          cluster_rows = input$doRowcluster,
          cluster_columns = input$doColumncluster,
          heatmap_legend_param = list(title = "Z-Score"),
          show_row_names = input$doLabelgene
        )
      }
    }
  })

  heatPlot <- reactive({
    if ((rv$toolarge == 0) | (rv$toolarge == 1 & rv$go == 2)) {
      heatPlot1()
    } else {
      if (rv$go <= 0) {
        showModal(modalWarn)
      }
      Heatmap(matrix(), column_title = "plotting cancelled", show_heatmap_legend = FALSE)
    }
  })

  heatPlotr <- reactive({
    output$heatPlot2 <- renderPlot(heatPlot())
    if (input$doTis + input$doBr == 2) {
      plotOutput("heatPlot2", width = plot_width * 100, height = plot_height * 100 * 2)
    } else {
      plotOutput("heatPlot2", width = plot_width * 100, height = plot_height * 100)
    }
  })

  # actually draw boxplot
  output$heatPlotUI <- renderUI({
    heatPlotr()
  })

  savePlot3 <- downloadHandler(
    filename = "heatplot.pdf",
    content = function(file) {
      h <- heatPlot1()
      if (input$doAutoresize) {
        pdf(file, width = (h@matrix %>% dim())[2], height = (h@matrix %>% dim())[1])
      } else {
        pdf(file, width = plot_width, height = plot_height)
      }
      print(h)
      dev.off()
    }
  )

  richtemp <- reactive({
    rv$line_refresh
    set.seed(1)
    if (length(historytablist) == 0) {
      return(data.frame())
    }
    genevec <- unique_to_clean(historytablist, namedvec) %>% str_to_upper()
    tops <- fisher(genevec, gmtlist, length_detected_genes)
    tops <- tops %>% dplyr::slice(1:max(min(which(tops$padj > 0.01)), 15))
    tops$pathway <- reorder(tops$pathway, tops$minuslog10)
    tops
  })

  richPlot1 <- reactive({
    tops <- richtemp()
    if (nrow(tops) == 0) {
      return(ggplot() +
        ggtitle("no genes loaded"))
    }
    ggplot(
      tops %>% dplyr::slice(1:15),
      aes(x = pathway, y = minuslog10, fill = -minuslog10, text = len)
    ) +
      geom_bar(stat = "identity") +
      xlab(paste0("enriched : ", str_remove(gmt_short, "_"))) +
      coord_flip() +
      cowplot::theme_minimal_vgrid() +
      theme(
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.position = "none"
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      geom_hline(yintercept = -log10(0.05))
  })

  output$richPlot <- renderPlotly({
    p <- ggplotly(richPlot1(),
      source = "richPlot", tooltip = "text", height = plot_height * 100, width = plot_width * 100
    ) %>%
      layout(autosize = F) %>%
      highlight()
    event_register(p, "plotly_click")
    p
  })

  savePlot2 <- downloadHandler(
    filename = "enriched.pdf",
    content = function(file) {
      ggplot2::ggsave(file, plot = richPlot1(), device = "pdf", width = plot_width, height = plot_height)
    }
  )

  observeEvent(event_data("plotly_click", source = "richPlot"), {
    tops <- richtemp()
    rv$richsel <- event_data("plotly_click", source = "richPlot")
    gene_vec <- tops[16 - rv$richsel$y, ] %>%
      pull(hits) %>%
      str_split(",") %>%
      unlist()
    carttablist <<- gene_vec
    rv$listn2 <- length(carttablist)
    rv$listn2renew <- rv$listn2renew + 1
    updateTabsetPanel(session,
      "side2",
      selected = "cart"
    )
  })

  saveEnrich <- downloadHandler("enrich_list.csv", content = function(file) {
    write_csv(richtemp(), file)
  })

  # kmer analysis
  kmertemp <- reactive({
    rv$line_refresh
    set.seed(1)
    if (length(historytablist) == 0) {
      return(data.frame())
    }
    genevec <- unique_to_clean(historytablist, namedvec) %>% str_to_upper()
    if (input$utrlen == "full length") {
      lenchoice <- 0
    } else {
      lenchoice <- input$utrlen
    }
    if (input$utr == "RBP") {
      topsk <- comp_motif(gene_vec = genevec)
    } else {
      if (input$utr == "3UTR") {
        utrchoice <- "utr3"
      } else {
        utrchoice <- "utr5"
        lenchoice <- -lenchoice
      }
      precal <- paste0(input$km, "mers_", utrchoice)
      topsk <- comp_kmer(
        gene_vec = genevec,
        bac = seqs_precal[[precal]],
        col = utrchoice,
        k = as.numeric(input$km)
      )
      if (nrow(topsk) == 0 | is.null(topsk)) {
        return(data.frame())
      }
      if (input$km == "5") {
        topsk <- topsk %>% left_join(fivemers, by = c("kmer" = "fivemer"))
      } else if (input$km == "6") {
        topsk <- topsk %>% left_join(sixmers, by = c("kmer" = "hexamer"))
      } else {
        topsk <- topsk %>%
          left_join(sevenmers, by = c("kmer" = "sevenmer")) %>%
          rename(RBP = "mir")
      }
    }
    topsk <- topsk %>%
      mutate(minuslog10 = -log10(adj.p.value), enrichment = log2(enrichment))
    topsk %>% replace_na(list(RBP = ""))
  })

  kmerPlot1 <- reactive({
    topsk <- kmertemp()
    if (nrow(topsk) == 0) {
      return(ggplot() +
        ggtitle("no genes loaded"))
    }
    if (input$kmlab == "RBP/mir") {
      topsk <- topsk %>%
        mutate(sig = factor(ifelse(minuslog10 >= -log10(0.05), "sig", "insig"))) %>%
        mutate(text2 = ifelse((sig == "sig" & row_number() <= 15), str_c(kmer, RBP, sep = "\n"), "")) %>%
        mutate(text1 = str_c(kmer, RBP, sep = "\n"))
    } else {
      topsk <- topsk %>%
        mutate(sig = factor(ifelse(minuslog10 >= -log10(0.05), "sig", "insig"))) %>%
        mutate(text2 = ifelse((sig == "sig" & row_number() <= 15), kmer, "")) %>%
        mutate(text1 = kmer)
    }
    if (input$kmlab == "none") {
      topsk <- topsk %>% mutate(text2 = "")
    }

    ggplot(topsk, aes(
      x = enrichment,
      y = minuslog10,
      text = text1,
      text2 = RBP,
      label = text2
    )) +
      geom_point(aes(color = sig), size = 0.5) +
      scale_color_manual(values = c("#B3B3B3", "#FC8D62")) + # RColorBrewer::brewer.pal(8, "Set2")
      geom_point(data = topsk %>%
        filter(str_detect(RBP, input$rbpterm)), color = "black", size = 0.5) +
      ggrepel::geom_text_repel(box.padding = 0.05, size = 3, aes(label = text2)) +
      xlab("log2enrichment") +
      labs(color = "")
  })

  kmerPlotr <- reactive({
    output$kmerPlot <- renderPlot(kmerPlot1())
    plotOutput("kmerPlot", width = plot_width * 100, height = plot_height * 100)
  })

  kmerPlotlyr <- reactive({
    g <- kmerPlot1()
    g2 <- ggplotly(g,
      source = "kmerPlotly",
      tooltip = "text", height = plot_height * 100, width = plot_width * 100
    ) %>%
      layout(autosize = F)
    output$kmerPlot2 <- renderPlotly(g2)
    plotlyOutput("kmerPlot2", width = plot_width * 100, height = plot_height * 100)
  })

  output$kmerPlotUI <- renderUI({
    if (input$doPlotly2 == FALSE) {
      kmerPlotr()
    } else {
      kmerPlotlyr()
    }
  })

  savePlot4 <- downloadHandler(
    filename = "kmers.pdf",
    content = function(file) {
      ggplot2::ggsave(file, plot = kmerPlot1(), device = "pdf", width = plot_width, height = plot_height)
    }
  )

  saveK <- downloadHandler("enrich_kmer.csv", content = function(file) {
    write_csv(kmertemp(), file)
  })

  # venn plotly
  venntemp <- reactive({
    rv$line_refresh
    temp <- list(
      `Set A` = gene_list[[input$seta]],
      `Set B` = gene_list[[input$setb]],
      `Set C` = gene_list[[input$setc]]
    )
    load_vec <- c(input$seta, input$setb, input$setc) == "_load_list"
    if (sum(load_vec) > 0) {
      for (element in names(temp)[load_vec]) {
        temp[[element]] <- historytablist
      }
    }
    load_vec <- c(input$seta, input$setb, input$setc) == "_cart_list"
    if (sum(load_vec) > 0) {
      for (element in names(temp)[load_vec]) {
        temp[[element]] <- carttablist
      }
    }
    if (input$doUpper) {
      temp <- sapply(temp, function(x) str_to_upper(x))
    }

    if (!is.list(temp)) {
      temp <- list(temp)
      names(temp) <- "Set A"
    }
    temp
  })

  vennPlot1 <- reactive({
    a <- venntemp()
    non_none <- !sapply(a, is.null) & !duplicated(c(input$seta, input$setb, input$setc))
    if (sum(non_none) > 1) {
      g <- ggvenn::ggvenn(a, names(a)[non_none], show_elements = TRUE)
      g2 <- ggvenn::ggvenn(a, names(a)[non_none],
        show_percentage = FALSE,
        set_name_size = 4, stroke_size = 0, text_size = 6
      )
      g2$layers[[3]]$data$text <- c(input$seta, input$setb, input$setc)[non_none]
      g2$layers[[4]]$data$text2 <- g$layers[[4]]$data$text
      rv$venntext <<- g$layers[[4]]$data$text
      g2$labels[["text2"]] <- "text2"
      g2$layers[[4]]$mapping <- aes(
        x = x,
        y = y,
        label = text,
        hjust = hjust,
        vjust = vjust,
        text = text2
      )
      g2 + theme_cowplot() + theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank()
      )
    } else if (sum(non_none) == 1) {
      g <- ggvenn::ggvenn(a, c(names(a)[non_none], NA), show_elements = TRUE)
      g2 <- ggvenn::ggvenn(a, c(names(a)[non_none], NA),
        show_percentage = FALSE,
        set_name_size = 4, stroke_size = 0, text_size = 6
      )
      g2$data <- g$data %>% filter(group == "A")
      g2$layers[[3]]$data <- g2$layers[[3]]$data[1, ]
      g2$layers[[3]]$data$text <- c(input$seta, input$setb, input$setc)[non_none]
      g2$layers[[4]]$data <- g2$layers[[4]]$data[1, ]
      g2$layers[[4]]$data$text2 <- g$layers[[4]]$data$text[1]
      rv$venntext <<- g$layers[[4]]$data$text
      g2$labels[["text2"]] <- "text2"
      g2$layers[[4]]$mapping <- aes(
        x = x,
        y = y,
        label = text,
        hjust = hjust,
        vjust = vjust,
        text = text2
      )
      g2 + theme_cowplot() + theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank()
      )
    } else {
      ggplot() +
        ggtitle("no gene lists loaded")
    }
  })

  output$vennPlot <- renderPlotly({
    p <- ggplotly(vennPlot1(),
      source = "vennPlot",
      tooltip = "text2", height = plot_height * 100, width = plot_width * 100
    ) %>%
      layout(autosize = F, showlegend = FALSE) %>%
      highlight()
    event_register(p, "plotly_click")
    p
  })

  savePlot6 <- downloadHandler(
    filename = "venn.pdf",
    content = function(file) {
      ggplot2::ggsave(file, plot = vennPlot1(), device = "pdf", width = plot_width, height = plot_height)
    }
  )

  observeEvent(event_data("plotly_click", source = "vennPlot"), {
    aa <- event_data("plotly_click", source = "vennPlot")
    if (!is.null(aa$pointNumber)) {
      gene_string <- rv$venntext[aa$pointNumber + 1]
      gene_vec <- str_split(gene_string, ",")[[1]]
      carttablist <<- gene_vec
      rv$listn2 <- length(carttablist)
      updateTabsetPanel(session,
        "side2",
        selected = "cart"
      )
      rv$listn2renew <- rv$listn2renew + 1
    }
  })

  onclick("Cart_all", {
    gene_vec <- unlist(venntemp()) %>% unique()
    carttablist <- gene_vec
    rv$listn2 <- length(carttablist)
    updateTabsetPanel(session,
      "side2",
      selected = "cart"
    )
    rv$listn2renew <- rv$listn2renew + 1
  })

  # list history genes as table
  output$historyl <- DT::renderDataTable({
    inid()
    if (length(historytab) > 0) {
      DT::datatable(data.table::as.data.table(list(historytab)),
        escape = FALSE,
        selection = "single",
        rownames = FALSE,
        colnames = "genes",
        options = list(searchable = FALSE, dom = "t", paging = FALSE, scrollY = TRUE)
      )
    } else {
      DT::datatable(data.table::as.data.table(list(c(""))),
        escape = FALSE,
        selection = "single",
        rownames = FALSE,
        colnames = "genes",
        options = list(searchable = FALSE, dom = "t", paging = FALSE, scrollY = TRUE)
      )
    }
  })

  observeEvent(input$historyl_rows_selected, {
    rv$run2 <- 1
    sel <- find_spelling(historytab[input$historyl_rows_selected], autocomplete_list)
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = sel,
      choices = autocomplete_list,
      server = T
    )
  })

  # find function on new input
  observeEvent(input$geneID != "", {
    if (rv$run2 == 1 & input$geneID != "" & !is.null(input$geneID) & input$tabMain == "plot") {
      rv$run2 <- 1
      click("Find")
    }
  })

  # find on clicking button
  onclick(
    "geneID",
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = "",
      choices = autocomplete_list,
      server = T
    )
  )

  # loading list and viewing
  observeEvent(input$file, {
    rv$listn <- 0
    path <- input$file$datapath
    if (str_detect(input$file$name, "\\.tsv")) {
      v_genes <- read_tsv(path, col_names = FALSE)
    } else {
      v_genes <- read_csv(path, col_names = FALSE)
    }
    if (nrow(v_genes) == 1) {
      v_genes <- v_genes[1, ] %>%
        unlist() %>%
        str_split(", ", simplify = FALSE) %>%
        unlist() %>%
        str_split(",", simplify = FALSE) %>%
        unlist() %>%
        str_trim()
    } else if (nrow(v_genes) == 0) {
      v_genes <- ""
    } else {
      v_genes <- v_genes %>% pull(1)
    }
    historytablist <<- autocomplete_list[str_to_upper(autocomplete_list) %in%
      str_to_upper(v_genes %>% unique())]
    rv$line_refresh <- rv$line_refresh + 1
  })

  # cart list
  onclick("Add", {
    carttablist <- unique(c(historytab[1], carttablist))
    rv$listn2 <- length(carttablist)
  })

  output$listn2 <- renderUI({
    HTML(str_c("<strong><h5> # in cart: ", rv$listn2, "</h5></strong>"))
  })

  output$saveList <- downloadHandler("cart.txt", content = function(file) {
    write_lines(carttablist, file)
  })

  onclick("Load", {
    historytablist <- carttablist
    rv$line_refresh <- rv$line_refresh + 1
  })

  # list cart genes as table
  output$tbllist2 <- DT::renderDataTable({
    rv$listn2
    rv$listn2renew
    if (length(carttablist) > 0) {
      DT::datatable(data.table::as.data.table(list(carttablist)),
        escape = FALSE,
        selection = "single",
        rownames = FALSE,
        colnames = "genes",
        options = list(searchable = FALSE, dom = "t", paging = FALSE, scrollY = TRUE)
      )
    } else {
      DT::datatable(data.table::as.data.table(list(c(""))),
        escape = FALSE,
        selection = "single",
        rownames = FALSE,
        colnames = "genes",
        options = list(searchable = FALSE, dom = "t", paging = FALSE, scrollY = TRUE)
      )
    }
  })

  observeEvent(input$tbllist2_rows_selected, {
    rv$run2 <- 1
    sel <- find_spelling(carttablist[input$tbllist2_rows_selected], autocomplete_list)
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = sel,
      choices = autocomplete_list,
      server = T
    )
  })

  onclick("Prev1", {
    rv$listn <- rv$listn - 1
    if (rv$listn <= 1) {
      rv$listn <- 1
      disable("Prev1")
      if (rv$listn < length(historytablist)) {
        enable("Next1")
      }
    } else {
      enable("Prev1")
      enable("Next1")
    }
    rv$run2 <- 1
    sel <- find_spelling(historytablist[rv$listn], autocomplete_list)
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = sel,
      choices = autocomplete_list,
      server = T
    )
  })

  onclick("Next1", {
    rv$listn <- rv$listn + 1
    if (rv$listn >= length(historytablist)) {
      rv$listn <- length(historytablist)
      disable("Next1")
      if (rv$listn > 1) {
        enable("Prev1")
      }
    } else {
      enable("Next1")
      enable("Prev1")
    }
    rv$run2 <- 1
    sel <- find_spelling(historytablist[rv$listn], autocomplete_list)
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = sel,
      choices = autocomplete_list,
      server = T
    )
  })

  output$listn <- renderUI({
    rv$line_refresh
    HTML(paste0("<strong><h5>", rv$listn, " of ", length(historytablist), "</h5></strong>"))
  })

  # list loaded genes as table
  output$tbllist <- DT::renderDataTable({
    rv$line_refresh
    if (length(historytablist) > 0) {
      enable("Next1")
      if (length(historytablist) > warningN) {
        rv$toolarge <- 1
        rv$go <- 0
      } else {
        rv$toolarge <- 0
        rv$go <- 0
      }
      DT::datatable(data.table::as.data.table(list(historytablist)),
        escape = FALSE,
        selection = "single",
        rownames = FALSE,
        colnames = "genes",
        options = list(searchable = FALSE, dom = "t", paging = FALSE, scrollY = TRUE)
      )
    } else {
      disable("Prev1")
      disable("Next1")
      rv$toolarge <- 0
      rv$go <- 0
      DT::datatable(data.table::as.data.table(list(c(""))),
        escape = FALSE,
        selection = "single",
        rownames = FALSE,
        colnames = "genes",
        options = list(searchable = FALSE, dom = "t", paging = FALSE, scrollY = TRUE)
      )
    }
  })

  observeEvent(input$tbllist_rows_selected, {
    rv$run2 <- 1
    sel <- find_spelling(historytablist[input$tbllist_rows_selected], autocomplete_list)
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = sel,
      choices = autocomplete_list,
      server = T
    )
  })

  orftbl <- reactive({
    if (input$doCollapse) {
      fulltbl_collapse
    } else {
      fulltbl
    }
  })

  # explore bed table
  output$tbl <- DT::renderDataTable({
    temp <- orftbl() %>%
      select(
        unique_gene_symbol,
        contains("cluster"),
        contains("LRT"),
        everything()
      )
    clustercol <- which(str_detect(colnames(temp), "cluster")) - 1
    padjcol <- which(str_detect(colnames(temp), "padj"))
    DT::datatable(
      temp,
      filter = "top",
      escape = FALSE,
      selection = "single",
      rownames = FALSE,
      extensions = "ColReorder",
      options = list(
        # searchHighlight = TRUE, # breaks width
        pageLength = pageN,
        columnDefs = list(list(width = "125px", visible = TRUE, targets = as.list(clustercol), className = "dt-center")),
        scrollX = TRUE,
        autoWidth = TRUE,
        colReorder = list(enable = !input$doLock)
      ),
      callback = JS(paste0("var tips = [", columns_tips, "],
                            firstRow = $('#tbl thead tr th');
                            for (var i = 0; i < tips.length; i++) {
                              $(firstRow[i]).attr('title', tips[i]);
                            }"))
    ) %>%
      DT::formatRound(columns = padjcol, digits = 4)
  })

  saveFiltered <- downloadHandler("filtré.csv", content = function(file) {
    s <- input$tbl_rows_all
    write_csv((orftbl() %>%
      select(
        unique_gene_symbol,
        contains("cluster"),
        contains("LRT"),
        everything()
      ) %>% select(unique_gene_symbol, everything()))[s, ], file)
  })

  onclick("loadtab", {
    s <- input$tbl_rows_all
    historytablist <- orftbl()[s, ] %>% pull(unique_gene_symbol)
    rv$line_refresh <- rv$line_refresh + 1
  })

  observeEvent(input$tbl_rows_selected, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = orftbl()[input$tbl_rows_selected, "unique_gene_symbol"],
      choices = autocomplete_list,
      server = T
    )
  })

  # explore majiq table
  majtbl <- reactive({
    if (input$doJoin) {
      maj %>% left_join(fulltbl_collapse)
    } else {
      maj
    }
  })

  output$alt <- DT::renderDataTable({
    DT::datatable(
      majtbl() %>%
        select(
          unique_gene_symbol,
          contains("significant"),
          LSV_ID, A5SS, A3SS, ES,
          everything()
        ),
      filter = "top",
      escape = FALSE,
      selection = "single",
      rownames = FALSE,
      extensions = "ColReorder",
      options = list(
        searchHighlight = TRUE,
        pageLength = pageN,
        scrollX = TRUE,
        autoWidth = TRUE,
        colReorder = list(enable = !input$doLock)
      )
    )
  })

  saveFilteredAS <- downloadHandler("filtré.csv", content = function(file) {
    s <- input$alt_rows_all
    write_csv((majtbl() %>%
      select(
        unique_gene_symbol,
        contains("significant"),
        LSV_ID, A5SS, A3SS, ES,
        everything()
      ))[s, ], file)
  })

  onclick("loadtab2", {
    s <- input$alt_rows_all
    historytablist <- majtbl()[s, ] %>%
      pull(unique_gene_symbol) %>%
      unique()
    rv$line_refresh <- rv$line_refresh + 1
  })

  observeEvent(input$alt_rows_selected, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = majtbl()[input$alt_rows_selected, "unique_gene_symbol"],
      choices = autocomplete_list,
      server = T
    )
  })

  # back to top
  observeEvent(input$back_to_top,
    {
      shinyjs::runjs("window.scrollTo(0, 0)")
    },
    ignoreNULL = T
  )

  # ABOUT
  output$intro <- renderUI({
    url <- str_c("manuscript link here")
    clean <- a("manuscript",
      href = url
    )
    tagList(tags$h6("RNA sequencing data and analysis for 13-lined ground squirrel brain samples from hibernation cycle. Please see ", clean, "for more details."))
  })

  output$track <- renderUI({
    url <- str_c(
      "http://genome.ucsc.edu/cgi-bin/hgTracks?db=",
      track_name,
      "&hubUrl=",
      track_url
    )
    clean <- a("UCSC browser",
      href = url
    )
    tagList(tags$h6("Newly assembled genome and annotated transcriptome are hosted on ", clean))
  })

  output$rawdata <- renderUI({
    url <- str_c(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geoN
    )
    clean <- a(geoN,
      href = url
    )
    tagList(tags$h6("Raw data is deposited on GEO, ", clean))
  })

  output$GOversion <- renderUI({
    url <- "https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp"
    clean <- a(gmt_file,
      href = url
    )

    tagList(tags$h6("GO database version: ", clean))
  })

  output$version <- renderUI({
    url <- "https://github.com/rnabioco/squirrelbox/"
    clean <- a(versionN,
      href = url
    )
    tagList(tags$h6("squirrelBox version: ", clean))
  })

  output$explain <- DT::renderDataTable({
    dfreg <- data.frame(
      region = region_order,
      short = region_short,
      letter = region_one
    )
    DT::datatable(dfreg,
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      options = list(
        searchable = FALSE,
        dom = "t",
        paging = FALSE,
        columnDefs = list(list(className = "dt-center", targets = 0:2))
      )
    )
  })

  output$explain2 <- DT::renderDataTable({
    dfreg2 <- data.frame(
      state = c(
        "SA",
        "IBA",
        "Ent",
        "LT",
        "Ar",
        "SpD"
      ),
      full_name = c(
        "Summer Active",
        "InterBout Arousal",
        "Entrance",
        "Late Torpor ",
        "Arousing",
        "Spring Dark"
      ),
      body_temp = c(
        "37°C",
        "34.1 ± 2.9°C",
        "25.4 ± 1.8°C",
        "&nbsp;5.9 ± 0.5°C",
        "&nbsp;8.7 ± 2.1°C",
        "37°C"
      )
    )
    DT::datatable(dfreg2,
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      options = list(
        searchable = FALSE,
        dom = "t",
        paging = FALSE,
        columnDefs = list(list(className = "dt-center", targets = 0:1))
      )
    )
  })

  # graying out buttons, warnings/hints
  observe({
    if (input$tabMain == "plot") {
      enable("savePlot")
      disable("saveTable")
      output$savePlot <- savePlot
      if (rv$tabinit_plot == 0) {
        showNotification("tabs, modules, and columns of tables can be dragged and rearranged",
          type = "message"
        )
        rv$tabinit_plot <<- 1
      }
    } else if (input$tabMain == "table_data") {
      disable("savePlot")
      enable("saveTable")
      output$saveTable <- saveFiltered
      if (rv$tabinit_data == 0) {
        showNotification("type `low...high` to input custom range for numeric filtering on column",
          type = "message"
        )
        rv$tabinit_data <<- 1
      }
    } else if (input$tabMain == "table_AS") {
      disable("savePlot")
      enable("saveTable")
      output$saveTable <- saveFilteredAS
    } else if (input$tabMain == "line_plot") {
      enable("savePlot")
      disable("saveTable")
      output$savePlot <- savePlot5
    } else if (input$tabMain == "enrichment_plot") {
      enable("savePlot")
      enable("saveTable")
      output$savePlot <- savePlot2
      output$saveTable <- saveEnrich
      if (rv$tabinit_enrich == 0) {
        showNotification("add corresponding genes to cart by clicking on GO term bar",
          type = "message"
        )
        rv$tabinit_enrich <<- 1
      }
    } else if (input$tabMain == "heat_plot") {
      enable("savePlot")
      disable("saveTable")
      output$savePlot <- savePlot3
    } else if (input$tabMain == "kmer_analysis") {
      enable("savePlot")
      enable("saveTable")
      output$savePlot <- savePlot4
      output$saveTable <- saveK
      if (rv$tabinit_kmer == 0) {
        showNotification("Annotations: 5mer - Ray2013 + Encode, 6mer - Transite R, 7mer TargetScan mir seed",
          type = "message"
        )
        rv$tabinit_kmer <<- 1
      }
    } else if (input$tabMain == "venn") {
      enable("savePlot")
      disable("saveTable")
      output$savePlot <- savePlot6
      if (rv$tabinit_venn == 0) {
        showNotification("click on venn numbers to load associated genes into cart",
          type = "message"
        )
        rv$tabinit_venn <<- 1
      }
    } else if (input$tabMain == "about") {
      disable("savePlot")
      disable("saveTable")
    }
  })

  observeEvent(input$showpanel, {
    if (input$showpanel != TRUE) {
      removeCssClass("MAIN", "col-sm-12")
      addCssClass("MAIN", "col-sm-9")
      removeCssClass("SIDE", "col-sm-6")
      addCssClass("SIDE", "col-sm-3")
      shinyjs::show(id = "SIDE")
    }
    else {
      removeCssClass("MAIN", "col-sm-9")
      addCssClass("MAIN", "col-sm-12")
      removeCssClass("SIDE", "col-sm-3")
      addCssClass("SIDE", "col-sm-6")
      shinyjs::hide(id = "SIDE")
    }
  })

  onclick("bsgo", {
    rv$go <- 2
    removeModal()
  })

  onclick("bscancel", {
    rv$go <- 1
    removeModal()
  })

  observeEvent(input$tabMain == "heat_plot", {
    if (rv$go < 2) {
      rv$go <- rv$go - 1
    }
  })

  observeEvent(input$tabMain == "line_plot", {
    if (rv$go < 2) {
      rv$go <- rv$go - 1
    }
  })

  modalWarn <- draggableModalDialog(
    id = "bsconfirm",
    icon("exclamation-triangle"),
    p("plotting large number of genes may be\n slow and hard to interpret"),
    br(),
    footer = NULL,
    size = "s",
    easyClose = FALSE,
    fade = TRUE,
    actionButton("bsgo", "Go"),
    actionButton("bscancel", "Cancel")
  )

  observeEvent(input$tutorial, {
    introjs(session, options = list(
      "nextLabel" = "next",
      "prevLabel" = "prev",
      "skipLabel" = "skip",
      "overlayOpacity" = -1
    ))
  })
}

# Run the application
shinyApp(ui = ui, server = server, enableBookmarking = "url")
