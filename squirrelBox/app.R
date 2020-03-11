library(shiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(readr)
library(stringr)
library(shinyjs)
library(visNetwork)
library(valr)
library(igraph)
library(tibble)
library(plotly)
library(tidyr)
library(shinyBS)
library(feather)
library(crosstalk)
library(purrr)
library(ComplexHeatmap)
library(transite)
library(ggrepel)
library(shinythemes)
library(shinycssloaders)

options(stringsAsFactors = FALSE)
options(spinner.type = 6)
theme_set(theme_cowplot())
# options(shiny.reactlog = TRUE)

### general data settings
versionN <- 0.97
geoN <- "G1234"
set_shinytheme <- "paper"
track_name <- "hub_1519131_KG_HiC"
track_url <- "http://squirrelhub.s3-us-west-1.amazonaws.com/hub/hub.txt"
gmt_file <- "c5.bp.v7.0.symbols.gmt"
gmt_short <- "GO_"
sig_cut <- 0.001
table_cols <- c(
  "gene_id", # comment out to remove from table outputs
  "unique_gene_symbol",
  "gene_symbol",
  "clean_gene_symbol",
  "original_gene_name",
  "source"
)
orf_cols <- c(
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
  "majiq_directed" # ,
  # "orf"
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
if (file.exists("combined2.feather")) {
  combined2 <- read_feather("combined2.feather")
  combined3 <- read_feather("combined3.feather")
} else if (file.exists("combined2.csv")) {
  combined2 <- fread("combined2.csv", nThread = nt)
  combined3 <- fread("combined3.csv")
  combined <- combined3 %>% inner_join(combined2, by = "gene_id")
} else if (file.exists("combined.tsv")) {
  combined <- fread("combined.tsv")
} else {
  combined <- read_tsv("combined.tsv.gz",
    col_types = cols()
  )
}

# read annotation file to find ucsc track
bed <- read_tsv("final_tx_annotations_20200201.tsv.gz",
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
) %>%
  select(-contains("X")) %>%
  mutate(majiq_directed = factor(ifelse(is.na(majiq_directed), 0, 1)))

# read modules/clusters
mod <- read_feather("clusters.feather")
eigen <- read_tsv("cluster_patterns_matrices/reference_patterns.tsv") %>%
  rename(state = X1) %>%
  mutate(state = factor(state,
    levels = state_order
  ))

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
eigen_gg[["empty"]] <- ggplot()

# query function
comb_fil_factor <- function(combined2, combined3, inid) {
  combined3 <- combined3 %>%
    filter(gene_id %in% inid | unique_gene_symbol %in% inid)
  combined2 <- combined2 %>%
    filter(gene_id %in% inid | gene_id %in% (combined3$gene_id %>%
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

gmt <- gmt_to_list(gmt_file, rm = gmt_short)

if (file.exists(paste0(gmt_file, ".rds"))) {
  gmtlist <- readRDS(paste0(gmt_file, ".rds"))
} else {
  gmtlist <- gmt_to_list(gmt_file, rm = gmt_short, per = FALSE)
  gmtlist <- sapply(gmtlist, function(x) {
    intersect(x, str_to_upper(bed$clean_gene_symbol %>% unique()))
  })
  gmtlist <- gmtlist[sapply(gmtlist, length) >= 5]
  saveRDS(gmtlist, paste0(gmt_file, ".rds"))
}


domains <- read_csv("novel_domains.csv", col_types = "cc")

br_expr <- combined2 %>%
  filter(region %in% region_main) %>%
  pull(gene_id) %>%
  unique()

nonbr_expr <- combined2 %>%
  filter(region %in% region_main2) %>%
  pull(gene_id) %>%
  unique()

# load orf predictions
orfs <- read_feather("padj_orf.feather") %>%
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
  left_join(orfs %>% select(orf_cols, contains("LRT")), by = "gene_id") %>%
  left_join(mod, by = c("unique_gene_symbol" = "gene")) %>%
  mutate(source = factor(source)) %>%
  mutate_at(vars(contains("cluster")), factor) %>%
  distinct()

fulltbl_collapse <- fulltbl %>%
  group_by(gene_id) %>%
  arrange(desc(orf_len), .by_group = TRUE) %>%
  dplyr::slice(1)

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
  # call = letters[1:length(sig_cols)],
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
    arrange(state) %>%
    group_by(state) %>%
    summarize(letter = str_c(letter, collapse = ""))
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
    pval <- phyper(hitInSample - 1, hitInPop, failInPop, sampleSize, lower.tail = FALSE)
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
    mutate(go_len = lapply(gmtlist, length)) %>%
    arrange(desc(minuslog10), desc(len)) %>%
    select(pathway, pval, padj, minuslog10, pval, hits, len, go_len)
}

maj <- read_tsv("MAJIQ_dpsi_summary_sig_squirrelBox.tsv.gz") %>%
  mutate(
    region = factor(region),
    comp = factor(comp)
  ) %>%
  rename(comp_pair = "comp") %>%
  left_join(orfs %>% select(gene_id, contains("LRT")), by = "gene_id") %>%
  select(-gene_id) %>%
  distinct()

# seqs for kmer
seqs <- read_feather("utrs_sq.feather") %>% filter(gene_id %in% combined3$gene_id)

if (file.exists("seqs_precal.rds")) {
  print("load")
  seqs_precal <- readRDS("seqs_precal.rds")
} else {
  seqs_precal <- list()
  seqs_precal[["5mers_utr3"]] <- generateKmers(seqs %>% filter(str_length(utr3) >= 50) %>% pull(utr3),
    k = 5
  )
  seqs_precal[["7mers_utr3"]] <- generateKmers(seqs %>% filter(str_length(utr3) >= 50) %>% pull(utr3),
    k = 7
  )
  seqs_precal[["5mers_utr5"]] <- generateKmers(seqs %>% filter(str_length(utr5) >= 50) %>% pull(utr5),
    k = 5
  )
  seqs_precal[["7mers_utr5"]] <- generateKmers(seqs %>% filter(str_length(utr5) >= 50) %>% pull(utr5),
    k = 7
  )
  saveRDS(seqs_precal, "seqs_precal.rds")
}

comp_kmer <- function(df = seqs,
                      gene_vec,
                      col = "utr3",
                      bac = seqs_precal[["5mers_utr3"]],
                      k = 5,
                      cutoff = 50) {
  enq <- df %>%
    filter(str_to_upper(unique_gene_symbol) %in% gene_vec) %>%
    pull(col) %>%
    na.omit()
  if (length(enq) == 0) {
    return(NA)
  }
  enq <- enq[str_length(enq) >= cutoff]

  enq_res <- generateKmers(enq, k)

  res <- computeKmerEnrichment(enq_res,
    bac,
    permutation = FALSE,
    chisq.p.value.threshold = 0,
    p.adjust.method = "BH"
  )
  res$kmer <- str_replace_all(names(bac), "T", "U")
  res %>% arrange(adj.p.value)
}

fivemers <- read_csv("RBP_5mer.csv")
sevenmers <- read_csv("mir_7mer.csv")

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
$(function() {
    setTimeout(function(){
      var vals = [0];
      var powStart = 3;
      var powStop = 0;
      for (i = powStart; i >= powStop; i--) {
        var val = Math.pow(10, -i);
        val = parseFloat(val.toFixed(8));
        vals.push(val);
      }
      $("#pvalue").data("ionRangeSlider").update({"values":vals})
}, 5)});
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
  useShinyjs(),
  tags$head(tags$script(HTML(jscode))),
  titlePanel(div(
    class = "header", img(src = "logo.png", style = "width : 4%;"),
    "13-lined ground squirrel gene-level RNA-seq expression", style = "font-size:23px"
  )),
  fixedPanel(
    style = "z-index:100;",
    actionButton("back_to_top", label = "back_to_top"),
    right = 10,
    bottom = 10
  ),
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:23%;margin-top: 60px;z-index:10;",
      width = 3,
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
      div(style = "display: inline-block;vertical-align:top; width: 10px;", actionButton("Find", "Find")),
      bsTooltip("Find", "gene id/symbols accepted"),
      # br(.noWS="outside"),
      tabsetPanel(
        tabPanel(
          "options",
          br(.noWS = "outside"),
          div(id = "doPlotlydiv", checkboxInput("doPlotly", "interactive plots", value = F, width = NULL)),
          div(id = "doPadjdiv", checkboxInput("doPadj", "indicate sig", value = T, width = NULL)),
          div(id = "doNamediv", checkboxInput("doName", "additional labels", value = F, width = NULL)),
          checkboxInput("doBr", "plot brain data", value = T, width = NULL),
          checkboxInput("doTis", "plot non-brain data", value = F, width = NULL),
          div(id = "doEigendiv", checkboxInput("doEigen", "plot model clusters", value = F, width = NULL)),
          checkboxInput("doUcsc", "pull track", value = T, width = NULL),
          checkboxInput("doMod", "find module", value = T, width = NULL),
          checkboxInput("doKegg", "GO terms", value = T, width = NULL),
          div(id = "doNormdiv", checkboxInput("doNorm", "line plot norm to SA", value = F, width = NULL)),
          checkboxInput("doTooltips", "show hover tips", value = T, width = NULL),
          bsTooltip("doPlotlydiv", "toggles main and kmer plot"),
          bsTooltip("doPadjdiv", str_c("p<=", sig_cut)),
          bsTooltip("doNamediv", "samples in main, gene in line plot"),
          bsTooltip("doEigendiv", "assigned per region"),
          bsTooltip("doNormdiv", "otherwise centered by mean")
        ),
        tabPanel(
          "links",
          uiOutput("conn"),
          uiOutput("tab"), uiOutput("blastlink"),
          uiOutput("tab2"), uiOutput("tab3"), uiOutput("tab4"),
          downloadButton("savePlot", label = "save plot")
        ),
        tabPanel(
          "hide",
        )
      ),
      tabsetPanel(
        tabPanel(
          span("load", title = "load list of genes for analysis from file or interactive table"),
          fileInput("file", label = NULL),
          actionButton("Prev", "Prev"),
          actionButton("Next", "Next"),
          bsTooltip("Prev", "query previous gene on loaded list"),
          bsTooltip("Next", "query next gene on loaded list"),
          uiOutput("listn"),
          DT::dataTableOutput("tbllist"),
          style = "height:300px; overflow-y: scroll;"
        ),
        tabPanel(
          span("history", title = "history list of query genes"),
          DT::dataTableOutput("historyl"),
          style = "height:300px; overflow-y: scroll;"
        ),
        tabPanel(
          span("cart", title = "cart list of genes to save and export"),
          uiOutput("listn2"),
          actionButton("Add", "Add"),
          actionButton("Load", "Load"),
          downloadButton(
            outputId = "saveList",
            label = "save"
          ),
          bsTooltip("Add", "add current query gene to cart"),
          bsTooltip("Load", "sent to loaded list in side panel"),
          DT::dataTableOutput("tbllist2"),
          style = "height:300px; overflow-y: scroll;"
        )
      )
    ),
    mainPanel(
      width = 9,
      style = "z-index:1;margin-top: 60px;",
      tabsetPanel(
        id = "tabMain",
        tabPanel(
          title = span("main", 
                       title= "Plot expression box plot and other info of query gene"),
          value = "plot",
          uiOutput("boxPlotUI") %>% withSpinner(),
          uiOutput("EigenPlot") %>% withSpinner(),
          tableOutput("results"),
          bsCollapse(
            id = "tabs", multiple = TRUE, open = "called_orfs",
            bsCollapsePanel(tableOutput("orfinfo") %>% withSpinner(),
              title = "called_orfs",
              style = "primary"
            ),
            bsCollapsePanel(tableOutput("majinfo") %>% withSpinner(),
              title = "majiq_alternative_splicing",
              style = "warning"
            ),
            bsCollapsePanel(htmlOutput("ucscPlot") %>% withSpinner(),
              title = "UCSC browser plot",
              style = "success"
            ),
            bsCollapsePanel(tableOutput("gotab") %>% withSpinner(),
              title = "go_terms/domains",
              style = "info"
            )
          )
        ),
        tabPanel(
          title = span("transcript_gene", 
                       title= "Table of expression and other info of all genes/transcripts"),
          value = "table_orf",
          div(id = "doCollapsediv",
            style = "display: inline-block;width: 160px;",
            checkboxInput("doCollapse",
              "longest transcript",
              value = T,
              width = NULL
            )
          ),
          bsTooltip("doCollapsediv", "only show longest orf transcript for each gene"),
          downloadButton(
            outputId = "saveFiltered",
            label = "save filtered data"
          ),
          actionButton("loadtab", "load"),
          bsTooltip("loadtab", "sent to loaded list in side panel"),
          DT::dataTableOutput("tbl")
        ),
        tabPanel(
          title = span("majiq_alt",
                       title= "Table of majiq output for alternative splicing events"),
          value = "table_maj",
          downloadButton(
            outputId = "saveFiltered4",
            label = "save filtered data"
          ),
          DT::dataTableOutput("alt")
        ),
        tabPanel(
          title = span("line_plot",
                       title= "Plot expression of loaded gene list"),
          value = "line_plot",
          plotlyOutput("linePlot") %>% withSpinner()
        ),
        tabPanel(
          title = span("heat_plot",
                       title= "Plot Z-Score of loaded gene list as heat map"),
          value = "heat_plot",
          br(),
          fluidRow(
            column(width = 2,
              downloadButton("savePlot3",
                            label = "save plot"
              ),
            ),
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
          plotOutput("heatPlot") %>% withSpinner()
        ),
        tabPanel(
          title = span("GO_enrichment",
                       title= "GO term enrichment for loaded gene list"),
          value = "enrichment_plot",
          downloadButton("savePlot2",
            label = "save plot"
          ),
          downloadButton(
            outputId = "saveEnrich",
            label = "save table"
          ),
          plotlyOutput("richPlot") %>% withSpinner()
        ),
        tabPanel(
          title = span("kmer",
                       title= "kmer enrichment analysis and annotation for loaded gene list"),
          value = "kmer_analysis",
          downloadButton("savePlot4",
            label = "save plot"
          ),
          downloadButton(
            outputId = "saveK",
            label = "save table"
          ),
          tags$style(HTML(".radio-inline {margin-left: 5px;margin-right: 25px;}")),
          div(
            style = "display: inline-block;vertical-align:top;",
            radioButtons("utr", "UTR choice", c("5UTR", "3UTR"), selected = "3UTR", inline = TRUE)
          ),
          div(
            style = "display: inline-block;vertical-align:top;",
            radioButtons("km", "kmer length", c("5", "7"), selected = "5", inline = TRUE)
          ),
          div(
            style = "display: inline-block;vertical-align:top;",
            radioButtons("kmlab", "annotate kmer", c("yes", "no"), selected = "yes", inline = TRUE)
          ),
          selectInput("utrlen", NULL, choices = c(200, 500, 1000, "full length"), selected = "full length"),
          plotOutput("kmerPlot") %>% withSpinner()
        ),
        tabPanel(
          span("about", 
               title = "View version and author info"),
          value = "about",
          uiOutput("intro"),
          uiOutput("track"),
          uiOutput("rawdata"),
          uiOutput("GOversion"),
          uiOutput("version"),
          uiOutput("GitHub"),
          uiOutput("contact"),
          column(width =4, DT::dataTableOutput("explain"))
          
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
  rv$mod_df <- data.frame()

  # hide some checkboxes
  hide("doKegg")
  hide("doMod")
  hide("doUcsc")
  hide("utrlen")

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
    outputtab <- outputtab()
    inid <- outputtab$unique_gene_symbol
    plot_temp <- comb_fil_factor(combined2, combined3, inid)

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
        separate(comp, into = c("region", "state1", NA, "state2")) %>%
        select(-padj, -call) %>%
        mutate(call1 = as.numeric(call1)) %>%
        mutate(region = as.character(region))
      temp2$region <- region_order[factor(temp2$region, level = region_short) %>% as.numeric()]
      temp3 <- groups_to_letters_igraph(temp2) %>% mutate(region = factor(region, level = region_order))

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
          )
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
      plotOutput("boxPlot", width = 800, height = 600)
    } else {
      plotOutput("boxPlot", width = 800, height = 300)
    }
  })

  # boxplot-plotly
  boxPlotlyr <- reactive({
    g <- boxPlot1()
    output$boxPlot2 <- renderPlotly(ggplotly(g + facet_wrap(~region), tooltip = "text") %>%
      layout(hovermode = "closest"))
    if (input$doTis + input$doBr == 2) {
      plotlyOutput("boxPlot2", width = 800, height = 600)
    } else {
      plotlyOutput("boxPlot2", width = 800, height = 300)
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
      plotOutput("boxPlot3", width = 800, height = 300)
    }
  })

  # boxplot - models
  output$boxPlot3 <- renderPlot({
    if (input$doEigen != T) {
      g <- ""
      g
    } else {
      outputtab <- outputtab()
      inid <- outputtab$unique_gene_symbol
      if (nrow(rv$mod_df) == 0) {
        mods <- c("empty", "empty", "empty")
      }
      mods <- (rv$mod_df[1, ] %>% unlist())[-1]
      cowplot::plot_grid(
        plotlist = map(mods, function(x) eigen_gg[[x]]),
        ncol = 3
      )
    }
  })

  # finding module/cluster info
  output$conn <- renderUI({
    if (input$doMod != T) {
      return()
    }

    if (nrow(rv$mod_df) == 0) {
      mod1 <- "low expression everywhere"
    } else {
      mod1 <- rv$mod_df %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        mutate(text = str_c(rowname, V1, sep = ":")) %>%
        pull(text) %>%
        str_c(collapse = "<br>")
    }

    HTML(str_c(
      mod1 # , r
    ))
  })

  # filter data
  outputtab <- reactive({
    inid <- inid()
    if ((inid %in% orfs$gene_id) | (inid %in% orfs$unique_gene_symbol)) {
      temp_orfs <- orfs %>% filter((gene_id == inid) | (unique_gene_symbol == inid))
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

    filtered <- comb_fil_factor(combined2, combined3, inid) %>%
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
  output$results <- renderTable(
    {
      outputtab()
    },
    digits = 0
  )

  # orf call table
  output$orfinfo <- renderTable(
    {
      outputtab()
      if (nrow(rv$temp_orfs) == 0) {
        return(rv$temp_orfs)
      }
      rv$temp_orfs %>% select(orf_cols)
    },
    digits = 0
  )

  # majik report table
  output$majinfo <- renderTable(
    {
      temp <- maj %>% filter(unique_gene_symbol == outputtab()$unique_gene_symbol[1])
      if (nrow(temp) == 0) {
        temp <- data.frame()
      }
      temp
    },
    digits = 0
  )

  # goterm table
  output$gotab <- renderTable({
    if (input$doKegg != T) {
      return()
    }
    outputtab <- outputtab()
    temp1 <- gmt %>% filter(genes == str_to_upper(outputtab$clean_gene_symbol))
    if (nrow(temp1) == 0) {
      temp1 <- domains %>% filter(gene_id %in% outputtab$gene_id)
    }
    temp1
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
    url <- a(outputtab$unique_gene_symbol,
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
    tagList("trackhub:", url)
  })

  # link genbank
  output$tab2 <- renderUI({
    outputtab <- outputtab()
    if (nrow(outputtab) > 1) {
      outputtab <- outputtab[1, ]
    }
    clean <- a(outputtab$unique_gene_symbol,
      href = str_c(
        "https://www.ncbi.nlm.nih.gov/gene/?term=",
        str_remove(outputtab$unique_gene_symbol, "_.+"),
        "[sym]+AND+human[ORGN]"
      )
    )
    tagList("genbank:", clean)
  })

  # link hgnc
  output$tab3 <- renderUI({
    outputtab <- outputtab()
    if (nrow(outputtab) > 1) {
      outputtab <- outputtab[1, ]
    }
    clean <- a(outputtab$unique_gene_symbol,
      href = str_c(
        "https://www.genenames.org/data/gene-symbol-report/#!/symbol/",
        str_remove(outputtab$unique_gene_symbol, "_.+")
      )
    )
    tagList("hgnc:", clean)
  })

  # link genecard
  output$tab4 <- renderUI({
    outputtab <- outputtab()
    if (nrow(outputtab) > 1) {
      outputtab <- outputtab[1, ]
    }
    clean <- a(outputtab$unique_gene_symbol,
      href = str_c(
        "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
        str_remove(outputtab$unique_gene_symbol, "_.+")
      )
    )
    tagList("genecard:", clean)
  })

  # link blast
  output$blastlink <- renderUI({
    if (rv$blast != "" & !(is.na(rv$blast))) {
      outputtab <- outputtab()
      if (nrow(outputtab) > 1) {
        outputtab <- outputtab[1, ]
      }
      orf <- rv$blast
      url <- a(outputtab$unique_gene_symbol,
        href = str_c(
          "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastp&DATABASE=nr&QUERY=",
          orf
        )
      )
      tagList("blast orf:", url)
    } else {
      return()
    }
  })

  # save plot as pdf
  output$savePlot <- downloadHandler(
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
        w <- 8
        h <- 6
      } else {
        w <- 8
        h <- 3
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
      summarize(counts = mean(2^log2_counts))
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
        # l$sample <- "0"
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
  output$linePlot <- renderPlotly({
    set.seed(1)
    linetemp()
    if (length(historytablist) == 0) {
      return(ggplotly(ggplot()))
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

    fac <- input$doTis + input$doBr
    if (input$doName == T) {

      # highlight(ggplotly(g, tooltip = "text"),"plotly_hover")
      ggplotly(g, tooltip = "text", height = 300 * fac, width = 800) %>% layout(
        autosize = FALSE,
        showlegend = TRUE
      )
    } else {
      # highlight(ggplotly(g, tooltip = "text"),"plotly_hover")
      ggplotly(g, tooltip = "text", height = 300 * fac, width = 800)
    }
  })

  # heatmap
  heatPlot1 <- reactive({
    set.seed(1)
    temp <- linetemp()
    if (length(historytablist) == 0) {
      return(NA)
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

  output$heatPlot <- renderPlot(heatPlot1())

  output$savePlot3 <- downloadHandler(
    filename = "heatplot.pdf",
    content = function(file) {
      h <- heatPlot1()
      if (input$doAutoresize) {
        pdf(file, width = (h@matrix %>% dim())[2], height = (h@matrix %>% dim())[1])
      } else {
        pdf(file, width = 8, height = 8)
      }
      print(h)
      dev.off()
    }
  )

  richPlot1 <- reactive({
    rv$line_refresh
    set.seed(1)
    if (length(historytablist) == 0) {
      return(ggplotly(ggplot()))
    }
    genevec <- unique_to_clean(historytablist, namedvec) %>% str_to_upper()
    tops <- fisher(genevec, gmtlist, length_detected_genes)
    tops <- tops %>% dplyr::slice(1:max(min(which(tops$padj > 0.01)), 15))
    tops <<- tops
    ggplot(tops %>% dplyr::slice(1:15), aes(x = reorder(pathway, minuslog10), y = minuslog10, fill = -minuslog10, text = len)) +
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
      scale_y_continuous(expand = c(0, 0))
  })

  output$richPlot <- renderPlotly({
    ggplotly(richPlot1(), source = "richPlot", tooltip = "text", height = 600, width = 800) %>%
      layout(autosize = F) %>%
      highlight()
  })

  output$savePlot2 <- downloadHandler(
    filename = "enriched.pdf",
    content = function(file) {
      ggplot2::ggsave(file, plot = richPlot1(), device = "pdf", width = 8, height = 6)
    }
  )

  observeEvent(event_data("plotly_click", source = "richPlot"), {
    rv$richsel <- event_data("plotly_click", source = "richPlot")
    carttablist <<- tops[16 - rv$richsel$y, ] %>%
      pull(hits) %>%
      str_split(",") %>%
      unlist()
    rv$listn2 <- length(carttablist)
  })

  output$saveEnrich <- downloadHandler("enrich_list.csv", content = function(file) {
    write_csv(tops, file)
  })

  # kmer analysis

  kmertemp <- reactive({
    rv$line_refresh
    set.seed(1)
    if (length(historytablist) == 0) {
      # return(ggplotly(ggplot()))
      return(data.frame())
    }
    genevec <- unique_to_clean(historytablist, namedvec) %>% str_to_upper()
    if (input$utrlen == "full length") {
      lenchoice <- 0
    } else {
      lenchoice <- input$utrlen
    }
    if (input$utr == "3UTR") {
      utrchoice <- "utr3"
    } else {
      utrchoice <- "utr5"
      lenchoice <- -lenchoice
    }
    precal <- paste0(input$km, "mers_", utrchoice)
    topsk <- comp_kmer(gene_vec = genevec,
                       bac = seqs_precal[[precal]], 
                       col = utrchoice, 
                       k = as.numeric(input$km))
    topsk <- topsk %>%
      mutate(minuslog10 = -log10(adj.p.value), enrichment = log2(enrichment))
    if (input$km == "5") {
      # topsk <- topsk %>%
      #   dplyr::slice(1:max(min(which(topsk$adj.p.value > 0.05)), 15)) %>%
      #   mutate(minuslog10 = -log10(adj.p.value), enrichment = log2(enrichment)) %>%
      #   left_join(fivemers, by = c("kmer" = "fivemer"))
      topsk <- topsk %>% left_join(fivemers, by = c("kmer" = "fivemer"))
    } else {
      # topsk <- topsk %>%
      #   dplyr::slice(1:max(min(which(topsk$adj.p.value > 0.05)), 15)) %>%
      #   mutate(minuslog10 = -log10(adj.p.value), enrichment = log2(enrichment)) %>%
      #   left_join(sevenmers, by = c("kmer" = "sevenmer")) %>%
      #   rename(RBP = "mir")
      topsk <- topsk %>%
        left_join(sevenmers, by = c("kmer" = "sevenmer")) %>%
        rename(RBP = "mir")
    }
    topsk
  })
  
  kmerPlot1 <- reactive({
    topsk <- kmertemp()
    if (nrow(topsk) == 0) {
      return(ggplot())
    }
    if (input$kmlab == "yes") {
      topsk <- topsk %>% mutate(sig = factor(ifelse(minuslog10 >= 2, "sig", "insig"))) %>%
        mutate(text2 = ifelse((sig =="sig" & row_number() <= 15), str_c(kmer, RBP, sep = "\n"), ""))
    } else {
      topsk <- topsk %>% mutate(sig = factor(ifelse(minuslog10 >= 2, "sig", "insig"))) %>%
        mutate(text2 = ifelse((sig =="sig" & row_number() <= 15), kmer, ""))
    }
    # ggplot(topsk %>% dplyr::slice(1:15), aes(x = reorder(kmer, minuslog10), y = minuslog10, text = RBP)) +
    #   geom_bar(stat = "identity", aes(fill = enrichment)) +
    #   scale_fill_gradient2(
    #     low = "blue",
    #     mid = "white",
    #     high = "red",
    #     midpoint = 0
    #   ) +
    #   coord_flip() +
    #   xlab("kmer") +
    #   labs(fill = "log2(enrichment)") +
    #   cowplot::theme_minimal_vgrid() +
    #   theme(
    #     axis.text.y = element_text(size = 6),
    #     axis.text.x = element_text(size = 10),
    #     axis.title.y = element_text(size = 10),
    #     axis.title.x = element_text(size = 10),
    #     legend.position = "right",
    #     legend.title = element_text(size = 6),
    #     legend.text = element_text(size = 6)
    #   ) +
    #   scale_y_continuous(expand = c(0, 0))
    ggplot(topsk, aes(x = enrichment,
                      y = minuslog10, 
                      text = str_c(kmer, RBP, sep = "\n"),
                      label = text2)) +
      geom_point(aes(color = sig)) +
      scale_color_manual(values = c("gray", "red")) +
      ggrepel::geom_text_repel(box.padding = 0.05, size = 3) +
      xlab("log2enrichment") +
      labs(color = "")
  })

  output$kmerPlot <- renderPlot({
    kmerPlot1()
    # ggplotly(kmerPlot1(), source = "kmerPlot", tooltip = "text", height = 600, width = 800) %>%
    #   layout(autosize = F)
  })
  
  output$savePlot4 <- downloadHandler(
    filename = "kmers.pdf",
    content = function(file) {
      ggplot2::ggsave(file, plot = kmerPlot1(), device = "pdf", width = 8, height = 6)
    }
  )

  output$saveK <- downloadHandler("enrich_kmer.csv", content = function(file) {
    write_csv(kmertemp(), file)
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
    historytablist <<- autocomplete_list[str_to_upper(autocomplete_list) %in% str_to_upper(v_genes %>%
      unique())]
    rv$line_refresh <- rv$line_refresh + 1
  })

  # cart list
  onclick("Add", {
    carttablist <- unique(c(historytab[1], carttablist))
    rv$listn2 <- length(carttablist)
  })

  output$listn2 <- renderUI({
    HTML(str_c("# in cart:", rv$listn2))
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

  onclick("Prev", {
    rv$listn <- rv$listn - 1
    if (rv$listn < 1) {
      rv$listn <- 1
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

  onclick("Next", {
    rv$listn <- rv$listn + 1
    if (rv$listn > length(historytablist)) {
      rv$listn <- length(historytablist)
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
    HTML(paste0(rv$listn, " of ", length(historytablist)))
  })

  # list loaded genes as table
  output$tbllist <- DT::renderDataTable({
    rv$line_refresh
    if (length(historytablist) > 0) {
      DT::datatable(data.table::as.data.table(list(historytablist)),
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
    DT::datatable(
      orftbl() %>%
        select(
          unique_gene_symbol,
          contains("LRT"),
          contains("cluster"),
          everything()
        ),
      filter = "top",
      escape = FALSE,
      selection = "single",
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(targets = c(7), visible = TRUE, width = "150px"),
          list(targets = c(8), visible = TRUE, width = "150px"),
          list(targets = c(9), visible = TRUE, width = "150px")
        ),
        scrollX = TRUE,
        autoWidth = TRUE
      )
    ) %>%
      DT::formatRound(columns = c(2, 3, 4, 5, 6, 7), digits = 4)
  })

  output$saveFiltered <- downloadHandler("filtré.csv", content = function(file) {
    s <- input$tbl_rows_all
    write_csv((orftbl() %>%
      select(
        unique_gene_symbol,
        contains("LRT"),
        contains("cluster"),
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
  output$alt <- DT::renderDataTable({
    DT::datatable(
      maj %>%
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
      options = list(
      scrollX = TRUE,
      autoWidth = TRUE)
    )
  })

  output$saveFiltered4 <- downloadHandler("filtré.csv", content = function(file) {
    s <- input$alt_rows_all
    write_csv((maj %>%
      select(
        unique_gene_symbol,
        contains("significant"),
        LSV_ID, A5SS, A3SS, ES,
        everything()
      ))[s, ], file)
  })

  observeEvent(input$alt_rows_selected, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = maj[input$alt_rows_selected, "unique_gene_symbol"],
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
      track_url)
    clean <- a("UCSC browser",
               href = url
    )
    tagList(tags$h6("Newly assembled genome and annotated transcriptome are hosted on ", clean))
  })
  
  output$rawdata <- renderUI({
    url <- str_c(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geoN)
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
    dfreg <- data.frame(region = region_order,
                     short = region_short,
                     letter = region_one)
    DT::datatable(dfreg,
                  escape = FALSE,
                  selection = "none",
                  rownames = FALSE,
                  options = list(searchable = FALSE, 
                                 dom = "t", 
                                 paging = FALSE,
                                 columnDefs = list(list(className = 'dt-center', targets = 0:2)))
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server, enableBookmarking = "url")
