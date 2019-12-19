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
options(stringsAsFactors = FALSE)
theme_set(theme_cowplot())
# options(shiny.reactlog = TRUE)


### general data settings
use_folder <- "wgcna11" # change to get to old version of data
track_name <- "hub_1519131_KG_HiC"
track_url <- "http://squirrelhub.s3-us-west-1.amazonaws.com/hub/hub.txt"
gmt_file <- "c5.all.v6.2.symbols.gmt"

### sample settings, define state colors and order, region order
### possibly read info out of csv file
state_cols <- c(
  SA = rgb(255, 0, 0, maxColorValue = 255),
  IBA = rgb(67, 205, 128, maxColorValue = 255),
  Ent = rgb(155, 48, 255, maxColorValue = 255),
  LT = rgb(25, 25, 112, maxColorValue = 255),
  Ar = rgb(0, 0, 255, maxColorValue = 255),
  SpD = rgb(255, 165, 0, maxColorValue = 255)
)
state_order <- c(
  "SA",
  "IBA",
  "Ent",
  "LT",
  "Ar",
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
region_short <- c(
  "fore",
  "hy",
  "med"
)
region_letter <- c(
  "B",
  "H",
  "M"
)


# read database
if (file.exists("combined2.feather")) {
  combined2 <- read_feather("combined2.feather")
  combined3 <- read_feather("combined3.feather")
} else if (file.exists("combined2.csv")) {
  # combined2 <- read_csv("combined2.csv", col_types = "ccncc")
  # combined3 <- read_csv("combined3.csv", col_types = "cccccc")
  combined2 <- fread("combined2.csv", nThread = nt)
  combined3 <- fread("combined3.csv")
  combined <- combined3 %>% inner_join(combined2, by = "gene_id")
} else if (file.exists("combined.tsv")) {
  # combined <- read_tsv("combined.tsv", col_types = cols())
  combined <- fread("combined.tsv")
} else {
  combined <- read_tsv("combined.tsv.gz",
    col_types = cols()
  )
}

comb_fil_factor <- function(combined2, combined3, inid) {
  combined3 <- combined3 %>% filter(gene_id == inid | unique_gene_symbol == inid)
  combined2 <- combined2 %>% filter(gene_id == inid | gene_id == combined3$gene_id[1]) %>%
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
# combined <- combined %>% mutate(
#   state = factor(state,
#     levels = state_order
#   ),
#   region = factor(region,
#     levels = region_order
#   )
# )

# lists
autocomplete_list <- str_sort(c(
  combined3$unique_gene_symbol,
  combined3$gene_id
) %>% unique(),
decreasing = T
)

# read annotation file to find ucsc track
bed <- read_tsv("final_annot_20191112.bed12",
  col_names = F,
  col_types = "cnncncnnnnncccc"
) %>%
  select(1:6, 13)

colnames(bed) <- c(
  "chrom",
  "start",
  "end",
  "unique_gene_symbol",
  "score",
  "strand",
  "gene_id"
)

# empty history list to start
historytab <- c()

# read modules
for (reg in region_short) {
  eval(parse(text = paste0(reg, '_modules <- suppressWarnings(read_csv(paste0(use_folder, \"/', reg, '_modules.csv\"), 
                                          col_types = \"ncncn\"))')))
}

# read node igraph object
for (reg in region_short) {
  eval(parse(text = paste0(reg, '_ig <- suppressWarnings(readRDS(paste0(use_folder, \"/', reg, '_conn_list\")))')))
}

for (reg in region_short) {
  eval(parse(text = paste0(reg, "_net <- suppressWarnings(toVisNetworkData(", reg, "_ig))")))
}

igraph_to_df <- function(ig, mod) {
  temp <- as_edgelist(ig)
  temp <- c(temp[, 1], temp[, 2]) %>%
    table() %>%
    as.tibble()
  colnames(temp) <- c("gene", "n")
  temp <- temp %>%
    left_join(mod, by = "gene") %>%
    select(-X1) %>%
    group_by(module_n) %>%
    mutate(rank = rank(-n)) %>%
    ungroup()
  temp
}

for (reg in region_short) {
  eval(parse(text = paste0(reg, "_temp4 <- suppressWarnings(igraph_to_df(", reg, "_ig, ", reg, "_modules))")))
}

# eigengene plots
for (reg in region_short) {
  eval(parse(text = paste0(reg, '_gg <- suppressWarnings(readRDS(paste0(use_folder, \"/', reg, '_gg\")))')))
}

# read go terms and TFs
gmt_to_list <- function(path,
                        cutoff = 0,
                        sep = "\thttp://www.broadinstitute.org/gsea/msigdb/cards/.*?\t",
                        rm = "REACTOME_") {
  df <- readr::read_csv(path,
    col_names = F, col_types = cols()
  )
  df <- tidyr::separate(df,
    X1,
    sep = sep,
    into = c("path", "genes")
  )
  df <- dplyr::mutate(df,
    path = stringr::str_replace(
      df$path,
      rm,
      ""
    ),
    genes = stringr::str_split(
      genes,
      "\t"
    )
  )
  tidyr::unnest_legacy(df, genes)
}

gmt <- gmt_to_list(gmt_file,
  rm = "^GO_"
)

refs <- combined3 %>% distinct(
  unique_gene_symbol,
  clean_gene_symbol
)
TFs <- gmt %>%
  filter(
    str_detect(
      path,
      str_to_upper("transcription_repressor_activity|transcription_factor_activity")
    ),
    path != "VIRAL_TRANSCRIPTION"
  ) %>%
  pull(genes) %>%
  unique()
refTFs <- refs %>%
  mutate(clean_gene_symbol = str_to_upper(clean_gene_symbol)) %>%
  filter(clean_gene_symbol %in% str_to_upper(TFs)) %>%
  pull(unique_gene_symbol)

# load orf predictions
orfs <- read_csv("padj_orf.csv") %>%
  select(gene_id,
    orf_len = len,
    exons,
    rna_len = transcript,
    orf,
    unique_gene_symbol,
    everything()
  )

starorfs <- orfs %>%
  group_by(unique_gene_symbol) %>%
  arrange(desc(orf)) %>%
  dplyr::slice(1) %>%
  filter(
    exons > 1,
    str_detect(unique_gene_symbol, "^G[0-9]+|_"),
    orf_len >= 100
  )
domains <- read_csv("novel_domains.csv", col_types = "cc")

# padj functions
find_padj <- function(region, state, tbl) {
  temp <- str_c(tbl[str_sub(tbl, 1, 1) == str_to_lower(str_sub(region, 1, 1)) &
    str_detect(tbl, state)],
  collapse = "\n"
  )
  if (length(temp) == 0) {
    temp <- "NA"
  }
  temp
}

sig_cols <- colnames(orfs)[colnames(orfs) %>% str_detect("vs")]
sig_sym <- data.frame(
  call = letters[1:length(sig_cols)],
  row.names = sig_cols
)

calls_sig <- function(padj, sig_sym) {
  temp <- cbind(padj, sig_sym)
  temp <- temp %>%
    rownames_to_column("comp") %>%
    mutate(call1 = ifelse(padj <= 0.05, 1, 0))
  temp
}

find_groups_igraph <- function(df) {
  df <- df %>% select(-region)
  g <- graph_from_data_frame(df, directed = FALSE)
  cg <- max_cliques(g)
  lapply(cg, names)
}

sort_groups <- function(groups) {
  all_groups <- state_order
  leftout <- list(setdiff(state_order, unlist(groups)))
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
    g3 <- sort_groups(g2)
    g3$region <- x
    return(g3)
  })
  do.call(rbind, g)
}

# read in rf importance
imp <- read_csv("rf_imp.csv")

# read in rf selection results
rf.vs1 <- readRDS("rfvs.rds")
trf <- rf.vs1$selec.history %>% 
  as_tibble() %>% 
  mutate(Vars.in.Forest = str_split(Vars.in.Forest, " \\+ "))
rfvars<- trf$Vars.in.Forest %>%
  unlist() %>%
  table() %>%
  sort(decreasing = TRUE) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>%
  mutate(rank = rank(-Freq, ties.method = "min")) %>%
  select(-2)
colnames(rfvars) <- c("gene", "rank")

vars_set <- function(rfvars, n) {
  rfvars %>% filter(rank <= n) %>% pull(gene)
}

# slightly faster than plot_grid?
gg_to_facet <- function(region_short, gg_list) {
  df_title <- lapply(gg_list, function(g) {
    g$labels$title
  })
  names(df_title) <- region_short
  
  df_list <- mapply(function(g, reg) {
    temp <- g$data
    if (is.null(nrow(temp))) {
      return(data.frame(eigengenes = NA, population = state_order[1], region = reg))
    }
    temp$region <- reg
    temp
  }, gg_list, region_short, SIMPLIFY = FALSE)
  
  temp2 <- do.call(rbind, df_list) %>% mutate(
    population = factor(population,
                   levels = state_order
    ),
    region = factor(region,
                    levels = region_short
    ))
  
  ggplot(temp2, aes(y = eigengenes, x = population)) + 
    geom_boxplot(aes(fill = population), outlier.shape = NA) +
    scale_fill_manual(values = state_cols) +
    facet_wrap(~region, scales = "free", labeller = as_labeller(unlist(df_title))) + 
    scale_x_discrete(limits = state_order) +
    # geom_point(aes(color = population), position = position_jitter(seed = 1)) + 
    # scale_color_manual(values = state_cols) + 
    theme(legend.position = "none")
}

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
  tags$style("
      .checkbox { /* checkbox is a div class*/
        line-height: 20px;
        margin-top: -15px;
        margin-bottom: -15px; /*set the margin, so boxes don't overlap*/
      }
  "),
  useShinyjs(),
  tags$head(tags$script(HTML(jscode))),
  titlePanel(div(img(src = "logo.png", style ="width : 4%; display: inline-block;"),
                 "13-lined ground squirrel gene-level RNA-seq expression")),
  fixedPanel(
    style="z-index:100;",
    actionButton("back_to_top", label = "back_to_top"),
    right = 10,
    bottom = 10
  ),
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:inherit;",
      width = 3,
      div(
        style = "display: inline-block;vertical-align:top; width: 200px;",
        tagAppendAttributes(selectizeInput("geneID",
          label = NULL,
          selected = "",
          choices = ""
        ),
        `data-proxy-click` = "Find"
        )
      ),
      div(style = "display: inline-block;vertical-align:top; width: 10px;", actionButton("Find", "Find")),
      tags$hr(style = "border-color: green;"),
      tabsetPanel(
        tabPanel(
          "options",
          br(),
          checkboxInput("doPlotly", "interactive padj", value = F, width = NULL),
          checkboxInput("doPadj", "indicate sig", value = F, width = NULL),
          checkboxInput("doName", "label by sample", value = F, width = NULL),
          checkboxInput("doTis", "plot non-brain", value = F, width = NULL),
          checkboxInput("doEigen", "plot eigengenes", value = T, width = NULL),
          checkboxInput("doUcsc", "download track", value = F, width = NULL),
          checkboxInput("doMod", "find module", value = T, width = NULL),
          checkboxInput("doNet", "plot network", value = F, width = NULL),
          checkboxInput("doKegg", "GO terms", value = T, width = NULL)
        ),
        tabPanel(
          "links",
          selectInput("region", label = NULL, choices = as.list(region_short), selected = region_short[1]),
          uiOutput("conn"),
          tags$hr(style = "border-color: green;"),
          uiOutput("tab"), uiOutput("blastlink"),
          uiOutput("tab2"), uiOutput("tab3"), uiOutput("tab4"),
          downloadButton("savePlot", label = "download plot")
        )
      ),
      # br(),
      tags$hr(style = "border-color: green;"),
      tabsetPanel(
        tabPanel(
          "file",
          fileInput("file", label = NULL),
          actionButton("Prev", "Prev"),
          actionButton("Next", "Next"),
          uiOutput("listn")
        ),
        tabPanel(
          "history",
          uiOutput("history1"),
          uiOutput("history2"),
          uiOutput("history3"),
          uiOutput("history4"),
          uiOutput("history5"),
          uiOutput("history6"),
          uiOutput("history7"),
          uiOutput("history8"),
          uiOutput("history9"),
          uiOutput("history10")
        ),
        tabPanel(
          "cart",
          uiOutput("listn2"),
          actionButton("Add", "Add"),
          downloadButton(
            outputId = "saveList",
            label = "cart to TXT"
          )
        )
      )
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabMain",
        tabPanel(
          title = "plot",
          value = "plot",
          uiOutput("boxPlotUI"),
          bsModal("modalPDF",
            title = "module-trait",
            trigger = "conn",
            size = "large",
            htmlOutput("pdfview")
          ),
          # plotOutput("boxPlot", width = 800, height = 600),
          uiOutput("EigenPlot"),
          tags$hr(style = "border-color: green;"),
          tableOutput("results"),
          tableOutput("orfinfo"),
          tags$hr(style = "border-color: green;"),
          htmlOutput("ucscPlot"),
          tags$hr(style = "border-color: green;"),
          visNetworkOutput("connPlot"),
          tags$hr(style = "border-color: green;"),
          tableOutput("gotab")
        ),
        tabPanel(
          title ="table_orf",
          value = "table_orf",
          downloadButton(
            outputId = "saveFiltered",
            label = "download filtered data"
          ),
          DT::dataTableOutput("tbl")
        ),
        tabPanel(
          title ="table_RF",
          value = "table_RF",
          downloadButton(
            outputId = "saveFiltered2",
            label = "download filtered data"
          ),
          DT::dataTableOutput("tbl2")
        ),
        tabPanel(
          title ="table_varsel",
          value = "table_varsel",
          div(plotlyOutput("oob", width = 400, height = 300, inline = TRUE),
          plotlyOutput("oob2", width = 400, height = 300, inline = TRUE),),
          div(
          downloadButton(
            outputId = "saveFiltered3",
            label = "download gene list"
          ),
          uiOutput("sel", inline = TRUE)),
          DT::dataTableOutput("tbl3")
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
  rv$old <- ""
  rv$blast <- ""
  rv$pval <- data.frame()
  rv$act_modules <- eval(parse(text = paste0(region_short[1], "_modules")))
  rv$act_temp4 <- eval(parse(text = paste0(region_short[1], "_temp4")))
  rv$act_net <- eval(parse(text = paste0(region_short[1], "_net")))
  rv$temp_orfs <- data.frame()
  rv$listn <- 1
  rv$listn2 <- 0
  rv$cart <- 0
  rv$xsel <- "NA"

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

  observeEvent(input$region, {
    rv$act_modules <- eval(parse(text = paste0(input$region, "_modules")))
    rv$act_temp4 <- eval(parse(text = paste0(input$region, "_temp4")))
    rv$act_net <- eval(parse(text = paste0(input$region, "_net")))
  })

  observeEvent(input$Find, {
    updateTabsetPanel(session, 
                      "tabMain",
                      selected = "plot")
  })
  
  inid <- eventReactive(input$Find,
    {
      rv$old <<- input$geneID
      shinyjs::runjs("window.scrollTo(0, 0)")
      if (rv$init == 1) {
        rv$init <<- 2
      } else {
        updateQueryString(paste0("?gene=", input$geneID), mode = "push")
      }
      input$geneID
    },
    ignoreNULL = T
  )

  historytab <- c()
  historytablist <- c()
  carttablist <- c()

  boxPlot1 <- reactive({
    outputtab <- outputtab()
    inid <- outputtab$unique_gene_symbol
    plot_temp <- comb_fil_factor(combined2, combined3, inid)
    mis <- setdiff(region_main, plot_temp$region %>% unique() %>% as.character())
    if (length(mis) > 0) {
      for (element in mis) {
        l <- as.list(plot_temp[1, 1:6])
        l$region <- element
        l$sample <- "0"
        l$log2_counts <- NA
        l$state <- state_order[1]
        plot_temp <- rbind(plot_temp, l)
      }
    }
    if (input$doTis != T) {
      plot_temp <- plot_temp %>% filter(region %in% region_main)
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
      # geom_point(aes(color = sample), position = position_jitter(seed = 1))
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
      temp2 <- temp2
      temp3 <- groups_to_letters_igraph(temp2)

      agg <- aggregate(log2_counts ~ state + region, plot_temp, max)
      agg_min <- aggregate(log2_counts ~ state + region, plot_temp, min)
      agg$min <- agg_min$log2_counts
      agg2 <- agg %>%
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

    g
  })

  boxPlotr <- reactive({
    g <- boxPlot1()
    output$boxPlot <- renderPlot(g)
    if (input$doTis == T) {
      plotOutput("boxPlot", width = 800, height = 600)
    } else {
      plotOutput("boxPlot", width = 800, height = 300)
    }
  })

  boxPlotlyr <- reactive({
    g <- boxPlot1()
    output$boxPlot2 <- renderPlotly(ggplotly(g + facet_wrap(~region), tooltip = "text") %>%
      layout(hovermode = "closest"))
    if (input$doTis == T) {
      plotlyOutput("boxPlot2", width = 800, height = 600)
    } else {
      plotlyOutput("boxPlot2", width = 800, height = 300)
    }
  })

  output$boxPlotUI <- renderUI({
    if (input$doPlotly == FALSE) {
      boxPlotr()
    } else {
      boxPlotlyr()
    }
  })

  output$EigenPlot <- renderUI({
    if (input$doEigen != T) {
      plotOutput("boxPlot3", height = 1)
    } else {
      plotOutput("boxPlot3", width = 800, height = 300)
    }
  })

  output$boxPlot3 <- renderPlot({
    if (input$doEigen != T) {
      g <- ""
      g
    } else {
      outputtab <- outputtab()
      inid <- outputtab$unique_gene_symbol
      for (reg in region_short) {
        eval(parse(text = paste0(reg, "_mod <- ", reg, "_modules %>% filter(gene == inid) %>% pull(module_n)")))
      }

      for (reg in region_short) {
        eval(parse(text = paste0(reg, "_fig <- tryCatch({", reg, "_gg[[as.numeric(", reg, "_mod) + 1]]}, 
                                          error = function(err) {
                                            return(ggplot() + theme_void())})")))
      }

      # eval(parse(text = paste0("cowplot::plot_grid(", str_c(region_short, "_fig", collapse = ","), ", ncol =", length(region_short), ")")))
      eval(parse(text = paste0("gg_to_facet(region_short, list(", str_c(region_short, "_fig", collapse = ","), "))")))
    }
  })

  output$connPlot <- renderVisNetwork({
    if (input$doMod != T | input$doNet != T) {
      return()
    }
    outputtab <- outputtab()
    queryid <- outputtab$unique_gene_symbol
    edgeq <- rv$act_net$edges %>% filter(from == queryid | to == queryid)
    nodeq <- c(edgeq$from, edgeq$to) %>% unique()
    edgeq2 <- tryCatch(
      {
        rv$act_net$edges %>%
          filter(from %in% nodeq | to %in% nodeq) %>%
          mutate(color = "gray", opacity = 0, width = 0)
      },
      error = function(err) {
        return(data.frame())
      }
    )
    nodeq2 <- tryCatch(
      {
        rv$act_net$nodes[nodeq, ] %>%
          left_join(table(c(edgeq2$from, edgeq2$to)) %>%
            as.data.frame(stringsAsFactors = F), by = c("id" = "Var1")) %>%
          mutate(
            value = Freq,
            color = ifelse(id == queryid, "red", "lightblue"),
            shape = ifelse(id %in% refTFs,
              "square",
              ifelse(id %in% starorfs$unique_gene_symbol,
                "star",
                "triangle"
              )
            ),
            border.color = "black"
          )
      },
      error = function(err) {
        return(data.frame())
      }
    )
    if (nrow(nodeq2) == 0) {
      return()
    }
    visNetwork(nodes = nodeq2, edges = edgeq2, height = "1200px") %>%
      visLayout(randomSeed = 23) %>%
      visNodes(
        borderWidth = 2, color = list(
          border = "green",
          highlight = "yellow"
        ),
        font = list(size = 9)
      ) %>%
      visPhysics(
        stabilization = F,
        solver = "repulsion",
        enabled = F
      ) %>%
      visEdges(smooth = F) %>%
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
  })

  output$conn <- renderUI({
    if (input$doMod != T) {
      return()
    }
    outputtab <- outputtab()
    inid <- outputtab$unique_gene_symbol

    edgeq <- rv$act_net$edges %>% filter(from == inid | to == inid)
    nodeq <- c(edgeq$from, edgeq$to) %>% unique()
    rv$conn <<- tryCatch(
      {
        length(nodeq) - 1
      },
      error = function(err) {
        return(0)
      }
    )

    rank <- rv$act_temp4 %>%
      filter(gene == inid) %>%
      pull(rank)
    if (length(rank) == 0) {
      rank <- "NA"
    }
    mod <- rv$act_modules %>%
      filter(gene == inid) %>%
      pull(module_n)
    if (length(mod) == 0) {
      mod <- "low expression"
    }
    r <- rv$act_modules %>%
      filter(gene == inid) %>%
      pull(r)
    if (length(r) == 0) {
      r <- ""
    } else {
      r <- as.character(round(r, digits = 2))
      r <- str_c(", r = ", r)
    }
    maxrank <- rv$act_modules %>%
      filter(module_n == mod) %>%
      nrow()
    HTML(str_c(
      "# of connections (",
      input$region, "): ",
      rv$conn, "<br>",
      rank, " out of ",
      maxrank,
      " in module ",
      mod,
      r
    ))
  })

  observeEvent(input$current_node_id, {
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = input$current_node_id,
      choices = autocomplete_list,
      server = T
    )
  })

  outputtab <- reactive({
    inid <- inid()
    
    if (inid %in% orfs$gene_id | inid %in% orfs$unique_gene_symbol) {
      temp_orfs <- orfs %>% filter(gene_id == inid | unique_gene_symbol == inid)
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
      select(1:6) %>%
      unique()
    
    if (nrow(filtered) == 0) {
      rv$blast <<- ""
      rv$pval <<- data.frame()
      rv$temp_orfs <<- data.frame()
      return(NULL)
    }
    
    filtered2 <- bed %>%
      filter(gene_id == filtered$gene_id[1]) %>%
      select(c(1, 2, 3, 6))
    if (nrow(filtered2) > 1) {
      filtered2 <- cbind(bed_merge(filtered2), filtered2[1, 4])
    }

    if (nrow(filtered2 > 0)) {
      tempvec <- unique(c(filtered$unique_gene_symbol, historytab))
      if (length(tempvec) > 10) {
        tempvec <- tempvec[1:10]
      }
      historytab <<- tempvec
    }

    out <- cbind(filtered, filtered2)

    out
  })

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
      rv$temp_orfs %>% select(1:6)
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
      temp1 <- domains %>% filter(gene_id == outputtab$gene_id)
    }
    temp1
  })

  # download ucscplot
  output$ucscPlot <- renderUI({
    if (input$doUcsc != T) {
      return()
    }
    outputtab <- outputtab()
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
      '">,',
      '<img src="',
      src,
      '">',
      "<a/>"
    ))
  })

  # various links in sidebar
  output$tab <- renderUI({
    outputtab <- outputtab()
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

  output$tab2 <- renderUI({
    outputtab <- outputtab()
    clean <- a(outputtab$unique_gene_symbol,
      href = str_c(
        "https://www.ncbi.nlm.nih.gov/gene/?term=",
        outputtab$clean_gene_symbol,
        "[sym]+AND+human[ORGN]"
      )
    )
    tagList("genbank:", clean)
  })

  output$tab3 <- renderUI({
    outputtab <- outputtab()
    clean <- a(outputtab$unique_gene_symbol,
      href = str_c(
        "https://www.genenames.org/data/gene-symbol-report/#!/symbol/",
        outputtab$clean_gene_symbol
      )
    )
    tagList("hgnc:", clean)
  })

  output$tab4 <- renderUI({
    outputtab <- outputtab()
    clean <- a(outputtab$unique_gene_symbol,
      href = str_c(
        "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
        outputtab$clean_gene_symbol
      )
    )
    tagList("genecard:", clean)
  })

  output$blastlink <- renderUI({
    if (rv$blast != "" & !(is.na(rv$blast))) {
      outputtab <- outputtab()
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
      sym <- tryCatch(outputtab()$unique_gene_symbol, error = function(err) {
        return("wrong")
      })
      paste0(sym, ".pdf", sep = "")
    },
    content = function(file) {
      if (input$doTis == T) {
        w <- 8
        h <- 6
      } else {
        w <- 8
        h <- 3
      }
      ggplot2::ggsave(file, plot = boxPlot1(), device = "pdf", width = w, height = h)
    }
  )

  # history
  output$history1 <- renderUI({
    outputtab()
    actionLink("history1", label = historytab[1])
  })
  output$history2 <- renderUI({
    outputtab()
    actionLink("history2", label = historytab[2])
  })
  output$history3 <- renderUI({
    outputtab()
    actionLink("history3", label = historytab[3])
  })
  output$history4 <- renderUI({
    outputtab()
    actionLink("history4", label = historytab[4])
  })
  output$history5 <- renderUI({
    outputtab()
    actionLink("history5", label = historytab[5])
  })
  output$history6 <- renderUI({
    outputtab()
    actionLink("history6", label = historytab[6])
  })
  output$history7 <- renderUI({
    outputtab()
    actionLink("history7", label = historytab[7])
  })
  output$history8 <- renderUI({
    outputtab()
    actionLink("history8", label = historytab[8])
  })
  output$history9 <- renderUI({
    outputtab()
    actionLink("history9", label = historytab[9])
  })
  output$history10 <- renderUI({
    outputtab()
    actionLink("history10", label = historytab[10])
  })
  observeEvent(input$history1, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[1],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history2, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[2],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history3, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[3],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history4, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[4],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history5, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[5],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history6, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[6],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history7, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[7],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history8, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[8],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history9, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[9],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$history10, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytab[10],
      choices = autocomplete_list,
      server = T
    )
  })
  observeEvent(input$geneID != "", {
    if (rv$run2 == 1 & input$geneID != "" & !is.null(input$geneID) & input$tabMain == "plot") {
      rv$run2 <- 0
      click("Find")
    }
  })

  onclick(
    "geneID",
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = "",
      choices = autocomplete_list,
      server = T
    )
  )

  # link to trait pdf
  onclick("conn", {
    filename <- str_c(input$region, "_trait.pdf")
    output$pdfview <- renderText({
      return(paste('<iframe style="height:800px; width:100%" src="',
        filename,
        '"></iframe>',
        sep = ""
      ))
    })
  })

  # loading list and viewing
  observeEvent(input$file, {
    rv$listn <- 0
    path <- input$file$datapath
    historytablist <<- read_csv(path, col_names = FALSE) %>% pull(1)
  })
  
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
  
  onclick("Prev", {
    rv$listn <- rv$listn - 1
    if (rv$listn < 1) {
      rv$listn <- 1
    }
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytablist[rv$listn],
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
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = historytablist[rv$listn],
      choices = autocomplete_list,
      server = T
    )
  })

  output$listn <- renderUI({
    HTML(str_c(rv$listn, "of ", length(historytablist)))
  })

  output$tbl <- DT::renderDataTable({
    DT::datatable(bed %>% select(
      unique_gene_symbol,
      everything()
    ),
    filter = "top",
    escape = FALSE,
    selection = "none",
    rownames = FALSE
    )
  })

  output$saveFiltered <- downloadHandler("filtré.csv", content = function(file) {
    s <- input$tbl_rows_all
    write_csv((bed %>% select(unique_gene_symbol, everything()))[s, ], file)
  })
  
  output$tbl2 <- DT::renderDataTable({
    DT::datatable(imp %>% select(
      unique_gene_symbol,
      contains("rank"),
      everything()
    ),
    filter = "top",
    escape = FALSE,
    selection = "single",
    rownames = FALSE
    )
  })
  
  output$saveFiltered2 <- downloadHandler("filtré2.csv", content = function(file) {
    s <- input$tbl_rows_all
    write_csv((imp %>% select(unique_gene_symbol, contains("rank"), everything()))[s, ], file)
  })
  
  observeEvent(input$tbl2_rows_selected, {
    rv$run2 <- 1
    updateSelectizeInput(session,
                         inputId = "geneID",
                         selected = imp[input$tbl2_rows_selected, "unique_gene_symbol"],
                         choices = autocomplete_list,
                         server = T
    )
  })
  
  output$tbl3 <- DT::renderDataTable({
    DT::datatable(
      rfvars,
      filter = "top",
      escape = FALSE,
      selection = "single",
      rownames = FALSE
    )
  })
  
  observeEvent(input$tbl3_rows_selected, {
    rv$run2 <- 1
    updateSelectizeInput(session,
                         inputId = "geneID",
                         selected = rfvars[input$tbl3_rows_selected, "gene"],
                         choices = autocomplete_list,
                         server = T
    )
  })
  
  output$oob <- renderPlotly({
    ggplot(trf, aes(x = Number.Variables, y = OOB)) +
      geom_point() +
      theme_cowplot()
  })

  output$oob2 <- renderPlotly({
    g <- ggplot(trf, aes(x = Number.Variables, y = OOB)) +
      geom_point() +
      theme_cowplot() +
      xlim(0,50) 
    ggplotly(g, source = "oob2", selectedpoints=list(9)) %>% 
      layout(xaxis = list(showspikes = TRUE))
  }) 
  
  plotlysel <- reactive({
    event_data("plotly_click", source = "oob2")
  })
  
  output$sel <- renderUI({
    rv$xsel <<- plotlysel()
    paste0("selected: ", as.character(rv$xsel$x))})
  
  output$saveFiltered3 <- downloadHandler("var_list.txt", content = function(file) {
    write_lines(vars_set(rfvars, rv$xsel$x), file)
  })
  
  observeEvent(input$back_to_top, {
    shinyjs::runjs("window.scrollTo(0, 0)")
  }, ignoreNULL = T)
}

# Run the application
shinyApp(ui = ui, server = server, enableBookmarking = "url")
