library(shiny)
library(dplyr)
library(ggplot2)
library(cowplot)
# library(magrittr)
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
options(stringsAsFactors = FALSE)
theme_set(theme_cowplot())

find_padj <- function(region, state) {
  temp <- str_c(padj2[str_sub(padj2,1,1) == str_to_lower(str_sub(region,1,1)) & str_detect(padj2, state)], collapse = "\n")
  if (length(temp) == 0){
    temp <- "NA"
  }
  temp
}

# define state colors and order, region order
state_cols <-  c(
  SA = rgb(255, 0, 0, maxColorValue=255),
  IBA = rgb(67, 205, 128, maxColorValue=255),
  Ent = rgb(155, 48, 255, maxColorValue=255),
  LT = rgb(25, 25, 112, maxColorValue=255),
  Ar = rgb(0, 0, 255, maxColorValue=255),
  SpD = rgb(255, 165, 0, maxColorValue=255) 
)

state_order = c(
  "SA", "IBA", "Ent", "LT", "Ar", "SpD"
)

region_order = c(
  "Forebrain", "Hypothalamus", "Medulla", "Adrenal", "Kidney", "Liver"
)

usevroom <- F #!require(vroom)
# read database

if (file.exists("combined2.csv")) {
  combined2 <- read_csv("combined2.csv", col_types = cols())
  combined3 <- read_csv("combined3.csv", col_types = cols())
  combined <- combined3 %>% left_join(combined2, by = "gene_id")
} else if (file.exists("combined.tsv")) {
  if (!usevroom) {
    combined <- read_tsv("combined.tsv", col_types = cols())
  } else {
    combined <- vroom::vroom("combined.tsv", altrep_opts = T)
  }
} else {
  combined <- read_tsv("combined.tsv.gz", col_types = cols())
}

# combined2 <- combined %>% select(gene_id, sample:region)
# combined3 <- combined %>% select(1:6) %>% distinct()
# c3 <- combined %>% unite(comb, sample, state, region) %>% spread(comb, log2_counts)
# write_csv(c3, "combined_wide.csv")
# combined <- read_csv("combined_wide.csv")
# combined <- combined %>% gather(key = "comb", value = "log2_counts", -(gene_id:source))
# combined <- combined %>% filter(!(is.na(log2_counts))) %>% separate("comb", into = c("sample", "state", "region"))

# reorder factors for correct display order
combined <- combined %>% mutate(state = factor(state, levels = state_order),
                     region = factor(region, levels = region_order))

# stash smaller files
# if (file.exists("env.rdata")) {
#   load("env.rdata")
# } else {
  
autocomplete_list <- str_sort(c(combined$unique_gene_symbol, combined$gene_id) %>%
                                unique(), decreasing = T)

# read annotation file to find ucsc track
if (!usevroom) {
  bed <- read_tsv("final_annot_bed12_20_sort.bed12", col_names = F, col_types = cols()) %>%
    select(1:6,13)
} else {
  bed <- vroom::vroom("final_annot_bed12_20_sort.bed12", col_names = F, altrep_opts = T)
  bed <- bed %>% select(1:6,13)
  vroom:::vroom_materialize(bed)
}

colnames(bed) <- c("chrom", "start", "end", "unique_gene_symbol", "score", "strand", "gene_id")

# empty history list to start
historytab <- c()

hy_modules <- suppressWarnings(read_csv("201907_hy_modules.csv", col_types = cols()))
med_modules <- suppressWarnings(read_csv("201907_med_modules.csv", col_types = cols()))
fore_modules <- suppressWarnings(read_csv("201907_fore_modules.csv", col_types = cols()))

# read node igraph object
hy_ig <- readRDS("201907hy_conn_list")
med_ig <-readRDS("201907med_conn_list")
fore_ig <- readRDS("201907fore_conn_list")
hy_net <- toVisNetworkData(hy_ig)
med_net <- toVisNetworkData(med_ig)
fore_net <- toVisNetworkData(fore_ig)
hy_temp4 <- as_edgelist(hy_ig) 
hy_temp4 <- c(hy_temp4[,1],hy_temp4[,2]) %>%
  table() %>%
  as.tibble()
colnames(hy_temp4) <- c("gene", "n")
med_temp4 <- as_edgelist(med_ig) 
med_temp4 <- c(med_temp4[,1],med_temp4[,2]) %>%
  table() %>%
  as.tibble()
colnames(med_temp4) <- c("gene", "n")
fore_temp4 <- as_edgelist(fore_ig) 
fore_temp4 <- c(fore_temp4[,1],fore_temp4[,2]) %>%
  table() %>%
  as.tibble()
colnames(fore_temp4) <- c("gene", "n")

hy_temp4 <- hy_temp4 %>%
  left_join(hy_modules, by = "gene") %>%
  select(-X1) %>% group_by(module_n) %>%
  mutate(rank = rank(-n)) %>%
  ungroup()
med_temp4 <- med_temp4 %>%
  left_join(med_modules, by = "gene") %>%
  select(-X1) %>% group_by(module_n) %>%
  mutate(rank = rank(-n)) %>%
  ungroup()
fore_temp4 <- fore_temp4 %>%
  left_join(fore_modules, by = "gene") %>%
  select(-X1) %>% group_by(module_n) %>%
  mutate(rank = rank(-n)) %>%
  ungroup()

# eigengene plots
hy_gg <- readRDS("hy_gg")
fore_gg <- readRDS("fore_gg")
med_gg <- readRDS("med_gg")

# read go terms
# gos <- read_tsv("/Users/rf/Downloads/goa_human.gaf", comment = "!", col_names = F) %>% dplyr::select(clean_gene_symbol = X3, go_id = X5)
# gos <- readRDS("sq_hm_mart") %>% mutate(clean_gene_symbol = str_to_title(hgnc_symbol))
gmt_to_list <- function(path,
                        cutoff = 0,
                        sep = "\thttp://www.broadinstitute.org/gsea/msigdb/cards/.*?\t",
                        rm = "REACTOME_") {
  df <- readr::read_csv(path,
                        col_names = F, col_types = cols())
  df <- tidyr::separate(df,
                        X1,
                        sep = sep,
                        into = c("path", "genes"))
  df <- dplyr::mutate(df, 
                      path = stringr::str_replace(df$path,
                                                      rm,
                                                      ""),
                      genes = stringr::str_split(genes,
                                                    "\t"))
  tidyr::unnest(df, genes)
}
gmt <- gmt_to_list("c5.all.v6.2.symbols.gmt", rm = "^GO_")

refs <- combined %>% distinct(unique_gene_symbol, clean_gene_symbol)
TFs <- gmt %>% filter(str_detect(path, str_to_upper("transcription_repressor_activity|transcription_factor_activity")), 
                      path != "VIRAL_TRANSCRIPTION") %>%
  pull(genes) %>%
  unique()

refTFs <- refs %>% mutate(clean_gene_symbol = str_to_upper(clean_gene_symbol)) %>%
  filter(clean_gene_symbol %in% str_to_upper(TFs)) %>%
  pull(unique_gene_symbol)

# load orf predictions
if (!usevroom) {
  orfs <- read_csv("padj_orf.csv", col_types = cols()) %>% select(gene_id, orf_len = len, exons, rna_len = transcript, orf, unique_gene_symbol, everything())
} else {
  orfs <- vroom::vroom("padj_orf.csv", altrep_opts = T) %>% select(gene_id, orf_len = len, exons, rna_len = transcript, orf, unique_gene_symbol, everything())
  vroom:::vroom_materialize(orfs)
}
starorfs <- orfs %>% group_by(unique_gene_symbol) %>% arrange(desc(orf)) %>% dplyr::slice(1) %>% filter(exons > 1, str_detect(unique_gene_symbol, "^G[0-9]+|_"), orf_len >= 100)
domains <- read_csv("novel_domains.csv", col_types = cols())

# lists <- ls(all.names = TRUE)
# save(list = lists[lists != "combined"], file = "env.rdata", envir = environment(), compress = F)
# }

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
  titlePanel("13-lined ground squirrel gene-level RNA-seq expression by tissue"),
  sidebarLayout(
    sidebarPanel(style = "position:fixed;width:inherit;",
                 div(style="display: inline-block;vertical-align:top; width: 200px;",tagAppendAttributes(selectizeInput("geneID", label = NULL, selected = "G8462", choices = "G8462"), `data-proxy-click` = "Find")), 
                 div(style="display: inline-block;vertical-align:top; width: 10px;",actionButton("Find", "Find")),
                 tags$hr(style="border-color: green;"),
                 checkboxInput("doPlotly", "interactive padj", value = F, width = NULL),
                 checkboxInput("doTis", "plot non-brain", value = F, width = NULL),
                 checkboxInput("doEigen", "plot eigengenes", value = F, width = NULL),
                 checkboxInput("doUcsc", "download track", value = F, width = NULL),
                 checkboxInput("doMod", "find module", value = T, width = NULL),
                 checkboxInput("doNet", "plot network", value = T, width = NULL),
                 checkboxInput("doKegg", "GO terms", value = T, width = NULL),
                 br(),
                 selectInput("region", label = NULL, choices = list("fore","hy","med"), selected = "hy"),
                 uiOutput("conn"),
                 tags$hr(style="border-color: green;"),
                 uiOutput("tab"), uiOutput("blastlink"), 
                 uiOutput("tab2"), uiOutput("tab3"), uiOutput("tab4"),
                 downloadButton('savePlot', label = "Download plot"),
                 # br(),
                 tags$hr(style="border-color: green;"),
                 uiOutput("history1"),uiOutput("history2"),uiOutput("history3"),uiOutput("history4"),uiOutput("history5"),uiOutput("history6"),uiOutput("history7"),uiOutput("history8"),uiOutput("history9"),uiOutput("history10")),
    mainPanel(uiOutput('boxPlotUI'),
              bsModal('modalPDF', title = "module-trait", trigger = "conn", size = "large", htmlOutput("pdfview")),
              #plotOutput("boxPlot", width = 800, height = 600),
              uiOutput('hyEigenPlot'),
              tags$hr(style="border-color: green;"),
              tableOutput("results"),
              tableOutput("orfinfo"),
              tags$hr(style="border-color: green;"),
              htmlOutput("ucscPlot"),
              tags$hr(style="border-color: green;"),
              visNetworkOutput("connPlot"),
              tags$hr(style="border-color: green;"),
              tableOutput("gotab"))
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
  rv$act_modules <- hy_modules
  rv$act_temp4 <- hy_temp4
  rv$act_net <- hy_net
  
  observeEvent(input$region, {
    if (input$region == "fore") {
      rv$act_modules <- fore_modules
      rv$act_temp4 <- fore_temp4
      rv$act_net <- fore_net
    } else if (input$region == "hy") {
      rv$act_modules <- hy_modules
      rv$act_temp4 <- hy_temp4
      rv$act_net <- hy_net
    } else if (input$region == "med") {
      rv$act_modules <- med_modules
      rv$act_temp4 <- med_temp4
      rv$act_net <- med_net
    }
  })
  
  inid <- eventReactive(input$Find, {
    #rv$init <- rv$init + 1
    rv$old <- input$geneID
    shinyjs::runjs("window.scrollTo(0, 0)")
    input$geneID
  }, ignoreNULL = F)
  
  historytab <- c()
  
  boxPlot1 <- reactive({
    inid <- isolate(inid())
    plot_temp <- combined %>% filter(gene_id == inid | unique_gene_symbol == inid)
    if (input$doTis != T) {
      plot_temp <- plot_temp %>% filter(region %in% c("Forebrain", "Hypothalamus", "Medulla"))
    }
    if (nrow(rv$pval) != 0) {
      padj2 <- rv$pval %>% select(ends_with("_wald_padj")) %>% t()
      padj2 <<- str_c(rownames(padj2), format(padj2[,1], digits = 2), sep = " : ")
      plot_temp <- plot_temp %>% mutate(text = mapply(find_padj, as.character(region), as.character(state)))
    } else {
      plot_temp <- plot_temp %>% mutate(text = "NA")
    }
    
    set.seed(1)
    g <- ggplot(plot_temp, aes(state, log2_counts, text = text)) +
      geom_boxplot(aes(fill = state), outlier.shape = NA) +
      geom_jitter() +
      scale_fill_manual(values = state_cols) +
      ylab("rlog(counts)") +
      facet_wrap(~region, scales = "free") +
      theme(legend.position = "none")
    g
  })
  
  boxPlotr <- reactive({
    g <- boxPlot1()
    output$boxPlot <- renderPlot(g)
    if (input$doTis == T) {
      plotOutput('boxPlot', width = 800, height = 600)
    } else {
      plotOutput('boxPlot', width = 800, height = 300)
    }
  })
  
  boxPlotlyr <- reactive({
    g <- boxPlot1()
    output$boxPlot2 <- renderPlotly(ggplotly(g + facet_wrap(~region), tooltip = "text"))
    if (input$doTis == T) {
      plotlyOutput('boxPlot2', width = 800, height = 600)
    } else {
      plotlyOutput('boxPlot2', width = 800, height = 300)
    }
  })
  
  output$boxPlotUI <- renderUI({
    if (input$doPlotly == T) {
      boxPlotlyr()
    } else {
      boxPlotr()
    }
  })
  
  output$hyEigenPlot <- renderUI({
    eigenplotr()
  })
  
  eigenplotr <- reactive({
    if (input$doEigen != T) {
      g <- ""
      output$boxPlot3 <- renderPlot(g)
      plotOutput("boxPlot3", height = 1)
    } else {
      outputtab <- outputtab()
      inid <- outputtab$unique_gene_symbol
      hy_mod <- hy_modules %>% filter(gene == inid) %>% pull(module_n)
      med_mod <- med_modules %>% filter(gene == inid) %>% pull(module_n)
      fore_mod <- fore_modules %>% filter(gene == inid) %>% pull(module_n)
      
      g <- cowplot::plot_grid(fore_gg[[as.numeric(fore_mod) + 1]], hy_gg[[as.numeric(hy_mod) + 1]], med_gg[[as.numeric(med_mod) + 1]], ncol = 3)
      output$boxPlot3 <- renderPlot(g)
      plotOutput("boxPlot3", width = 800, height = 300)
    }
  })
  
  output$connPlot <- renderVisNetwork({
    if (input$doMod != T) {
      return()
    }
    outputtab <- outputtab()
    queryid <- outputtab$unique_gene_symbol
    edgeq <- rv$act_net$edges %>% filter(from == queryid | to == queryid)
    nodeq <- c(edgeq$from, edgeq$to) %>% unique()
    edgeq2 <- tryCatch({rv$act_net$edges %>% filter(from %in% nodeq | to %in% nodeq) %>% mutate(color = "gray", opacity = 0, width = 0)}, error = function(err) {
      return(data.frame())
    })
    nodeq2 <- tryCatch({rv$act_net$nodes[nodeq,] %>% left_join(table(c(edgeq2$from, edgeq2$to)) %>% as.data.frame(stringsAsFactors = F), by = c("id" = "Var1")) %>% mutate(value = Freq, color = ifelse(id == queryid, "red", "lightblue"), shape = ifelse(id %in% refTFs, "square", ifelse(id %in% starorfs$unique_gene_symbol, "star", "triangle")), border.color = "black")}, error = function(err) {
    return(data.frame())
  })
    rv$conn <<- tryCatch({nodeq2 %>% filter(id == queryid) %>% pull(value)}, error = function(err) {
      return(0)
    })
    if (input$doNet != T | nrow(nodeq2) == 0) {
      return()
    }
    visNetwork(nodes = nodeq2, edges = edgeq2, height = "1200px") %>%
      visLayout(randomSeed = 23) %>% 
      visNodes(borderWidth = 2, color = list(border = "green", highlight = "yellow"), font = list(size = 9)) %>% 
      visPhysics(stabilization = F, solver = "repulsion", enabled = F) %>%
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
    rank <- rv$act_temp4 %>% filter(gene == inid) %>% pull(rank)
    if (length(rank) == 0) {
      rank <- "NA"
    }
    mod <- rv$act_modules %>% filter(gene == inid) %>% pull(module_n)
    if (length(mod) == 0) {
      mod <- "low expression"
    }
    maxrank <- rv$act_modules %>% filter(module_n == mod) %>% nrow()
    HTML(str_c("# of connections (", input$region, "): ", rv$conn, "<br>", rank, " out of ", maxrank, " in module ", mod))
  })
  
  observeEvent(input$current_node_id, {
    updateSelectizeInput(session, inputId = "geneID", selected = input$current_node_id, choices = autocomplete_list, server = T)
  })
  
  outputtab <- reactive({
    inid <- inid()
    filtered <- combined %>% filter(gene_id == inid | unique_gene_symbol == inid) %>% 
      select(1:6) %>%
      unique()
    filtered2 <- bed %>% filter(gene_id == filtered$gene_id[1]) %>%
      select(c(1,2,3,6))
    if (nrow(filtered2) > 1) {
      filtered2 <- cbind(bed_merge(filtered2), filtered2[1,4])
    }
    
    tempvec <- unique(c(filtered$unique_gene_symbol, historytab))
    if (length(tempvec) > 10) {
      tempvec <- tempvec[1:10]
    }
    historytab <<- tempvec
    
    cbind(filtered, filtered2)
  })
  
  output$results <- renderTable({
    outputtab()
  }, digits = 0)
  
  # orf call table
  output$orfinfo <- renderTable({
    inid <- outputtab()$gene_id
    if (inid %in% orfs$gene_id) {
      temp_orfs <- orfs %>% filter(gene_id == inid)
      if (nrow(temp_orfs) == 0) {
        rv$blast <<- ""
        rv$pval <<- data.frame()
        temp_orfs <- data.frame()
        return(temp_orfs)
      } else {
        rv$pval <<- temp_orfs
        rv$blast <<- temp_orfs$orf[1]
      }
    } else {
      rv$blast <<- ""
      rv$pval <<- data.frame()
      temp_orfs <- data.frame()
      return(temp_orfs)
    }
    temp_orfs %>% select(1:6)
  }, digits = 0)
  
  # goterm table
  output$gotab <- renderTable({
    if (input$doKegg != T) {
      return()
    }
    outputtab <- outputtab()
    temp1 <- gmt %>% filter(genes == str_to_upper(outputtab$clean_gene_symbol)) #%>% arrange(namespace_1003) %>% select(clean_gene_symbol, go_term = name_1006, go_type = namespace_1003, go_id)
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
    url <- str_c("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hub_209779_KG_HiC&hubUrl=http://amc-sandbox.ucdenver.edu/User11/tls10k/hub.txt&position=", outputtab$chrom, ":", outputtab$start - l, "-", outputtab$end + l)
    src <- str_c("http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hub_209779_KG_HiC&hubUrl=http://amc-sandbox.ucdenver.edu/User11/tls10k/hub.txt&ignoreCookie=1&position=", outputtab$chrom, ":", outputtab$start - l, "-", outputtab$end + l)
    HTML(str_c('<a href ="', url, '">,', '<img src="',src,'">', '<a/>'))
  })
  
  # various links in sidebar
  output$tab <- renderUI({
    outputtab <- outputtab()
    url <- a(outputtab$unique_gene_symbol, href=str_c("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hub_209779_KG_HiC&hubUrl=http://amc-sandbox.ucdenver.edu/User11/tls10k/hub.txt&position=", outputtab$chrom, ":", outputtab$start, "-", outputtab$end))
    tagList("trackhub:", url)
  })
  
  output$tab2 <- renderUI({
    outputtab <- outputtab()
    clean <- a(outputtab$unique_gene_symbol, href=str_c("https://www.ncbi.nlm.nih.gov/gene/?term=", outputtab$clean_gene_symbol, "[sym]+AND+human[ORGN]"))
    tagList("genbank:", clean)
  })
  
  output$tab3 <- renderUI({
    outputtab <- outputtab()
    clean <- a(outputtab$unique_gene_symbol, href=str_c("https://www.genenames.org/data/gene-symbol-report/#!/symbol/", outputtab$clean_gene_symbol))
    tagList("hgnc:", clean)
  })
  
  output$tab4 <- renderUI({
    outputtab <- outputtab()
    clean <- a(outputtab$unique_gene_symbol, href=str_c("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", outputtab$clean_gene_symbol))
    tagList("genecard:", clean)
  })
  
  output$blastlink <- renderUI({
    if (rv$blast != "" & !(is.na(rv$blast))) {
      outputtab <- outputtab()
      orf <- rv$blast
      url <- a(outputtab$unique_gene_symbol, href=str_c("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastp&DATABASE=nr&QUERY=", orf))
      tagList("blast orf:", url)
      } else {return()}
  })
  
  # save plot as pdf
  output$savePlot <- downloadHandler(
    filename = function() {
      sym <- tryCatch(outputtab()$unique_gene_symbol, error = function(err) {return("wrong")})
      paste0(sym,".pdf",sep = "")
    },
    content = function(file) {
      if (input$doTis == T) {
        w = 8
        h = 6
      } else {
        w = 8
        h= 3
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
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[1], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history2, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[2], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history3, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[3], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history4, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[4], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history5, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[5], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history6, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[6], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history7, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[7], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history8, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[8], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history9, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[9], choices = autocomplete_list, server = T)
  })
  observeEvent(input$history10, {
    rv$run2 <- 1
    updateSelectizeInput(session, inputId = "geneID", selected = historytab[10], choices = autocomplete_list, server = T)
  })
  observeEvent(input$geneID, {
    if (rv$run2 == 1 & input$geneID != "" & !is.null(input$geneID)) {
      click("Find")
      rv$run2 <- 0
    }
  })
  
  # handling enter and initialization
  observeEvent(rv$init, {
    if (rv$init == 0) {
      updateSelectizeInput(session, inputId = "geneID", selected = "G8462", choices = autocomplete_list, server = T)
      rv$init <- 2
    }
  })
  onclick("geneID", 
          updateSelectizeInput(session, inputId = "geneID", selected = "", choices = autocomplete_list, server = T)
  )
  
  # link to trait pdf
  onclick("conn",{
    filename <- str_c("201907_", input$region, "_trait.pdf")
    # output$pdfview <- renderUI({
    #   tags$iframe(style="height:800px; width:100%", src=filename)
    # })
    output$pdfview <-renderText({
      return(paste('<iframe style="height:800px; width:100%" src="',filename, '"></iframe>', sep = ""))
    })
  }
          )
}

# Run the application 
shinyApp(ui = ui, server = server)