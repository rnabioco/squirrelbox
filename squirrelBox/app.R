library(shiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(magrittr)
library(readr)
library(stringr)
library(shinyjs)
library(visNetwork)
library(valr)
library(igraph)
library(tibble)

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

# read database
if (file.exists("combined.tsv")) {
  combined <- read_tsv("combined.tsv")
} else {
  combined <- read_tsv("combined.tsv.gz")
}

# reorder factors for correct display order
combined %<>% mutate(state = factor(state, levels = state_order), region = factor(region, levels = region_order))
autocomplete_list <- str_sort(c(combined$unique_gene_symbol, combined$gene_id) %>% unique(), decreasing = T)

# read annotation file to find ucsc track
bed <- read_tsv("final_annot_bed12_20_sort.bed12", col_names = F) %>% select(1:6)
colnames(bed) <- c("chrom", "start", "end", "unique_gene_symbol", "score", "strand")

# empty history list to start
historytab <- c()

# read node igraph object
hy_ig <- readRDS("201906hy_conn_list")
hy_net <- toVisNetworkData(hy_ig)
temp2 <- as_edgelist(hy_ig) 
temp3 <- c(temp2[,1],temp2[,2]) %>% table() %>% as.tibble()
colnames(temp3) <- c("gene", "n")
modules <- read.csv("201905_hy_modules.csv")
temp4 <- temp3 %>% left_join(modules, by = "gene") %>% select(-X) %>% group_by(module) %>% mutate(rank = rank(-n)) %>% ungroup()

# read go terms
gos <- readRDS("sq_hm_mart")
refs <- combined %>% distinct(unique_gene_symbol, clean_gene_symbol)
TFs <- gos %>% filter(str_detect(name_1006, "DNA-binding transcription activator activity")) %>% pull(hgnc_symbol) %>% unique()
refTFs <- refs %>% mutate(clean_gene_symbol = str_to_upper(clean_gene_symbol)) %>% filter(clean_gene_symbol %in% TFs) %>% pull(unique_gene_symbol)

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
  useShinyjs(),
  tags$head(tags$script(HTML(jscode))),
  titlePanel("13-lined ground squirrel gene-level RNA-seq expression by tissue"),
  sidebarLayout(
    sidebarPanel(style = "position:fixed;width:inherit;",
                 div(style="display: inline-block;vertical-align:top; width: 200px;",tagAppendAttributes(selectizeInput("geneID", label = NULL, selected = "ENSSTOG00000002411", choices = NULL), `data-proxy-click` = "Find")), 
                 div(style="display: inline-block;vertical-align:top; width: 10px;",actionButton("Find", "Find")),
                 tags$hr(style="border-color: green;"),
                 checkboxInput("doNet", "plot network", value = T, width = NULL),
                 uiOutput("conn"),
                 tags$hr(style="border-color: green;"),
                 uiOutput("tab"), uiOutput("tab2"), br(), br(),
                 uiOutput("history1"),uiOutput("history2"),uiOutput("history3"),uiOutput("history4"),uiOutput("history5"),uiOutput("history6"),uiOutput("history7"),uiOutput("history8"),uiOutput("history9"),uiOutput("history10")),
    mainPanel(plotOutput("boxPlot", width = 800, height = 600),
              tags$hr(style="border-color: green;"),
              tableOutput("results"),
              tags$hr(style="border-color: green;"),
              visNetworkOutput("connPlot"))
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
  
  inid <- eventReactive(input$Find, {
    if (rv$init == 0) {
      rv$init <- 1
    }
    rv$old <- input$geneID
    shinyjs::runjs("window.scrollTo(0, 0)")
    input$geneID
  }, ignoreNULL = F)
  
  historytab <- c()
  
  output$boxPlot <- renderPlot({
    set.seed(1)
    ggplot(combined %>% filter(gene_id == inid() | unique_gene_symbol == inid()), aes(state, log2_counts)) +
      geom_boxplot(aes(fill = state), outlier.shape = NA) +
      geom_jitter() +
      scale_fill_manual(values = state_cols) +
      ylab("rlog(counts)") +
      facet_wrap(~region, scales = "free") +
      theme(legend.position = "none")
  })
  
  output$connPlot <- renderVisNetwork({
    outputtab <- outputtab()
    queryid <- outputtab$unique_gene_symbol
    edgeq <- hy_net$edges %>% filter(from == queryid | to == queryid)
    nodeq <- c(edgeq$from, edgeq$to) %>% unique()
    edgeq2 <- tryCatch({hy_net$edges %>% filter(from %in% nodeq | to %in% nodeq) %>% mutate(color = "gray", opacity = 0, width = 0)}, error = function(err) {
      return(data.frame())
    })
    nodeq2 <- tryCatch({hy_net$nodes[nodeq,] %>% left_join(table(c(edgeq2$from, edgeq2$to)) %>% data.frame(), by = c("id" = "Var1")) %>% mutate(value = Freq, color = ifelse(id == queryid, "red", "lightblue"), shape = ifelse(id %in% refTFs, "square", "triangle"), border.color = "black")}, error = function(err) {
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
    outputtab <- outputtab()
    inid <- outputtab$unique_gene_symbol
    rank <- temp4 %>% filter(gene == inid) %>% pull(rank)
    if (length(rank) == 0) {
      rank <- "NA"
    }
    mod <- modules %>% filter(gene == inid) %>% pull(module)
    if (length(mod) == 0) {
      mod <- "low expression"
    }
    maxrank <- modules %>% filter(module == mod) %>% nrow()
    HTML(str_c("# of connections (hy): ", rv$conn, "<br>", rank, " out of ", maxrank, " in module ", mod))
  })
  
  observeEvent(input$current_node_id, {
    updateSelectizeInput(session, inputId = "geneID", selected = input$current_node_id, choices = autocomplete_list, server = T)
  })
  
  outputtab <- reactive({
    filtered <-
      combined %>%
      filter(gene_id == inid() | unique_gene_symbol == inid()) %>% 
      select(1:6) %>%
      unique()
    filtered2 <-
      bed %>%
      filter(unique_gene_symbol == filtered$unique_gene_symbol[1] | unique_gene_symbol == inid()) %>%
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
  })
  
  output$tab <- renderUI({
    outputtab <- outputtab()
    url <- a(outputtab$unique_gene_symbol, href=str_c("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hub_209779_KG_HiC&position=", outputtab$chrom, ":", outputtab$start, "-", outputtab$end))
    tagList("trackhub:", url)
  })
  
  output$tab2 <- renderUI({
    outputtab <- outputtab()
    clean <- a(outputtab$unique_gene_symbol, href=str_c("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", outputtab$clean_gene_symbol))
    tagList("genecard:", clean)
  })
  
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
    if (rv$run2 == 1 & input$geneID != "") {
      click("Find")
      rv$run2 <- 0
    }
  })
  observeEvent(rv$init == 1, {
    updateSelectizeInput(session, inputId = "geneID", selected = "Zfp36", choices = autocomplete_list, server = T)
    rv$run2 <- 1
  })
  onclick("geneID", 
          updateSelectizeInput(session, inputId = "geneID", selected = "", choices = autocomplete_list, server = T)
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
