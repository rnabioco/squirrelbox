library(shiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(magrittr)
library(readr)
library(stringr)
library(shinyjs)

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
combined <- read_tsv("combined.tsv.gz")

# reorder factors for correct display order
combined %<>% mutate(state = factor(state, levels = state_order), region = factor(region, levels = region_order))

# read annotation file to find ucsc track
bed <- read_tsv("final_annot_bed12_20_sort.bed12", col_names = F) %>% select(1:6)
colnames(bed) <- c("chrom", "start", "end", "unique_gene_symbol", "score", "strand")

# empty history list to start
historytab <- c()

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
    sidebarPanel(div(style="display: inline-block;vertical-align:top; width: 200px;",tagAppendAttributes(textInput("geneID", label = NULL, value = "ENSSTOG00000002411"), `data-proxy-click` = "Find")), 
                 div(style="display: inline-block;vertical-align:top; width: 10px;",actionButton("Find", "Find")),
                 uiOutput("tab"), uiOutput("tab2"), br(), br(),
                 uiOutput("history1"),uiOutput("history2"),uiOutput("history3"),uiOutput("history4"),uiOutput("history5"),
                 tableOutput("print")),
    mainPanel(plotOutput("boxPlot", width = 800, height = 600),
              br(), br(),
              tableOutput("results"))
  )
)

# Define server logic required to draw the boxplot and render metadata table
server <- function(input, output, session) {
  
  inid <- eventReactive(input$Find, {
    input$geneID
  }, ignoreNULL = FALSE)
  
  rv <- reactiveValues()
  rv$run2 <- 0
  
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
    
    tempvec <- unique(c(filtered$unique_gene_symbol, historytab))
    if (length(tempvec) > 5) {
      tempvec <- tempvec[1:5]
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

  observeEvent(input$history1, {
    updateTextInput(session, inputId = "geneID", value = historytab[1])
    rv$run2 <- 1
  })
  observeEvent(input$history2, {
    updateTextInput(session, inputId = "geneID", value = historytab[2])
    rv$run2 <- 1
  })
  observeEvent(input$history3, {
    updateTextInput(session, inputId = "geneID", value = historytab[3])
    rv$run2 <- 1
  })
  observeEvent(input$history4, {
    updateTextInput(session, inputId = "geneID", value = historytab[4])
    rv$run2 <- 1
  })
  observeEvent(input$history5, {
    updateTextInput(session, inputId = "geneID", value = historytab[5])
    rv$run2 <- 1
  })
  observeEvent(input$geneID, {
    if (rv$run2 == 1) {
      click("Find")
      rv$run2 <- 0
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
