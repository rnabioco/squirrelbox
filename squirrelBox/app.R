library(shiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(magrittr)
library(readr)
library(stringr)

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

# Define UI for application that draws the boxplot
ui <- fluidPage(
  
  titlePanel("13-lined ground squirrel gene-level RNA-seq expression by tissue"),
  sidebarLayout(
    sidebarPanel(textInput("geneID",label = "Gene ID", value = "ENSSTOG00000002411"),
                 uiOutput("tab")),
    mainPanel(plotOutput("boxPlot", width = 800, height = 600),
              br(), br(),
              tableOutput("results"))
  )
)

# Define server logic required to draw the boxplot and render metadata table
server <- function(input, output) {
  
  output$boxPlot <- renderPlot({
    ggplot(combined %>% filter(gene_id == input$geneID | unique_gene_symbol == input$geneID), aes(state, log2_counts)) +
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
      filter(gene_id == input$geneID | unique_gene_symbol == input$geneID) %>% 
      select(1:6) %>%
      unique()
    filtered2 <-
      bed %>%
      filter(unique_gene_symbol == filtered$unique_gene_symbol[1] | unique_gene_symbol == input$geneID) %>%
      select(c(1,2,3,6))
    cbind(filtered, filtered2)
  })
  
  output$results <- renderTable({
    outputtab()
  })
  
  output$tab <- renderUI({
    outputtab <- outputtab()
    url <- a(outputtab$unique_gene_symbol, href=str_c("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hub_209779_KG_HiC&position=", outputtab$chrom, ":", outputtab$start, "-", outputtab$end))
    tagList("trackhub link:", url)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
