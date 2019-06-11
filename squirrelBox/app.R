#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(magrittr)
library(readr)

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

# Define UI for application that draws the boxplot
ui <- fluidPage(

    titlePanel("13-lined ground squirrel gene-level RNA-seq expression by tissue"),
    sidebarLayout(
        sidebarPanel(textInput("geneID",label = "Gene ID", value = "ENSSTOG00000002411")),
        mainPanel(plotOutput("boxPlot", width = 800, height = 600),
                  br(), br(),
                  tableOutput("results"))
        )
    )

# Define server logic required to draw the boxplot and render metadata table
server <- function(input, output) {
    
    output$boxPlot <- renderPlot({
        ggplot(combined %>% filter(gene_id == input$geneID), aes(state, log2_counts)) +
            geom_boxplot(aes(fill = state)) +
            scale_fill_manual(values = state_cols) +
            ylab("rlog(counts)") +
            facet_wrap(~region, scales = "free") +
            theme(legend.position = "none")
    })
    
    output$results <- renderTable({
        filtered <-
            combined %>%
            filter(gene_id == input$geneID) %>% 
            select(1:6) %>%
            unique()
        filtered
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
