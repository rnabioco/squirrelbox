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
        width: calc(100% - 70px);
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
        }",
      ".btn {
  text-transform: unset !important;
}"
    )
  )),
  useShinyjs(),
  introjsUI(),
  use_bs_tooltip(),
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
  tags$script(src = "shinyBS_mod.js"),
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
    style = "font-size:23px; background:white;text-decoration: underline; background-clip: inherit;"
  )),
  fixedPanel(
    style = "z-index:100;",
    right = 10,
    top = 22,
    tags$head(
      tags$style(HTML('#tutorial{background-color:gold;z-index:100;}'))
    ),
    introBox(
      actionButton("tutorial", "", icon = icon("question")) %>%
        bs_embed_tooltip("Take a tour through the app!", placement = "bottom"),
      data.step = 1,
      data.intro = "Welcome to the squirrelBox.<br><br>
      Please note that most buttons, tabs, and table columns have hover-over tips.",
      data.position = "left"
    )
  ),
  fixedPanel(
    style = "z-index:100;",
    actionButton("back_to_top", label = "to Top", icon = icon("angle-double-up")) %>% bs_embed_tooltip("scroll back to the top of the page"),
    bsButton("showpanel", "Sidebar", type = "toggle", value = FALSE, icon = icon("bars")) %>% bs_embed_tooltip("turn sidebar on/off"),
    right = 10,
    bottom = 10
  ),
  sidebarLayout(
    sidebarPanel(
      id = "SIDE",
      style = "position:fixed;width:23%;margin-top: 60px;z-index:50;",
      width = 3,
      div(
        id = "sideall",
        introBox(
          div(
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
              actionButton("Find", "Find", icon = icon("search")) %>% bs_embed_tooltip("gene id/symbols accepted", placement = "right")
            )
          ),
          data.step = 2,
          data.intro = "Query individual genes by id or symbol.",
          data.position = "bottom"
        ),
        div(
          id = "sidediv",
          tabsetPanel(
            id = "side1",
            tabPanel(
              span(icon("link", class = NULL, lib = "font-awesome"), "Gene_links", title = "additional external links for querying a specific gene"),
              br(.noWS = "outside"),
              introBox(
                fluidRow(
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
                data.step = 5,
                data.intro = "Other external links for the query gene.",
                data.position = "right"
              ),
              br(.noWS = "outside"),
              fluidRow(
                column(
                  width = 4,
                  downloadButton("savePlot", label = "PLOT", class = "download_this") %>%
                    bs_embed_tooltip("save current plot as .pdf", placement = "bottom")
                ),
                column(
                  width = 1
                ),
                column(
                  width = 4,
                  actionButton("Add", "to Cart") %>% bs_embed_tooltip("add current query gene to cart", placement = "bottom")
                )
              ),
              br(.noWS = "outside"),
              style = "height:150px;"
            ),
            tabPanel(
              span(icon("cogs", class = NULL, lib = "font-awesome"), "Plot_options", title = "options specific to each main tab"),
              br(.noWS = "outside"),
              div(id = "doEigendiv", checkboxInput("doEigen", "plot model clusters", value = T, width = NULL)),
              checkboxInput("doUcsc", "pull track", value = T, width = NULL),
              checkboxInput("doMod", "find module", value = T, width = NULL),
              checkboxInput("doKegg", "GO terms", value = T, width = NULL),
              checkboxInput("doTooltips", "show hover tips", value = T, width = NULL),
              div(id = "doBrdiv", checkboxInput("doBr", "plot brain data", value = T, width = NULL)),
              div(id = "doTisdiv", checkboxInput("doTis", "plot non-brain data", value = F, width = NULL)),
              div(
                id = "doLockdiv",
                checkboxInput("doLock", "lock main panel order", value = T, width = NULL) %>% 
                  bs_embed_tooltip("if unlocked, tabs, sections, and table columns can be dragged and reordered",
                                   placement = "right")
              ),
              tags$table(
                tags$head(
                  tags$style(HTML('#pval{margin-top: 0px; margin-bottom: -10px; font-size:12px;}'))
                ),
                tags$head(
                  tags$style(HTML('#plotw{margin-top: 0px; margin-bottom: -10px; font-size:12px;}'))
                ),
                tags$head(
                  tags$style(HTML('#ploth{margin-top: 0px; margin-bottom: -10px; font-size:12px;}'))
                ),
                tags$tr(width = "100%",
                        tags$td(width = "50%", div(style = "font-size:12px;", "p-value cutoff")),
                        tags$td(width = "50%", textInput("pval", NULL, value = sig_cut) %>%
                                  bs_embed_tooltip("p-value cut off used for all plotting/analyses", placement = "right"))),
                tags$tr(width = "100%",
                        tags$td(width = "50%", tags$div(style = "font-size:12px;", "plot height")),
                        tags$td(width = "50%", textInput("ploth", NULL, value = plot_height, width = "100px") %>%
                                  bs_embed_tooltip("plot height for app (px) and saved pdf (in), halved for box and line plots", placement = "right"))),
                tags$tr(width = "100%",
                        tags$td(width = "50%", tags$div(style = "font-size:12px;", "plot width")),
                        tags$td(width = "50%", textInput("plotw", NULL, value = plot_width, width = "100px") %>%
                                  bs_embed_tooltip("plot width for app (px) and saved pdf (in)", placement = "right"))),
                
              ),
              fluidRow(
                column(
                  width = 4,
                  downloadButton("savePlot2", label = "PLOT", class = "download_this") %>%
                    bs_embed_tooltip("save current plot as .pdf", placement = "bottom")
                ),
                column(
                  width = 1
                ),
                column(
                  width = 4,
                  downloadButton(outputId = "saveTable", label = "TABLE", class = "download_this") %>%
                    bs_embed_tooltip("save output/filtered table as .csv", placement = "bottom")
                )
              ),
              style = "height:150px;"
            )
          ),
          tabsetPanel(
            id = "side2",
            tabPanel(
              introBox(
                span(icon("file-alt", class = NULL, lib = "font-awesome"), "Genelist", title = "load list of genes for analysis from file or interactive table"),
                data.step = 8,
                data.intro = "Genelist can be loaded from external file, or passed from the tables/cart.<br><br>
                The other multi-gene analysis tabs, Lineplot/Heatmap/GO/Kmer, all use genes from this list.",
                data.position = "top"
              ),
              value = "load",
              div(id = "filediv", fileInput("file", label = NULL) %>%
                    bs_embed_tooltip("expects gene symbols as first column, or comma separated")),
              div(
                uiOutput("listn"),
                style = "display: inline-block;vertical-align:middle;"
              ),
              div(
                style = "display: inline-block;float:right;vertical-align:middle;",
                disabled(actionButton("Prev1", "Prev", icon = icon("angle-up")) %>%
                           bs_embed_tooltip("query previous gene on loaded list", placement = "bottom")),
                disabled(actionButton("Next1", "Next", icon = icon("angle-down")) %>%
                           bs_embed_tooltip("query next gene on loaded list", placement = "bottom"))
              ),
              
              DT::dataTableOutput("tbllist"),
              style = "height:300px; overflow-y: scroll;"
            ),
            tabPanel(
              value = "cart",
              introBox(
                span(icon("shopping-cart", class = NULL, lib = "font-awesome"), "Cart", title = "cart list of genes to save and export"),
                data.step = 9,
                data.intro = "A cart list stores query genes that were added (see button above), which can be exported to .txt file, or moved to the loaded Genelist.<br><br>
                Interactive clicking on the GO_enrich and Venn plots also put the corresponding genes into the cart.",
                data.position = "top"
              ),
              uiOutput("listn2"),
              fluidRow(
                column(
                  width = 4,
                  downloadButton(
                    outputId = "saveList",
                    label = "CART"
                  ) %>% bs_embed_tooltip("save genes in cart as .txt", placement = "bottom")
                ),
                column(
                  width = 1
                ),
                column(
                  width = 4,
                  actionButton("Load", "to Genelist") %>% bs_embed_tooltip("send genes in Cart to loaded Genelist in side panel", placement = "bottom")
                )
              ),
              DT::dataTableOutput("tbllist2"),
              style = "height:300px; overflow-y: scroll;"
            ),
            tabPanel(
              introBox(
                span(icon("history", class = NULL, lib = "font-awesome"), "History", title = "history list of query genes"),
                data.step = 10,
                data.intro = "A list of genes previously queried is also documented for review.",
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
          introBox(
            span(icon("pencil-ruler", class = NULL, lib = "font-awesome"),
                 "Gene_query",
                 title = "Plot expression box plot and other info of query gene"
            ),
            data.step = 3,
            data.intro = "For the query gene, this tab displays the expression boxplot, as well as other annotations and analyses.",
            data.position = "top"
          ),
          value = "Gene_query",
          div(
            id = "sorted",
            DT::dataTableOutput("results"),
            div(
              div(
                style = "display: inline-block;vertical-align:top;",
                introBox(
                  dropdownButton(
                    circle = FALSE, status = "options", icon = icon("gear"), width = "200px", size = "sm",
                    tooltip = tooltipOptions(title = "plotting options"), margin = "20px",
                    br(),
                    div(id = "doPlotlydiv", checkboxInput("doPlotly", "interactive plots", value = F, width = NULL) %>%
                          bs_embed_tooltip("display interactive plot with additional info on hover", placement = "right")),
                    div(id = "doPadjdiv", checkboxInput("doPadj", "indicate sig", value = T, width = NULL) %>%
                          bs_embed_tooltip("label groups by nonsignficance", placement = "right")),
                    div(id = "doNamediv", checkboxInput("doName", "additional labels", value = F, width = NULL) %>%
                          bs_embed_tooltip("label points by sample", placement = "right"))
                  ),
                  data.step = 4,
                  data.intro = "Additional plotting options, for interactivity and labels, can be accessed here.",
                  data.position = "left"
                ),
              ),
              div(
                style = "display: inline-block;vertical-align:top;",
                uiOutput("boxPlotUI") %>% withLoader()
                # conditionalPanel(
                #   condition = 'input.doPlotly == true',
                #   plotlyOutput('boxPlot_ly') %>% withLoader()
                # ),
                # 
                # conditionalPanel(
                #   condition = 'input.doPlotly == false',
                #   plotOutput('boxPlot_g') %>% withLoader()
                # )
              )
            ),
            bsCollapse(
              id = "tabs", multiple = TRUE, open = NULL, #open = "cluster_assignments",
              bsCollapsePanel(
                uiOutput("EigenPlot") %>% withLoader(),
                title = "cluster_assignments",
                style = "danger"
              )
            ),
            bsCollapse(
              id = "tabs2", multiple = TRUE, open = NULL, #open = "called_orfs",
              bsCollapsePanel(DT::dataTableOutput("orfinfo") %>% withLoader(),
                              title = "called_orfs",
                              style = "primary"
              )
            ),
            bsCollapse(
              id = "tabs3", multiple = TRUE, open = NULL, #open = "majiq_alternative_splicing",
              bsCollapsePanel(
                DT::dataTableOutput("majinfo") %>% withLoader(),
                title = "majiq_alternative_splicing",
                style = "warning"
              )
            ),
            bsCollapse(
              id = "tabs4", multiple = TRUE, open = NULL, #open = "UCSC browser plot",
              bsCollapsePanel(htmlOutput("ucscPlot") %>% withLoader(),
                              title = "UCSC browser plot",
                              style = "success"
              )
            ),
            bsCollapse(
              id = "tabs5", multiple = TRUE, open = NULL, #open = "go_terms/domains",
              bsCollapsePanel(DT::dataTableOutput("gotab") %>% withLoader(),
                              title = "go_terms/domains",
                              style = "info"
              )
            )
          )
        ),
        tabPanel(
          introBox(
            span(icon("table", class = NULL, lib = "font-awesome"),
                 "Transcript/Gene",
                 title = "Table of expression and other info of all genes/transcripts"
            ),
            data.step = 6,
            data.intro = "Here we summarize the genes in this study.<br><br>
            The table can be filtered, exported as .csv, and passed to Genelist for additional on-the-fly analyses.<br><br>
            Hover over column names for additional descriptions.",
            data.position = "top"
          ),
          value = "table_data",
          div(
            id = "doCollapsediv",
            style = "display: inline-block;width: 160px;",
            checkboxInput("doCollapse",
                          "longest transcript",
                          value = T,
                          width = NULL
            ) %>% bs_embed_tooltip("only show longest orf transcript for each gene", placement = "bottom")
          ),
          actionButton("loadtab", "to Genelist") %>%
            bs_embed_tooltip("send filtered results to loaded Genelist in side panel", placement = "right"),
          DT::dataTableOutput("genes")
        ),
        tabPanel(
          introBox(
            span(icon("list", class = NULL, lib = "font-awesome"), "Majiq_alt",
                 title = "Table of majiq output for alternative splicing events"
            ),
            data.step = 7,
            data.intro = "Similarly, splicing analysis via MAJIQ is presented as a table.<br><br>
            The table can be filtered, exported as .csv, and passed to Genelist for additional on-the-fly analyses.<br><br>
            Hover over column names for additional descriptions.",
            data.position = "top"
          ),
          value = "table_AS",
          div(
            id = "doJoindiv",
            style = "display: inline-block;width: 160px;",
            checkboxInput("doJoin",
                          "gene info",
                          value = FALSE,
                          width = NULL
            ) %>% bs_embed_tooltip("bring in gene info as last columns", placement = "bottom")
          ),
          actionButton("loadtab2", "to Genelist") %>%
            bs_embed_tooltip("send filtered results to loaded Genelist in side panel", placement = "right"),
          DT::dataTableOutput("alt")
        ),
        tabPanel(
          introBox(
            span(icon("chart-line", class = NULL, lib = "font-awesome"), "Line_plot",
                 title = "Plot expression of loaded Genelist"
            ),
            data.step = 11,
            data.intro = "Visualize loaded Genelist as line plot, also supports summarized line.",
            data.position = "top"
          ),
          value = "line_plot",
          div(
            style = "display: inline-block;vertical-align:top;",
            dropdownButton(
              circle = FALSE, status = "options", icon = icon("gear"), width = "200px", size = "sm",
              tooltip = tooltipOptions(title = "plotting options"), margin = "20px",
              br(),
              div(id = "doName2div", checkboxInput("doName2", "additional labels", value = F, width = NULL) %>%
                    bs_embed_tooltip("show toggleable legend", placement = "right")),
              div(id = "doNormdiv", checkboxInput("doNorm", "normalize to SA", value = F, width = NULL) %>%
                    bs_embed_tooltip("otherwise centered by mean expression", placement = "right")),
              div(id = "doSummaryiv", checkboxInput("doSummary", "summary line", value = F, width = NULL) %>%
                    bs_embed_tooltip("summarize instead of individual lines", placement = "right"))
            )
          ),
          div(
            style = "display: inline-block;vertical-align:top;",
            plotlyOutput("linePlot") %>% withLoader()
          )
        ),
        tabPanel(
          introBox(
            span(icon("th", class = NULL, lib = "font-awesome"), "Heatmap",
                 title = "Plot Z-Score of loaded Genelist as heat map"
            ),
            data.step = 12,
            data.intro = "Similar to the lineplot, visualize loaded Genelist as heatmap.",
            data.position = "top"
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
          introBox(
            span(icon("chart-bar", class = NULL, lib = "font-awesome"), "GO_terms",
                 title = "GO term enrichment for loaded Genelist (slow)"
            ),
            data.step = 13,
            data.intro = "GO term enrichment of loaded Genelist by fisher exact test.<br><br>
            Top 15 results are plotted, while full table can be exported.<br><br>
            Clicking on bar loads the corresponding genes into Cart.",
            data.position = "top"
          ),
          value = "enrichment_plot",
          plotlyOutput("richPlot") %>% withLoader()
        ),
        tabPanel(
          introBox(
            span(icon("kickstarter-k", class = NULL, lib = "font-awesome"), "Kmer",
                 title = "Kmer enrichment analysis and annotation for loaded Genelist (slow)"
            ),
            data.step = 14,
            data.intro = "Kmer analysis of loaded Genelist, with option to annotate known RBP motifs or mir seeds.<br><br>
            Note that this process may take ~30 seconds.",
            data.position = "top"
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
            radioButtons("kmlab", "annotate kmer", c("RBP/mir", "seq", "none"), selected = "RBP/mir", inline = TRUE) %>%
              bs_embed_tooltip("Annotations: 5mer - Ray2013 + Encode, 6mer - Transite R, 7mer TargetScan mir seed")
          ),
          div(
            style = "display: inline-block;vertical-align:top;",
            selectInput("utrlen", NULL, choices = c(200, 500, 1000, "full length"), selected = "full length")
          ),
          div(
            id = "rbptermdiv",
            style = "display: inline-block;vertical-align:top;width:250px",
            textInput("rbpterm", "highlight annotation", value = "MEX3C") %>%
              bs_embed_tooltip("highlights annotation in black, case insensitive", placement = "bottom")
          ),
          div(
            id = "doPlotly2div",
            checkboxInput("doPlotly2", "interactive plot", value = T, width = NULL) %>%
              bs_embed_tooltip("display interactive plot with additional info on hover", placement = "right"),
            style = "width:200px",
          ),
          uiOutput("kmerPlotUI") %>% withLoader(loader = "pacman", proxy.height = paste0(plot_height * 100 / 2, "px"))
        ),
        tabPanel(
          introBox(
            span(icon("circle-notch", class = NULL, lib = "font-awesome"), "Genes_venn",
                 title = "visualize gene overlap between regions by venn diagram, and retrieve lists"
            ),
            data.step = 15,
            data.intro = "Use venn diagram to visualize documented and loaded Genelist/Cart.<br><br>
            Clicking on numbers moves genes of that category to Cart",
            data.position = "top"
          ),
          value = "venn",
          div(
            style = "display: inline-block;vertical-align:top; width: 160px;",
            selectizeInput("seta", "geneset_A",
                           choices = c(
                             "_none", "_Gene_list", "_Cart_list",
                             names(gene_list)
                           ),
                           selected = "fb_sig"
            )
          ),
          div(
            style = "display: inline-block;vertical-align:top; width: 160px;",
            selectizeInput("setb", "geneset_B",
                           choices = c(
                             "_none", "_Gene_list", "_Cart_list",
                             names(gene_list)
                           ),
                           selected = "hy_sig"
            )
          ),
          div(
            style = "display: inline-block;vertical-align:top; width: 160px;",
            selectizeInput("setc", "geneset_C",
                           choices = c(
                             "_none", "_Gene_list", "_Cart_list",
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
            ) %>%
              bs_embed_tooltip("coerce all gene symbols to upper case", placement = "bottom"),
            style = "display: inline-block; width: 100px;"
          ),
          div(
            id = "loadalldiv",
            actionButton("Cart_all", "all genes to Cart") %>%
              bs_embed_tooltip("add all genes from these sets into `Cart` side panel", placement = "bottom"),
            style = "display: inline-block"
          ),
          plotlyOutput("vennPlot") %>% withLoader()
        ),
        tabPanel(
          introBox(
            span(icon("info", class = NULL, lib = "font-awesome"),
                 "About",
                 title = "View version and author info",
            ),
            data.step = 16,
            data.intro = "Additional information on the study and authors.<br><br>
            Hope squirrelBox will be informative for your data explorations.",
            data.position = "top"
          ),
          value = "about",
          uiOutput("intro"),
          uiOutput("track"),
          uiOutput("rawdata"),
          uiOutput("bsgenome"),
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