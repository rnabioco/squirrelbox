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
  rv$listn3 <- ""
  rv$listn4 <- ""
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
  rv$tabinit_as <- 0
  rv$tabinit_enrich <- 0
  rv$tabinit_kmer <- 0
  rv$tabinit_venn <- 0
  rv$starttutorial <- 0
  rv$region_main <- region_main
  rv$region_main2 <- region_main2
  rv$region_order <- region_order
  rv$region_short <- region_short
  rv$region_short_main <- region_short_main
  rv$region_one <- region_one
  
  # hide some checkboxes
  removeModal()
  hide("doKegg")
  hide("doMod")
  hide("doUcsc")
  hide("doEigendiv")
  hide("doTooltips")
  hide("utrlen")
  hide("utr")
  hide("doBr")
  hide("doTis")

  # empty history list to start
  historytab <- c()
  
  # reordering
  observeEvent(input$bsselectconfirm, {
    rv$region_main <- input$bsshow_order
    rv$region_main2 <- input$bshide_order
    rv$region_order <- c(rv$region_main, rv$region_main2)
    print(rv$region_order)
    rv$region_short <- region_short[match(rv$region_order, region_order)]
    rv$region_short_main <- region_short[match(rv$region_main, region_order)]
    rv$region_one <- region_one[match(rv$region_order, region_order)]
  })
  observeEvent(input$bsselectdefault, {
    rv$region_main <- region_main
    rv$region_main2 <- region_main2
    rv$region_order <- region_order
    print(rv$region_order)
    rv$region_short <- region_short
    rv$region_short_main <- region_short_main
    rv$region_one <- region_one
  })
  
  # init
  observeEvent(rv$init, {
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
      jqui_sortable(ui = "#sidediv", operation = "enable", options = list(cancel = ".form-control"))
      jqui_sortable(ui = "#tabMain", operation = "enable")
    }
  })

  # jump to plot
  observeEvent(input$Find, {
    if ((input$geneID != "") & (input$geneID != rv$old)) {
      updateTabsetPanel(session,
        "tabMain",
        selected = "Gene_query"
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

    t1 <- Sys.time()
    if (nrow(plot_temp) == 0) {
      return(ggplot())
    }

    if (input$doTis & input$doBr) {
      mis <- setdiff(rv$region_order, plot_temp$region %>% unique() %>% as.character())
    } else if (!(input$doTis) & input$doBr) {
      mis <- setdiff(rv$region_main, plot_temp$region %>% unique() %>% as.character())
    } else {
      mis <- setdiff(rv$region_main2, plot_temp$region %>% unique() %>% as.character())
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
      plot_temp <- plot_temp %>% filter(!(region %in% rv$region_main2))
    }
    if (!input$doBr) {
      plot_temp <- plot_temp %>% filter(!(region %in% rv$region_main))
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
    plot_temp <- plot_temp %>% mutate(region = factor(region, levels = rv$region_order))
    # plot_temp2 <<- plot_temp

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

    if (input$doPadj == T & nrow(rv$pval) != 0) { #  & input$doPlotly == F
      t2 <- Sys.time()
      padj2 <<- padj
      padj <- padj[str_detect(rownames(padj), paste(rv$region_short_main, collapse = "|")), , drop = FALSE]
      sig_sym <- sig_sym[str_detect(rownames(sig_sym), paste(rv$region_short_main, collapse = "|")), , drop = FALSE]
      temp2 <- calls_sig(padj, sig_sym, as.numeric(input$pval))
      temp2 <- temp2 %>%
        replace_na(list(call1 = list(0))) %>%
        separate(comp, into = c("region", "state1", NA, "state2"), extra = "drop") %>%
        mutate(call1 = as.numeric(call1))
      temp2 <- temp2[, !(names(temp2) %in% c("padj", "call")), drop = F]

      temp2$region <- rv$region_order[factor(temp2$region, level = rv$region_short) %>% as.numeric()]
      temp3 <- groups_to_letters_igraph(temp2) %>%
        mutate(region = factor(region, level = rv$region_order))

      if (verbose_bench) {
        print(paste0("boxPlot1 agg step2: ", Sys.time() - t2))
      }

      if (!(is.na(plot_temp$log2_counts) %>% all())) {
        agg <- aggregate(log2_counts ~ state + region, plot_temp, max)
        agg_min <- aggregate(log2_counts ~ state + region, plot_temp, min)
        agg$min <- agg_min$log2_counts
        agg2 <- agg %>%
          mutate(region = factor(region, level = rv$region_order)) %>%
          group_by(region) %>%
          mutate(
            maxy = max(log2_counts),
            miny = min(min),
            nudgey = (maxy - miny) * 0.1
          ) %>%
          ungroup()
        agg3 <- agg2 %>%
          left_join(temp3[, c("region", "state", "letter"), drop = F], by = c("state", "region")) %>%
          replace_na(list(letter = list("")))

        g <- g +
          geom_text(data = agg3, aes(
            text = letter,
            label = letter,
            y = log2_counts + nudgey,
            x = state,
            group = "text"
          ))
      }
    }

    if (verbose_bench) {
      print(paste0("boxPlot1 step: ", Sys.time() - t1))
    }
    g
  })

  # boxplot size
  boxPlotr <- reactive({
    g <- boxPlot1()
    output$boxPlot <- renderPlot(g)
    if (input$doTis + input$doBr == 2) {
      plotOutput("boxPlot", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100)
    } else {
      plotOutput("boxPlot", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100 / 2)
    }
  })

  # boxplot-plotly
  boxPlotlyr <- reactive({
    g <- boxPlot1()
    output$boxPlot2 <- renderPlotly(ggplotly(g + facet_wrap(~region), tooltip = "text") %>%
      layout(hovermode = "closest"))
    if (input$doTis + input$doBr == 2) {
      plotlyOutput("boxPlot2", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100)
    } else {
      plotlyOutput("boxPlot2", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100 / 2)
    }
  })

  # actually draw boxplot
  output$boxPlotUI <- renderUI({
    if (rv$init != 0 & rv$run2 != 0) {
      if (input$doPlotly == FALSE) {
        suppressWarnings(boxPlotr())
      } else {
        boxPlotlyr()
      }
    }
  })

  # actually draw model module/cluster
  output$EigenPlot <- renderUI({
    if (input$doEigen != T) {
      plotOutput("boxPlot3", height = 1)
    } else {
      plotOutput("boxPlot3", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100 / 2)
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
        title <- ggdraw() +
          draw_label(
            "cluster model expression (empty for genes with insignificant change or filtered out for low expression)",
            fontface = "bold",
            x = 0,
            hjust = 0
          ) +
          theme(plot.margin = margin(0, 0, 0, 7))
        plots <- cowplot::plot_grid(
          plotlist = map(mods, function(x) eigen_gg[[x]]),
          labels = str_c(str_remove(reg, "cluster_"), mods, sep = ": "),
          ncol = 3,
          label_x = .3, hjust = 0
        )
        cowplot::plot_grid(
          title,
          plots,
          ncol = 1,
          rel_heights = c(0.1, 1)
        )
      }
    },
    cacheKeyExpr = {
      tryCatch(rv$mod_df %>% select(-1),
        error = function(e) {
          "error"
        }
      )
    },
    sizePolicy = sizeGrowthRatio(width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100 / 2, growthRate = 1.2)
  )

  # filter data
  outputtab <- reactive({
    inid <- inid()

    t1 <- Sys.time()

    rv$plot_temp <- comb_fil_factor(combined2, combined3, inid)

    temp_orfs <- orfs %>% filter(unique_gene_symbol == rv$plot_temp$unique_gene_symbol[1])
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

    filtered <- rv$plot_temp %>%
      select(any_of(table_cols)) %>%
      unique()

    if (nrow(filtered) == 0) {
      rv$blast <<- ""
      rv$pval <<- data.frame()
      rv$temp_orfs <<- data.frame()
      return(NULL)
    }

    filtered2 <- bed %>%
      filter(gene_id == filtered$gene_id[1])
    filtered2 <- filtered2[, c("chrom", "start", "end", "gene_id"), drop = F]
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

    out <- merge(filtered, filtered2)
    # clusters
    mod1 <- mod %>%
      filter(gene %in% out$unique_gene_symbol)
    if (length(mod1) == 0) {
      rv$mod_df <<- data.frame()
    } else {
      rv$mod_df <<- mod1
    }

    if (verbose_bench) {
      print(paste0("outputtab step: ", Sys.time() - t1))
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
      temp <- rv$temp_orfs %>% select(any_of(orf_cols))
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
      temp <- temp %>% select(any_of(maj_cols))
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
        w <- as.numeric(input$plotw)
        h <- as.numeric(input$ploth)
      } else {
        w <- as.numeric(input$plotw)
        h <- as.numeric(input$ploth) / 2
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
      ungroup() %>% 
      mutate(region = factor(region, levels = rv$region_order))
    
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
    mis <- setdiff(rv$region_main, plot_temp$region %>% unique() %>% as.character())
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
      plot_temp <- plot_temp %>% filter(!(region %in% rv$region_main2))
    }
    if (!input$doBr) {
      plot_temp <- plot_temp %>% filter(!(region %in% rv$region_main))
    }
    plot_temp
  })

  # plotly interative parts of line plot
  observeEvent(linetemp(), {
    rv$line <- rv$line + 1
    shared_d <<- SharedData$new(linetemp, ~unique_gene_symbol, as.character(rv$line))
    shared_d$clearSelection()
  })

  # actually draw line plot
  linePlot1 <- reactive({
    set.seed(1)
    temp <- linetemp()
    if (length(historytablist) == 0) {
      return(ggplotly(ggplot() +
        ggtitle("no genes loaded")))
    }

    # temp2 <<- temp
    if (input$doSummary) {
      temp <- temp %>%
        group_by(region, state) %>%
        summarize(
          log2_counts = log2(mean(2^log2_counts)),
          sem = log2(2^log2_counts) / sqrt(n())
        ) %>%
        mutate(unique_gene_symbol = "summary") %>%
        ungroup()

      g <- ggplot(temp, aes(state, log2_counts,
        group = unique_gene_symbol,
        text = unique_gene_symbol
      )) +
        ylab("log2fold") +
        facet_wrap(~region) +
        geom_point(aes(color = unique_gene_symbol)) +
        geom_line(aes(color = unique_gene_symbol)) +
        geom_errorbar(aes(ymin = log2_counts - sem, ymax = log2_counts + sem), width = .05, size = 0.5)
    } else {
      g <- ggplot(shared_d, aes(state, log2_counts,
        group = unique_gene_symbol,
        text = unique_gene_symbol
      )) +
        ylab("log2fold") +
        facet_wrap(~region) +
        geom_point(aes(color = unique_gene_symbol)) +
        geom_line(aes(color = unique_gene_symbol))
    }
    g
  })

  output$linePlot <- renderPlotly({
    g <- linePlot1() + theme(legend.position = "none")
    fac <- input$doTis + input$doBr
    if ((rv$toolarge == 0) | (rv$toolarge == 1 & rv$go == 2) | input$doSummary) {
      if (input$doName2) {
        ggplotly(g, tooltip = "text", height = as.numeric(input$ploth) * 100 * fac / 2, width = as.numeric(input$plotw) * 100) %>%
          layout(
            autosize = FALSE,
            showlegend = TRUE
          )
      } else {
        ggplotly(g, tooltip = "text", height = as.numeric(input$ploth) * 100 * fac / 2, width = as.numeric(input$plotw) * 100)
      }
    } else {
      if (rv$go <= 0) {
        showModal(modalWarn_line)
      }
      ggplotly(ggplot() +
        ggtitle("plotting cancelled"))
    }
  })

  savePlot5 <- downloadHandler(
    filename = "lineplot.pdf",
    content = function(file) {
      fac <- input$doTis + input$doBr
      g <- linePlot1() + theme(legend.title = element_blank())
      if (!input$doName2) {
        g <- g + theme(legend.position = "none")
      }
      ggplot2::ggsave(file, plot = g, device = "pdf", width = as.numeric(input$plotw), height = as.numeric(input$ploth) / 2 * fac)
    }
  )

  # heatmap
  heatPlot1 <- reactive({
    set.seed(1)
    temp <- linetemp() %>% mutate(region = factor(region, levels = rv$region_order))
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
          column_split = factor(str_remove(colnames(temp2), ":.+"), levels = rv$region_order),
          column_labels = str_remove(colnames(temp2), "^.+:"),
          show_column_names = input$doLabelgene,
          heatmap_legend_param = list(title = "Z-Score")
        )
      } else {
        Heatmap(temp2,
          cluster_rows = input$doRowcluster,
          cluster_columns = input$doColumncluster,
          row_split = factor(str_remove(rownames(temp2), ":.+"), levels = rv$region_order),
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
      plotOutput("heatPlot2", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100 * 2)
    } else {
      plotOutput("heatPlot2", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100)
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
        pdf(file, width = as.numeric(input$plotw), height = as.numeric(input$ploth))
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
    if (input$background == "brain highly expressed") {
      gmtlist1 <- gmtlist_br
    } else if (input$background == "all squirrel genes") {
      gmtlist1 <- gmtlist_sq
    } else if (input$background == "all human genes") {
      gmtlist1 <- gmtlist_human
    }
    tops <- fisher(genevec, gmtlist1[[input$gocat]], gmtlist1[["n"]])
    tops$pathway <- reorder(tops$pathway, tops$minuslog10)
    tops
  })

  richPlot1 <- reactive({
    tops <- richtemp()
    if (nrow(tops) == 0) {
      return(ggplot() +
        ggtitle("no genes loaded"))
    }
    tops <- tops %>%
      dplyr::slice(1:max(min(which(tops$padj > as.numeric(input$pval))), 15)) %>%
      mutate(pathway = str_to_lower(pathway))

    g <- ggplot(
      tops %>% dplyr::slice(1:15),
      aes(x = pathway, y = minuslog10, fill = -minuslog10, text = len)
    ) +
      geom_bar(stat = "identity") +
      xlab(paste0("enriched : ", str_remove(gmt_short, "_"))) +
      scale_x_discrete(limits = rev(tops$pathway[1:15])) +
      coord_flip() +
      cowplot::theme_minimal_vgrid() +
      theme(
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.position = "none"
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      ylab("-log10(FDR)")

    if (input$doPline) {
      g <- g + geom_hline(yintercept = -log10(as.numeric(input$pval)))
    }
    g
  })

  output$richPlot <- renderPlotly({
    p <- ggplotly(richPlot1(),
      source = "richPlot", tooltip = "text", height = as.numeric(input$ploth) * 100, width = as.numeric(input$plotw) * 100
    ) %>%
      layout(autosize = F) %>%
      highlight()
    p
  })

  savePlot2 <- downloadHandler(
    filename = "enriched.pdf",
    content = function(file) {
      ggplot2::ggsave(file, plot = richPlot1(), device = "pdf", width = as.numeric(input$plotw), height = as.numeric(input$ploth))
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

    t1 <- Sys.time()

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
      topsk <<- comp_kmer(
        gene_vec = genevec,
        bac = seqs_precal[[precal]],
        col = utrchoice,
        k = as.numeric(input$km)
      )

      if (verbose_bench) {
        print(paste0("comp_kmer step: ", Sys.time() - t1))
      }

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
      mutate(minuslog10 = -log10(adj_p_value), enrichment = log2(enrichment))
    res <- topsk %>% replace_na(list(RBP = ""))
    res
  })

  kmerPlot1_pre_debounce <- reactive({
    topsk <- kmertemp()
    t1 <- Sys.time()

    de_rbpterm <- input$rbpterm
    if (nrow(topsk) == 0) {
      return(data.frame())
    }
    topsk <- topsk %>%
      mutate(sig = factor(ifelse(minuslog10 >= -log10(as.numeric(input$pval)), "sig", "insig"),
        levels = c("sig", "insig")
      ))
    if (input$kmlab == "RBP/mir") {
      topsk <- topsk %>%
        mutate(text2 = ifelse((sig == "sig" & row_number() <= 15), str_c(kmer, RBP, sep = "\n"), "")) %>%
        mutate(text1 = str_c(kmer, RBP, sep = "\n"))
    } else if (input$kmlab == "none") {
      topsk <- topsk %>% mutate(text2 = "", text1 = "")
    } else {
      topsk <- topsk %>%
        mutate(text2 = ifelse((sig == "sig" & row_number() <= 15), kmer, "")) %>%
        mutate(text1 = kmer)
    }

    if (verbose_bench) {
      print(paste0("kmerplot1 step1: ", Sys.time() - t1))
    }

    topsk
  })

  kmerPlot1_debounce <- reactive({
    de_rbpterm <- input$rbpterm
    # de_rbpterm <- debounce(input$rbpterm, 100)
    topsk <- kmerPlot1_pre_debounce()
    if (nrow(topsk) == 0) {
      return(data.frame())
    }
    topsk <- topsk %>%
      mutate(found = ifelse(str_detect(str_to_upper(RBP), str_to_upper(de_rbpterm)) |
        str_detect(kmer, str_to_upper(de_rbpterm)),
      1,
      0
      ))
    topsk
  })

  kmerPlot1 <- reactive({
    topsk <- kmerPlot1_debounce()
    if (nrow(topsk) == 0) {
      return(ggplot() +
        ggtitle("no genes loaded"))
    }
    g <- ggplot(topsk, aes(
      x = enrichment,
      y = minuslog10,
      text = text1,
      text2 = RBP,
      label = text2
    )) +
      geom_point(aes(color = sig), size = 0.5) +
      scale_color_manual(values = c("#FC8D62", "#B3B3B3")) +
      geom_point(data = topsk %>%
        filter(found == 1), color = "black", size = 0.5) +
      ggrepel::geom_text_repel(box.padding = 0.05, size = 3, aes(label = text2)) +
      xlab("log2enrichment") +
      ylab("-log10(FDR)") +
      labs(color = "") +
      scale_y_continuous(expand = c(0, 0))

    if (input$doPline2) {
      g <- g + geom_hline(yintercept = -log10(as.numeric(input$pval)))
    }

    if (verbose_bench) {
      print(paste0("kmerplot1 step2: ", Sys.time() - t1))
    }

    g
  })

  kmerPlotr <- reactive({
    output$kmerPlot <- renderPlot(kmerPlot1())
    plotOutput("kmerPlot", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100)
  })

  kmerPlotlyr <- reactive({
    g <- kmerPlot1()
    g2 <- suppressWarnings(ggplotly(g,
      source = "kmerPlotly",
      tooltip = "text", height = as.numeric(input$ploth) * 100, width = as.numeric(input$plotw) * 100
    )) %>%
      layout(autosize = F)
    output$kmerPlot2 <- renderPlotly(g2)
    plotlyOutput("kmerPlot2", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100)
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
      ggplot2::ggsave(file, plot = kmerPlot1(), device = "pdf", width = as.numeric(input$plotw), height = as.numeric(input$ploth))
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
    load_vec <- c(input$seta, input$setb, input$setc) == "_Gene_list"
    if (sum(load_vec) > 0) {
      for (element in names(temp)[load_vec]) {
        temp[[element]] <- historytablist
      }
    }
    load_vec <- c(input$seta, input$setb, input$setc) == "_Cart_list"
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
    non_none <- !sapply(a, is.null) & !duplicated(c(input$seta, input$setb, input$setc)) & !(c(input$seta, input$setb, input$setc) == "_none")
    if (sum(non_none) > 1) {
      g <- ggvenn(a, names(a)[non_none], show_elements = TRUE)
      g2 <- ggvenn(a, names(a)[non_none],
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
      g <- ggvenn(a, c(names(a)[non_none], NA), show_elements = TRUE)
      g2 <- ggvenn(a, c(names(a)[non_none], NA),
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
        ggtitle("no Genelist loaded")
    }
  })

  output$vennPlot <- renderPlotly({
    p <- ggplotly(vennPlot1(),
      source = "vennPlot",
      tooltip = "text2", height = as.numeric(input$ploth) * 100, width = as.numeric(input$plotw) * 100
    ) %>%
      layout(autosize = F, showlegend = FALSE) %>%
      highlight()
    p
  })

  savePlot6 <- downloadHandler(
    filename = "venn.pdf",
    content = function(file) {
      ggplot2::ggsave(file, plot = vennPlot1(), device = "pdf", width = as.numeric(input$plotw), height = as.numeric(input$ploth))
    }
  )

  observeEvent(event_data("plotly_click", source = "vennPlot"), {
    aa <- event_data("plotly_click", source = "vennPlot")
    if (!is.null(aa$pointNumber)) {
      gene_string <- rv$venntext[aa$pointNumber + 1]
      gene_vec <- tryCatch(str_split(gene_string, ",")[[1]],
        error = function(e) {
          ""
        }
      )
      carttablist <<- gene_vec
      rv$listn2 <- length(carttablist)
      updateTabsetPanel(session,
        "side2",
        selected = "cart"
      )
      rv$listn2renew <- rv$listn2renew + 1
    }
  })

  # circos
  circostrack <- reactive({
    rv$line_refresh
    if (input$guse == "_none") {
      rv$bed <<- ""
      rv$listn3 <<- "Displaying 0 genes"
      return(NA)
    } else if (input$guse == "_Gene_list") {
      templist <- historytablist
    } else if (input$guse == "_Cart_list") {
      templist <- carttablist
    } else {
      templist <- gene_list[[input$guse]]
    }
    if (length(templist) == 0) {
      rv$bed <<- ""
      rv$listn3 <<- "Displaying 0 genes"
      return(NA)
    }

    bed_f <- bed %>% filter((unique_gene_symbol %in% templist) | (str_to_upper(unique_gene_symbol) %in% templist))

    listn3_1 <- length(templist)
    listn3_2 <- bed_f %>%
      filter(as.numeric(str_remove(chrom, "Itri")) <= chrlimit) %>%
      pull(unique_gene_symbol) %>%
      unique() %>%
      length()
    rv$listn3 <<- paste0("Displaying ", listn3_2, " out of ", listn3_1, " genes")

    bed_f <- bed_f %>% filter(as.numeric(str_remove(chrom, "Itri")) <= chrlimit)
    bed_f <- bed_f %>%
      group_by(unique_gene_symbol) %>%
      bed_merge(max_dist = 10000000)
    rv$bed <<- bed_f
    bed_f <- bed_f %>% bed_slop(sq_g, both = 2500000, trim = TRUE)

    arcs_chromosomes <- str_remove(bed_f$chrom, "Itri") # Chromosomes on which the arcs should be displayed
    arcs_begin <- bed_f$start
    arcs_end <- bed_f$end
    arcs_lab <- bed_f$unique_gene_symbol

    tracklist <- BioCircosArcTrack("myArcTrack", arcs_chromosomes, arcs_begin, arcs_end,
      labels = arcs_lab,
      minRadius = 0.45, maxRadius = 0.65, colors = "blue", opacities = 0.4
    )
    tracklist
  })

  circostrack2 <- reactive({
    rv$line_refresh
    if (input$guse2 == "_none") {
      rv$listn4 <<- "Displaying 0 genes"
      return(NA)
    } else if (input$guse2 == "_Gene_list") {
      templist <- historytablist
    } else if (input$guse2 == "_Cart_list") {
      templist <- carttablist
    } else {
      templist <- gene_list[[input$guse2]]
    }
    if (length(templist) == 0) {
      rv$listn4 <<- "Displaying 0 genes"
      return(NA)
    }

    bed_f <- bed %>% filter((unique_gene_symbol %in% templist) | (str_to_upper(unique_gene_symbol) %in% templist))

    listn4_1 <- length(templist)
    listn4_2 <- bed_f %>%
      filter(as.numeric(str_remove(chrom, "Itri")) <= chrlimit) %>%
      pull(unique_gene_symbol) %>%
      unique() %>%
      length()
    rv$listn4 <<- paste0("Displaying ", listn4_2, " out of ", listn4_1, " genes")

    bed_f <- bed_f %>% filter(as.numeric(str_remove(chrom, "Itri")) <= chrlimit)
    bed_f <- bed_f %>%
      group_by(unique_gene_symbol) %>%
      bed_merge(max_dist = 1000000000)
    rv$bed <<- bed_f
    bed_f <- bed_f %>% bed_slop(sq_g, both = 2500000, trim = TRUE)

    arcs_chromosomes <- str_remove(bed_f$chrom, "Itri") # Chromosomes on which the arcs should be displayed
    arcs_begin <- bed_f$start
    arcs_end <- bed_f$end
    arcs_lab <- bed_f$unique_gene_symbol

    tracklist <- BioCircosArcTrack("myArcTrack", arcs_chromosomes, arcs_begin, arcs_end,
      labels = arcs_lab,
      minRadius = 0.2, maxRadius = 0.4, colors = "green", opacities = 0.4
    )
    tracklist
  })

  circostrack3 <- reactive({
    if (input$guse3 == "_none") {
      return(NA)
    } else if (input$guse3 == "Editing_Riemondy2018") {
      tracklist <- BioCircosBarTrack(
        trackname = "fc",
        chromosomes = str_remove(edits$chrom, "Itri"),
        starts = edits$start,
        ends = edits$end,
        values = edits$max,
        labels = edits$unique_gene_symbol,
        range = c(0, 1),
        color = "#000000",
        maxRadius = 0.87,
        minRadius = 0.7
      )
      return(tracklist)
    } else {
      bed_temp <- bed_fc %>%
        filter(region == input$guse3) %>%
        filter(as.numeric(str_remove(chrom, "Itri")) <= chrlimit) %>%
        filter(abs(fold) >= 0.1)
    }
    temp3 <- bed_temp
    tracklist <- BioCircosBarTrack(
      trackname = "fc",
      chromosomes = str_remove(temp3$chrom, "Itri"),
      starts = temp3$start,
      ends = temp3$end,
      values = temp3$fold,
      labels = temp3$unique_gene_symbol,
      range = c(0, log2(ceiling(max(bed_fc$fold)))),
      color = "#000000",
      maxRadius = 0.87,
      minRadius = 0.7
    )
    tracklist
  })

  circosPlot1 <- reactive({
    tracks <- BioCircosTracklist()
    if (!(is.na(circostrack()))) {
      tracks <- tracks + circostrack()
    }
    if (!(is.na(circostrack2()))) {
      tracks <- tracks + circostrack2()
    }
    if (!(is.na(circostrack3()))) {
      tracks <- tracks + circostrack3()
    }
    if (input$guse3 == "Editing_Riemondy2018") {
      BioCircos(tracks,
        genome = sq1 %>% setNames(names(sq1) %>% str_remove("Itri")),
        ARCMouseOverTooltipsHtml02 = "<!––",
        ARCMouseOverTooltipsHtml03 = "",
        ARCMouseOverTooltipsHtml04 = "––><br/>",
        chrPad = 0.023,
        genomeLabelTextSize = "9pt",
        BARMouseOverTooltipsHtml02 = "<!––",
        BARMouseOverTooltipsHtml03 = "",
        BARMouseOverTooltipsHtml04 = "––><br/>",
        BARMouseOverTooltipsHtml05 = "<br/>max_fraction: "
      )
    } else {
      BioCircos(tracks,
        genome = sq1 %>% setNames(names(sq1) %>% str_remove("Itri")),
        ARCMouseOverTooltipsHtml02 = "<!––",
        ARCMouseOverTooltipsHtml03 = "",
        ARCMouseOverTooltipsHtml04 = "––><br/>",
        chrPad = 0.023,
        genomeLabelTextSize = "9pt",
        BARMouseOverTooltipsHtml02 = "<!––",
        BARMouseOverTooltipsHtml03 = "",
        BARMouseOverTooltipsHtml04 = "––><br/>",
        BARMouseOverTooltipsHtml05 = "<br/>max_log2FC: "
      )
    }
  })

  circosPlotr <- reactive({
    output$circos <- renderBioCircos({
      circosPlot1()
    })
    BioCircosOutput("circos", width = as.numeric(input$plotw) * 100, height = as.numeric(input$ploth) * 100)
  })

  output$circosUI <- renderUI({
    circosPlotr()
  })

  savePlot7 <- downloadHandler(
    filename = "genome.html",
    content = function(file) {
      htmlwidgets::saveWidget(circosPlot1(), file = file)
    }
  )

  saveFilteredbed <- downloadHandler("filtré.csv", content = function(file) {
    write_csv(rv$bed, file)
  })

  output$listn3 <- renderUI({
    HTML(str_c("<strong><h5>", rv$listn3, "</h5></strong>"))
  })

  output$listn4 <- renderUI({
    HTML(str_c("<strong><h5 style='color:green'>", rv$listn4, "</h5></strong>"))
  })

  # misc

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
  observeEvent(input$geneID, {
    if (rv$run2 == 1 & input$geneID != "" & !is.null(input$geneID) & input$tabMain == "Gene_query") {
      rv$run2 <- 1
      shinyjs::click("Find")
    }
  })

  observeEvent(input$geneID2, {
    rv$run2 <- 1
    if (input$geneID2 %in% autocomplete_list) {
      updateSelectizeInput(session,
        inputId = "geneID",
        selected = input$geneID2,
        choices = autocomplete_list,
        server = T
      )
    }
  })

  observeEvent(input$geneID3, {
    rv$run2 <- 1
    if (input$geneID3 %in% autocomplete_list) {
      updateSelectizeInput(session,
        inputId = "geneID",
        selected = input$geneID3,
        choices = autocomplete_list,
        server = T
      )
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
    historytablist <- historytablist %>% unique()
    rv$line_refresh <- rv$line_refresh + 1
  })

  # cart list
  onclick("Add", {
    carttablist <- unique(c(historytab[1], carttablist))
    rv$listn2 <- length(carttablist)
  })

  onclick("Empty", {
    carttablist <- c()
    rv$listn2 <- length(carttablist)
  })


  output$listn2 <- renderUI({
    HTML(str_c("<strong><h5> # in Cart: ", rv$listn2, "</h5></strong>"))
  })

  output$saveList <- downloadHandler("cart.txt", content = function(file) {
    write_lines(carttablist, file)
  })

  onclick("Load", {
    historytablist <- carttablist
    historytablist <- historytablist %>% unique()
    rv$line_refresh <- rv$line_refresh + 1
    updateTabsetPanel(session,
      "side2",
      selected = "load"
    )
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
  output$genes <- DT::renderDataTable({
    temp <- orftbl() %>%
      select(
        unique_gene_symbol,
        contains("cluster"),
        contains("LRT"),
        everything()
      )
    if (!(file.exists(paste0(annotpath, "/genes.csv")))) {
      write_lines(colnames(temp), paste0(annotpath, "/genes.csv"))
    }
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
      callback = JS(paste0("var tips = [", genes_tips, "],
                            firstRow = $('#genes thead tr th');
                            for (var i = 0; i < tips.length; i++) {
                              $(firstRow[i]).attr('title', tips[i]);
                            }"))
    ) %>%
      DT::formatRound(columns = padjcol, digits = 4)
  })

  saveFiltered <- downloadHandler("filtré.csv", content = function(file) {
    s <- input$genes_rows_all
    write_csv((orftbl() %>%
      select(
        unique_gene_symbol,
        contains("cluster"),
        contains("LRT"),
        everything()
      ) %>% select(unique_gene_symbol, everything()))[s, ], file)
  })

  onclick("loadtab", {
    s <- input$genes_rows_all
    if (length(s) == 0) {
      return()
    }
    historytablist <- orftbl()[s, ] %>% pull(unique_gene_symbol)
    historytablist <- historytablist %>% unique()
    rv$line_refresh <- rv$line_refresh + 1
    updateTabsetPanel(session,
      "side2",
      selected = "load"
    )
  })

  observeEvent(input$genes_rows_selected, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = orftbl()[["unique_gene_symbol"]][input$genes_rows_selected],
      choices = autocomplete_list,
      server = T
    )
  })

  # explore majiq table
  majtbl <- reactive({
    if (input$doJoin) {
      maj %>% left_join(fulltbl_collapse, by = c("gene_id", "unique_gene_symbol"))
    } else {
      maj
    }
  })

  output$alt <- DT::renderDataTable({
    temp <- majtbl() %>%
      select(
        unique_gene_symbol,
        contains("significant"),
        LSV_ID, A5SS, A3SS, ES,
        everything()
      )
    if (!(file.exists(paste0(annotpath, "/alt.csv")))) {
      write_lines(colnames(temp), paste0(annotpath, "/alt.csv"))
    }
    DT::datatable(
      temp,
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
      ),
      callback = JS(paste0("var tips = [", maj_tips, "],
                            firstRow = $('#alt thead tr th');
                            for (var i = 0; i < tips.length; i++) {
                              $(firstRow[i]).attr('title', tips[i]);
                            }"))
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
    if (length(s) == 0) {
      return()
    }
    historytablist <- majtbl()[s, ] %>%
      pull(unique_gene_symbol)
    historytablist <- historytablist %>% unique()
    rv$line_refresh <- rv$line_refresh + 1
    updateTabsetPanel(session,
      "side2",
      selected = "load"
    )
  })

  observeEvent(input$alt_rows_selected, {
    rv$run2 <- 1
    updateSelectizeInput(session,
      inputId = "geneID",
      selected = majtbl()[["unique_gene_symbol"]][input$alt_rows_selected],
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

  output$bsgenome <- renderUI({
    url <- str_c(
      "https://bioconductor.org/packages/release/data/annotation/html/", bsgenomeL
    )
    clean <- a(bsgenomeL,
      href = url
    )
    tagList(tags$h6("Full genome sequences are deposited to ___ and stored in Biostring/BSgenome format, ", clean))
  })

  output$GOversion <- renderUI({
    url <- "https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp"
    clean <- a(gmt_file,
      href = url
    )

    tagList(tags$h6("GO database version: ", clean))
  })

  output$version <- renderUI({
    clean <- a(versionN,
      href = url
    )
    s3link <- a("S3",
      href = s3
    )
    docklink <- a("docker_hub",
      href = docker
    )
    tagList(tags$h6(icon("github-square"), "squirrelBox version: ", clean, "; or available to download and run locally: ", s3link," and ", docklink))
  })

  output$explain <- DT::renderDataTable({
    dfreg <- data.frame(
      region = region_order,
      short = region_short,
      letter = region_one
    ) %>% filter(region %in% (combined2$region %>% unique()))
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

  output$explain3 <- DT::renderDataTable({
    dfreg3 <- data.frame(
      file = c(
        "brain_editing_site_proportions.bed",
        "clusters.feather",
        "combined2.feather",
        "combined3.feather",
        "fc.rds",
        "MAJIQ_dpsi_summary_sig_squirrelBox.tsv.gz",
        "MiPepid_pred.csv",
        "novel_domains.csv",
        "padj_orf.feather",
        "seqs_precal_noG.rds",
        "SmProt_blast.csv",
        "utrs_sq_noG.feather",
        "utrs_sq.feather"
      ),
      desc = c(
        "RNA editing events found from seq data (Riemondy2018)",
        "cluster assignments to each tissue",
        "log expression values per gene per animal",
        "gene id, symbol, version info",
        "gene interval and max fold change",
        "alternative splicing summary",
        "micropeptide prediction",
        "novel domain scanning",
        "scanned orfs and pair-wise differential expression p-values",
        "precalculated kmer info of all sequences (non-novel)",
        "homologous sequences to SmProt micropeptide database",
        "all utr5/3/cds sequences (non-novel)",
        "all utr5/3/cds sequences"
      )
    )
    DT::datatable(dfreg3,
      escape = FALSE,
      selection = "none",
      rownames = FALSE,
      options = list(
        searchable = FALSE,
        dom = "t",
        paging = FALSE
      )
    )
  })

  # graying out buttons, warnings/hints
  observe({
    if (input$tabMain == "Gene_query") {
      enable("savePlot")
      enable("savePlot2")
      disable("saveTable")
      disable("saveTable2")
      output$savePlot <- savePlot
      output$savePlot2 <- savePlot
    } else if (input$tabMain == "table_data") {
      disable("savePlot")
      disable("savePlot2")
      enable("saveTable")
      enable("saveTable2")
      output$saveTable <- saveFiltered
      output$saveTable2 <- saveFiltered
      if (rv$tabinit_data == 0) {
        showNotification("Type `low...high` to input custom range for numeric filtering on column, and upper or lower limits can be omitted. For instance, ...0.001 filters for significant p value, and 200...2000 filters RNA length to that range.",
          type = "message"
        )
        rv$tabinit_data <<- 1
      }
    } else if (input$tabMain == "table_AS") {
      disable("savePlot")
      disable("savePlot2")
      enable("saveTable")
      enable("saveTable2")
      output$saveTable <- saveFilteredAS
      output$saveTable2 <- saveFilteredAS
      if (rv$tabinit_as == 0) {
        showNotification("Type `low...high` to input custom range for numeric filtering on column, and upper or lower limits can be omitted. For instance, ...0.001 filters for significant p value, and 200...2000 filters RNA length to that range.",
          type = "message"
        )
        rv$tabinit_as <<- 1
      }
    } else if (input$tabMain == "line_plot") {
      enable("savePlot")
      enable("savePlot2")
      disable("saveTable")
      disable("saveTable2")
      output$savePlot <- savePlot5
      output$savePlot2 <- savePlot5
    } else if (input$tabMain == "enrichment_plot") {
      enable("savePlot")
      enable("savePlot2")
      enable("saveTable")
      enable("saveTable2")
      output$savePlot <- savePlot2
      output$savePlot2 <- savePlot2
      output$saveTable <- saveEnrich
      output$saveTable2 <- saveEnrich
    } else if (input$tabMain == "heat_plot") {
      enable("savePlot")
      enable("savePlot2")
      disable("saveTable")
      disable("saveTable2")
      output$savePlot <- savePlot3
      output$savePlot2 <- savePlot3
    } else if (input$tabMain == "kmer_analysis") {
      enable("savePlot")
      enable("savePlot2")
      enable("saveTable")
      enable("saveTable2")
      output$savePlot <- savePlot4
      output$savePlot2 <- savePlot4
      output$saveTable <- saveK
      output$saveTable2 <- saveK
    } else if (input$tabMain == "venn") {
      enable("savePlot")
      enable("savePlot2")
      disable("saveTable")
      disable("saveTable2")
      output$savePlot <- savePlot6
    } else if (input$tabMain == "about") {
      disable("savePlot")
      disable("savePlot2")
      disable("saveTable")
      disable("saveTable2")
    } else if (input$tabMain == "genome") {
      enable("savePlot")
      enable("savePlot2")
      enable("saveTable")
      enable("saveTable2")
      output$savePlot <- savePlot7
      output$saveTable <- saveFilteredbed
      output$saveTable2 <- saveFilteredbed
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

  onclick("bssum", {
    updateCheckboxInput(session = session, inputId = "doSummary", value = T, label = "summary line")
    rv$go <- 1
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

  modalWarn_line <- draggableModalDialog(
    id = "bsconfirm",
    icon("exclamation-triangle"),
    p("plotting large number of genes may be\n slow and hard to interpret"),
    br(),
    footer = NULL,
    size = "s",
    easyClose = FALSE,
    fade = TRUE,
    actionButton("bsgo", "Go"),
    actionButton("bssum", "Summary only"),
    actionButton("bscancel", "Cancel")
  )

  observeEvent(input$tutorial, {
    introjs(session,
      options = list(
        "nextLabel" = ">",
        "prevLabel" = "<",
        "skipLabel" = "skip",
        "overlayOpacity" = -1
      ) # ,events = list("onexit" = I("alert('abc')"))
    )
  })

  observeEvent(rv$run2, {
    if (start_tutorial & rv$starttutorial == 0) {
      rv$starttutorial <- 1
      introjs(session, options = list(
        "nextLabel" = ">",
        "prevLabel" = "<",
        "skipLabel" = "skip",
        "overlayOpacity" = -1
      )) # events = list("onexit" = I("document.getElementsByClassName('introjs-nextbutton').blur()")))
    }
  })
}
