#' Server
#'
#' Query strings are supported as follows:
#' ?genes=[list of genes]&taxa=[human|mouse|rat]&...
server <- function(input, output, session) {
  # On connect, observe the query string. Search if any genes are specified
  observeEvent(input$LOAD, {
    output$results_header <- renderUI(generateResultsHeader(
      HTML('<div style="margin-bottom: 10px"><h2 style="display: inline">No enrichments yet</h2></div>')))
    
    query <- getQueryString(session)
    options <- getConfig()
    
    genes <- query$genes %>% jsonify
    
    sig <- switch(is.null(query$sig) + 1, jsonify(query$sig), options$sig$value)
    fc <- switch(is.null(query$fc) + 1, jsonify(query$fc) %>% as.numeric, options$fc$value)
    
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    session$sendCustomMessage('querySet', genes)
    
    updateTextInput(session, 'sig', value = sig)
    updateSliderInput(session, 'fc', value = fc)
    
    for(mName in Filter(function(x) !(x %in% c('sig', 'fc')), names(options))) {
      do.update(session, options[[mName]], ifelse(is.null(query[[mName]]), options[[mName]]$value, query[[mName]]))
    }

    # Search if any genes are specified. Specify parameters to make sure no reactives.
    # Run it next ms to give the server time to propagate layout changes.
    if(!is.null(query$genes))
      shinyjs::delay(1, searchGenes(genes, sig, query$taxa, update = F))
  }, once = T, ignoreNULL = T)
  
  observeEvent(input$genes.csv, {
    if(!is.null(input$genes.csv))
      session$sendCustomMessage('fileUpload', T)
  })
  
  output$genes.csv.ui <- renderUI({
    input$reset
    fileInput('genes.csv', 'Or Upload CSV', accept = 'text/csv', placeholder = 'N/A')
  })
  
  # Reset the search bar
  observeEvent(input$reset, {
    updateTextInput(session, 'genes', value = '')
    session$sendCustomMessage('fileUpload', F)
    updateSelectizeInput(session, 'taxa', selected = NULL)
    updateSelectizeInput(session, 'method', selected = NULL)
    updateTextInput(session, 'sig', value = NULL)
    
    options <- getConfig()
    for(mName in Filter(function(x) !(x %in% c('taxa', 'method', 'sig')), names(options))) {
      do.update(session, options[[mName]], options[[mName]]$value)
    }
  })
  
  # Search
  observeEvent(input$search, {
    genes <- input$genes
    if(!is.null(input$genes.csv))
      genes <- read.csv(input$genes.csv$datapath, header = F)$V1 %>% as.character
    
    searchGenes(genes)
  })
  
  # Open other tabs
  observeEvent(input$tabs, {
    if(input$tabs == 'Gene Info' && !is.null(session$userData$plotData) && is.null(session$userData$genesRendered)) {
      session$userData$genesRendered <- T
      
      synchronise({
        geneEvidence(session$userData$plotData$gene.ID,
                     session$userData$plotData$options$taxa$value)$
          then(function(evidence) {
            output$results_genes <- generateGenePage(evidence)
          })
      })
    } else if(input$tabs == 'GO Enrichment' && !is.null(session$userData$plotData) && is.null(session$userData$goRendered)) {
      session$userData$goRendered <- T
      
      output$results_go <- generateGOPage(ora(hitlist = session$userData$plotData$gene.Name,
                                              annotation = paste0('Generic_', session$userData$plotData$options$taxa$value)))
    } else if(input$tabs == 'Hierarchical View (beta)' && !is.null(session$userData$plotData) && is.null(session$userData$treeRendered)) {
      session$userData$treeRendered <- T
      
      output$results_tree <- generateResultsTree(session$userData$plotData$conditions, session$userData$plotData$options)
    } else if(input$tabs == 'Word Cloud' && !is.null(session$userData$plotData) && is.null(session$userData$cloudRendered)) {
      session$userData$cloudRendered <- T
      
      output$results_cloud <- generateResultsCloud(session$userData$plotData$conditions, session$userData$options)
    }
  })
  
  # Force an updated of the gene view when we get new search results
  observeEvent(input$UPDATED, {
    session$userData$genesRendered <- NULL
    session$userData$goRendered <- NULL
    session$userData$treeRendered <- NULL
    session$userData$cloudRendered <- NULL

    if(input$tabs == 'Gene Info') {
      shinyjs::delay(100, {
        if(input$tabs == 'Gene Info' && !is.null(session$userData$plotData) && is.null(session$userData$genesRendered)) {
          session$userData$genesRendered <- T
          
          synchronise({
            geneEvidence(session$userData$plotData$gene.ID, # TODO for "any"
                         session$userData$plotData$options$taxa$value)$
              then(function(evidence) {
                output$results_genes <- generateGenePage(evidence)
              })
          })
        }
      })
    } else if(input$tabs == 'GO Enrichment') {
      session$userData$goRendered <- T
      
      output$results_go <- generateGOPage(ora(hitlist = session$userData$plotData$gene.Name, # TODO for "any"
                                              annotation = paste0('Generic_', session$userData$plotData$options$taxa$value)))
    } else if(input$tabs == 'Hierarchical View (beta)') {
      session$userData$treeRendered <- T
      
      output$results_tree <- generateResultsTree(session$userData$plotData$conditions, session$userData$plotData$options)
    } else if(input$tabs == 'Word Cloud' && !is.null(session$userData$plotData) && is.null(session$userData$cloudRendered)) {
      session$userData$cloudRendered <- T
      
      output$results_cloud <- generateResultsCloud(session$userData$plotData$conditions, session$userData$options)
    }
  })
  
  observeEvent(input$RANDOM_GENES, {
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    while(length(genes <- DATA.HOLDER[[input$taxa]]@gene.meta[sample(1:.N, sample(1:10, 1))] %>%
      as.data.frame %>% .[, sample(c(1, 3, 4, 6), 1)] %>%
      .[. != 'None']) == 0) {
    }
    
    genes <- genes %>% {
      paste0('[', paste0(., collapse = ','), ']')
    } %>% jsonify
    
    session$sendCustomMessage('queryReset', genes)
  })
  
  # Open plot
  observeEvent(input$plotData, {
    shiny::showModal(modalDialog(
      easyClose = T,
      size = 'l',
      footer = NULL,
      fluidRow(style = 'height: 85vh;',
               column(10, plotlyOutput('plot') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase'))))),
               column(2, style = 'padding-top: 15px; padding-right: 30px;',
                      fluidRow(style = 'display: flex; flex-direction: row; justify-content: space-evenly; margin-bottom: 15px;',
                               downloadButton('plot_save_jpg', 'Save JPG'), downloadButton('plot_save_pdf', 'Save PDF')),
                      selectInput('plot_type', 'Type', list('Boxplot', 'Violin plot', 'Heatmap', 'Scatterplot', 'Jitterplot')),
                      selectInput('plot_data', 'Data', list('Gene Expression')),
                      pickerInput('plot_genes', 'Genes', list('Loading...'), multiple = T),
                      pickerInput('plot_conditions', 'Contrasts', list('Loading...'), multiple = T))
      )
    ), session)
    
    generatePlot <- function(value) {
      session$userData$plotData$exprInfo <- value
      
      updatePickerInput(session, 'plot_genes',
                        choices = as.list(intersect(session$userData$plotData$gene.Name, rownames(value$expr))),
                        selected = intersect(session$userData$plotData$gene.Name, rownames(value$expr)))
      updatePickerInput(session, 'plot_conditions', # TODO Selection strategy
                        choices = as.list(session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)]),
                        selected = data.table::first(session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)]))
      
      output$plot <- renderPlotly({
        generateResultsPlot(session$userData$plotData$gene.Name,
                            session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)],
                            value,
                            session$userData$plotData$options,
                            input$plot_genes, input$plot_conditions, input$plot_type, input$plot_data)
      })
    }
    
    if(is.null(session$userData$plotData$exprInfo)) {
      # If we didn't already process we need to give the UI a tick to update and then fetch the data.
      if(!is.null(session$userData$plotData$conditions)) {
        shinyjs::delay(1, { 
          synchronise(getTags(session$userData$plotData$options$taxa$value,
                              session$userData$plotData$experiments) %>% # TODO Selection strategy
                        .[cf.Cat %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.Cat)]) &
                            cf.BaseLongUri %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.BaseLongUri)]) &
                            cf.ValLongUri %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.ValLongUri)])] %>%
                        unique %>% # TODO for "any"
                        merge(DATA.HOLDER[[session$userData$plotData$options$taxa$value]]@experiment.meta[rsc.ID %in% session$userData$plotData$experiments, .(rsc.ID, ee.ID, ee.NumSamples)],
                              by = 'rsc.ID', sort = F, allow.cartesian = T) %>% {
                          if(nrow(.) == 0) return(NULL)
                          
                          # Select some samples to queue for gene expression visualization
                          while((tmp <- .[sample(1:nrow(.))])[1, ee.NumSamples] > getOption('max.gemma')) {
                          }
                          
                          tmp[cumsum(ee.NumSamples) <= getOption('max.gemma')]
                        } %>% {
                          geneExpression(.[, unique(ee.ID)],
                                         .[, unique(rsc.ID)], # TODO for "any"
                                         session$userData$plotData$options$taxa$value,
                                         session$userData$plotData$gene.ID)$then(generatePlot)
                        })
        })
      }
    } else {
      # If we already processed then we can just display it
      output$plot <- renderPlotly({
        generateResultsPlot(session$userData$plotData$gene.Name,
                            session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)],
                            session$userData$plotData$exprInfo,
                            session$userData$plotData$options,
                            input$plot_genes, input$plot_conditions, input$plot_type, input$plot_data)
      })
    }
  })
  
  #' Display a failure message
  #'
  #' Ends the search protocol with a failure message
  endFailure <- function() {
    setProgress(environment())
    output$results_header <- renderUI({ # TODO Tooltips aren't beautiful
      generateResultsHeader(HTML('<h2 data-toggle="tooltip" data-placement="top" title="Relax thresholds or modify gene set.">Invalid search.</h2>'))
    })
    output$results <- NULL
  }
  
  #' Display an empty message
  #'
  #' Ends the search protocol with no results
  endEmpty <- function() {
    setProgress(environment())
    output$results_header <- renderUI({ generateResultsHeader('No conditions found in scope.') })
    output$results <- NULL
  }
  
  #' Display a success message
  #'
  #' Ends the search protocol successfully
  #' 
  #' @param genes The genes that were searched
  #' @param experiments The experiment rankings that were obtained
  #' @param conditions The condition rankings that were obtained
  #' @param options Any additional options that were passed
  endSuccess <- function(genes, experiments, conditions, options) {
    # Generate the results header
    if(exists('output'))
      output$results_header <- renderUI({
        generateResultsHeader(HTML(paste0('<div style="margin-bottom: 10px"><h2 style="display: inline">Enriched ',
                                          nrow(experiments), ' experiment', ifelse(nrow(experiments) > 1, 's', ''), ' for ',
                                          ifelse(which(colnames(experiments) == 'score') == 2, colnames(experiments)[1],
                                                 paste0('<span data-toggle="tooltip" data-placement="top" title="',
                                                        paste0(colnames(experiments)[1:(which(colnames(experiments) == 'score') - 1)], collapse = ', '), '">',
                                                        which(colnames(experiments) == 'score') - 1, ' gene', ifelse(which(colnames(experiments) == 'score') > 2, 's', ''), '</span>')),
                                          '</h2><span class="timestamp">in ',
                                          format(difftime(session$userData$endTime, session$userData$startTime), digits = 3), '.</span></span>')))
      })
    
    conditions[, Contrast := paste0('<b>', cf.BaseLongUri, '</b> vs. <b>', cf.ValLongUri, '</b>')]
    
    advanceProgress('Cross linking')
    
    conditions <- conditions %>%
      .[abs(stat) >= 1] %>%
      .[, stat := round(stat, 3)] %>%
      setorder(distance) %>%
      .[, head(.SD, 1), stat]
    
    tmp <- getTags(options$taxa$value, experiments$rn) %>%
      .[cf.Cat %in% conditions[, unique(cf.Cat)] &
          cf.BaseLongUri %in% conditions[, unique(cf.BaseLongUri)] &
          cf.ValLongUri %in% conditions[, unique(cf.ValLongUri)]] %>% # TODO for "any"
      merge(DATA.HOLDER[[options$taxa$value]]@experiment.meta[rsc.ID %in% experiments$rn, .(rsc.ID, ee.ID, ee.Name)],
            by = 'rsc.ID', sort = F, allow.cartesian = T) %>%
      .[, N := length(unique(ee.ID)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
    
    # This chunk is costly, but associates all Gemma linkages. It would be great if we could defer this
    # but I don't know how to :(
    # It's much faster to not query ee.Name though
    if(options$gemmaLink$value) {
      tmp[, Evidence :=
          paste0('<span data-toggle="popover" title="Experiments" data-html="true" data-content="',
                 paste0('<a target=_blank href=https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=', unique(ee.ID), '>', unique(ee.ID), '</a>') %>% paste0(collapse = ', '),
                 '">',
                 paste(N, paste0('Experiment', ifelse(N > 1, 's', '')), '<i class="fas fa-question-circle" style="cursor: pointer;"></i>'),
                 '</span>'),
        .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
    } else
      tmp[, Evidence := paste(N, paste0('Experiment', ifelse(N > 1, 's', ''))), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
    
    tmp %>%
      .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, N, Evidence)] %>%
      unique %>%
      merge(conditions, by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), sort = F) %>%
      setorder(-stat, na.last = T) %>%
      setnames(c('stat', 'C', 'distance'), c('Test Statistic', 'Observations', 'Ontology Steps')) %>%
      .[, Observations := round(Observations, 2)]
  }
  
  #' Handle Search
  #' 
  #' Call the processing function to handle searching and update the display.
  #'
  #' @param genes A character vector of entrez IDs
  #' @param options The search options
  handleSearch <- function(genes, options) {
    experiments <- search(genes, options, verbose = T)
    
    if(is.null(experiments))
      endFailure()
    else {
      conditions <- enrich(experiments, options, verbose = T, inprod = T)
      session$userData$endTime <- Sys.time()
      
      if(is.null(conditions))
        endEmpty()
      else {
        conditions <- endSuccess(genes, experiments, conditions, options)
        
        # TODO for "any"
        geneInfo <- DATA.HOLDER[[options$taxa$value]]@gene.meta[entrez.ID %in% genes, .(entrez.ID, gene.Name)]
        
        output$results <- generateResults(conditions)
        
        # Prepare some plotting information.
        session$userData$plotData <- list(
          options = options,
          gene.ID = geneInfo$entrez.ID,
          gene.Name = geneInfo$gene.Name,
          experiments = experiments$rn,
          conditions = conditions[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, `Test Statistic`, `Ontology Steps`)]
        )
      }
    }
  }
  
  #' Search Genes
  #' 
  #' Begin the gene search process. Cleans the parameters as necessary and calls a function to
  #' perform the actual search and update the display.
  #'
  #' @param genes A character vector of genes (entrez ID, ensembl ID, name or keyword)
  #' @param signature A DE signature to search
  #' @param taxa The taxon ([human|mouse|rat])
  #' @param update Whether or not to update the browser's query string.
  searchGenes <- function(genes = input$genes,
                          signature = input$signature,
                          taxa = input$taxa,
                          update = T) {
    session$userData$startTime <- as.POSIXct(Sys.time())
    setProgress(environment(), 0, 'Validating input', n.steps = getOption('max.progress.steps'))
    
    signature <- unlist(strsplit(gsub(' ', ',', gsub('(;|,( ?))', ',', signature)), ','))
    
    options <- lapply(names(getConfig()), function(x) {
      ret <- getConfig(x)
      if(x == 'sig')
        v <- signature
      else
        v <- input[[x]]
      
      if(!is.null(v))
        ret$value <- v
      ret
    }) %>% `names<-`(names(getConfig()))
    
    # Update the query string
    if(update) {
      query <- paste0('?genes=', switch(min(2, length(genes)), genes, paste0('[', paste0(genes, collapse = ','), ']')))
      
      extras <- lapply(names(getConfig()), function(x) {
        if(x == 'sig')
          v <- signature
        else
          v <- input[[x]]
        
        if(!is.null(v) && v != '' && (!isTRUE(all.equal(v, getConfig(key = x)$value)) || length(v) != length(getConfig(key = x)$value))) {
          if(length(v) > 1)
            paste0(x, '=', paste0('[', paste0(v, collapse = ','), ']'))
          else
            paste0(x, '=', v)
        }
      }) %>% unlist %>% paste0(collapse = '&')
      
      if(nchar(extras) > 0)
        query <- paste0(query, '&', extras)
      
      updateQueryString(query, 'push')
    }
    
    if(is.null(taxa)) taxa <- options$taxa$value
    
    tidyGenes <- function(genes, taxa) {
      if(taxa == 'any') {
        taxOptions <- getConfig(key = 'taxa')$core
        taxIDs <- getConfig(key = 'taxa')$mapping
        
        orthologs <- lapply(taxOptions, function(i) {
          lapply(taxOptions, function(j) {
            homologene(genes, inTax = unname(taxIDs[i]), outTax = unname(taxIDs[j])) %>% {
              if(nrow(.) > 0) {
                if(i == j) {
                  .[, !duplicated(colnames(.))] %>% .[, grepl('_ID', colnames(.))] %>% {
                    data.frame(.) %>% `colnames<-`(paste0(unname(taxIDs[i]), '_ID')) %>%
                      rbind(
                        data.frame(
                          DATA.HOLDER[[i]]@gene.meta[gene.Name %in% genes | entrez.ID %in% genes, entrez.ID]
                        ) %>% `colnames<-`(paste0(unname(taxIDs[i]), '_ID'))
                      )
                  }
                } else
                  .[, grepl('_ID', colnames(.))]
              }
            }
          }) %>% rbindlist(fill = T) %>% {
            if(nrow(.) > 0) .
          }
        }) %>% rbindlist(fill = T) %>% .[, names(.) := lapply(.SD, as.character)] %>%
          .[, key := fcoalesce(.)] %>%
          group_by(key) %>% summarise_all(~na.omit(unique(.))[1]) %>%
          as.data.table
        
        return(orthologs)
      }
      
      oGenes <- genes
      # Clean numerics (interpreted as entrez IDs) and remove them from further processing.
      cleanGenes <- suppressWarnings(Filter(function(x) !is.na(as.integer(x)), genes))
      idMap <- sapply(as.character(cleanGenes), function(x) which(goGenes == x), USE.NAMES = F)
      
      genes <- genes[!(genes %in% cleanGenes)]
      
      # If it matches (ENSG|ENSMUS|ENSRNO)\d{11}, it's an Ensembl ID (for human, mouse or rat).
      if(length(genes) > 0) {
        ensembl <- grep('(ENSG|ENSMUS|ENSRNO)\\d{11}', genes, value = T)
        
        if(length(ensembl) != 0) {
          # Extract genes with a matching Ensembl ID and clean them too.
          ensembls <- DATA.HOLDER[[taxa]]@gene.meta[ensembl.ID %in% ensembl, .(entrez.ID, ensembl.ID)]
          cleanGenes <- c(cleanGenes, ensembls[, entrez.ID])
          idMap <- c(idMap, sapply(ensembls[, ensembl.ID], function(x) which(oGenes == x), USE.NAMES = F))
          
          genes <- genes[!(genes %in% ensembls[, ensembl.ID])]
        }
      }
      
      if(length(genes > 0)) {
        go <- grep('(GO:)\\d{7}', genes, value = T)
        
        if(length(go) != 0) {
          gos <- DATA.HOLDER[[taxa]]@go[id %in% go, .(id, entrez.ID)]
          cleanGenes <- c(cleanGenes, gos[, entrez.ID])
          idMap <- c(idMap, sapply(gos[, id], function(x) which(oGenes == x), USE.NAMES = F))
          
          genes <- genes[!(genes %in% gos[, id])]
        }
      }
      
      # Try to match to gene names and descriptions.
      if(length(genes) > 0) { # TODO this could grepl, and likely the idMap is bad
        descriptors <- DATA.HOLDER[[taxa]]@gene.meta[gene.Name %in% genes | gene.Desc %in% genes,
                                                     .(entrez.ID, gene.Name, gene.Desc)]
        if(nrow(descriptors) != 0) {
          cleanGenes <- c(cleanGenes, descriptors[, entrez.ID])
          idMap <- c(idMap, sapply(descriptors[, gene.Name], function(x) which(oGenes == x), USE.NAMES = F))
          idMap <- c(idMap, sapply(descriptors[, gene.Desc], function(x) which(oGenes == x), USE.NAMES = F))

          genes <- genes[!(genes %in% descriptors[, c(gene.Name, gene.Desc)])]
        }
      }
      
      # If anything is left, try to match it to gene aliases.
      if(length(genes) > 0) { # TODO this can multimatch and shouldn't, also no idMap yet
        aliases <- DATA.HOLDER[[taxa]]@gene.meta[, parseListEntry(alias.Name), entrez.ID] %>%
          .[grepl(paste0(genes, collapse = '|'), V1)]
        if(nrow(aliases) > 0)
          cleanGenes <- c(cleanGenes, aliases[, unique(entrez.ID)])
      }
      
      cleanGenes[order(unlist(idMap))]
    }
    
    # Done processing, handle the search.
    # TODO should prevent signature from working for taxa == 'all'
    handleSearch(tidyGenes(genes, taxa), options)
  }
}
