#' Server
#'
#' Query strings are supported as follows:
#' ?genes=[list of genes]&taxa=[human|mouse|rat]&...
server <- function(input, output, session) {
  session$onSessionEnded(function() {
    session$userData$INTERRUPT <- T
  })
  
  advanceProgress <- function(detail) {
    setProgress(session$userData$progress + 1, detail)
  }
  
  setProgress <- function(progress = NULL, detail = '', n.steps = NULL) {
    if(is.null(progress)) {
      if(is.null(session$userData$progress.bar)) return(NULL)
      progress <- session$userData$progress.bar$getMax()
      n.steps <- progress
    }
    
    session$userData$progress <- progress
    
    if(progress == 0) {
      session$userData$progress.bar <- shiny::Progress$new(min = 0, max = n.steps)
      session$userData$progress.bar$set(message = 'Working...', detail = detail)
    } else {
      session$userData$progress.bar$set(value = progress, detail = detail)
      
      if(progress >= session$userData$progress.bar$getMax()) {
        session$userData$progress.bar$close()
        shinyjs::enable('search')
        shinyjs::enable('reset')
      }
    }
  }
  
  # On connect, observe the query string. Search if any genes are specified
  observeEvent(input$LOAD, {
    output$results_header <- renderUI(generateResultsHeader(
      HTML('<div style="margin-bottom: 10px"><h2 style="display: inline">No enrichments yet</h2></div>')))
    
    query <- getQueryString(session)
    options <- getConfig()
    
    genes <- query$genes %>% jsonify
    taxa <- query$taxa %>% jsonify
    
    sig <- switch(is.null(query$sig) + 1, jsonify(query$sig), options$sig$value)
    fc <- switch(is.null(query$fc) + 1, jsonify(query$fc) %>% as.numeric, options$fc$value)
    categories <- switch(is.null(query$categories) + 1, jsonify(query$categories), options$categories$value)
    
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    updatePickerInput(session, 'taxa', selected = taxa)
    session$sendCustomMessage('querySet', genes)
    
    updateTextInput(session, 'sig', value = sig)
    updateSliderInput(session, 'fc', value = fc)
    updateSelectizeInput(session, 'categories', selected = categories)
    
    for(mName in Filter(function(x) !(x %in% c('sig', 'fc', 'categories', 'taxa')), names(options))) {
      do.update(session, options[[mName]], ifelse(is.null(query[[mName]]), options[[mName]]$value, query[[mName]]))
    }
    
    # Search if any genes are specified. Specify parameters to make sure no reactives.
    # Run it next ms to give the server time to propagate layout changes.
    if(!is.null(query$genes))
      shinyjs::delay(1, searchGenes(genes, sig, jsonify(query$taxa), update = F))
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
    updatePickerInput(session, 'taxa', selected = getConfig('taxa')$value)
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
        geneEvidence(session$userData$plotData$genes,
                     session$userData$plotData$options$taxa$value)$
          then(function(evidence) {
            output$results_genes <- generateGenePage(evidence)
          })
      })
    } else if(input$tabs == 'GO Enrichment' && !is.null(session$userData$plotData) && is.null(session$userData$goRendered)) {
      session$userData$goRendered <- T
      
      # TODO Only doing this on first taxon
      if(length(session$userData$plotData$options$taxa$value) > 1)
        mGenes <- mGenes[taxon == session$userData$plotData$options$taxa$value[1],
                         unique(gene.Name)]
      else
        mGenes <- DATA.HOLDER[[session$userData$plotData$options$taxa$value]]@gene.meta[
          entrez.ID %in% unique(session$userData$plotData$genes), gene.Name]
      
      output$results_go <- generateGOPage(ora(hitlist = mGenes,
                                              annotation = paste0('Generic_', session$userData$plotData$options$taxa$value[1])))
    } else if(input$tabs == 'Word Cloud' && !is.null(session$userData$plotData) && is.null(session$userData$cloudRendered)) {
      session$userData$cloudRendered <- T
      
      output$results_cloud <- generateResultsCloud(session$userData$plotData$conditions, session$userData$options)
    } else if(input$tabs == 'Gene Contributions' && !is.null(session$userData$plotData) && is.null(session$userData$contribsRendered)) {
      session$userData$contribsRendered <- T
      
      output$results_contribs <- generateGeneContribs(session$userData$plotData$conditions, session$userData$plotData$options)
    }
  })
  
  # Force an updated of the gene view when we get new search results
  observeEvent(input$UPDATED, {
    session$userData$genesRendered <- NULL
    session$userData$goRendered <- NULL
    session$userData$cloudRendered <- NULL
    session$userData$contribsRendered <- NULL
    
    if(input$tabs == 'Gene Info') {
      shinyjs::delay(100, {
        if(input$tabs == 'Gene Info' && !is.null(session$userData$plotData) && is.null(session$userData$genesRendered)) {
          session$userData$genesRendered <- T
          
          synchronise({
            geneEvidence(session$userData$plotData$genes,
                         session$userData$plotData$options$taxa$value)$
              then(function(evidence) {
                output$results_genes <- generateGenePage(evidence)
              })
          })
        }
      })
    } else if(input$tabs == 'GO Enrichment') {
      session$userData$goRendered <- T
      
      # TODO Only doing this on first taxon
      if(length(session$userData$plotData$options$taxa$value) > 1)
        mGenes <- mGenes[taxon == session$userData$plotData$options$taxa$value[1],
                         unique(gene.Name)]
      else
        mGenes <- DATA.HOLDER[[session$userData$plotData$options$taxa$value]]@gene.meta[
          entrez.ID %in% unique(session$userData$plotData$genes), gene.Name]
      
      output$results_go <- generateGOPage(ora(hitlist = mGenes,
                                              annotation = paste0('Generic_', session$userData$plotData$options$taxa$value[1])))
    } else if(input$tabs == 'Word Cloud' && !is.null(session$userData$plotData) && is.null(session$userData$cloudRendered)) {
      session$userData$cloudRendered <- T
      
      output$results_cloud <- generateResultsCloud(session$userData$plotData$conditions, session$userData$options)
    } else if(input$tabs == 'Gene Contributions' && !is.null(session$userData$plotData) && is.null(session$userData$contribsRendered)) {
      session$userData$contribsRendered <- T
      
      output$results_contribs <- generateGeneContribs(session$userData$plotData$conditions, session$userData$options)
    }
  })
  
  observeEvent(input$RANDOM_GENES, {
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    while(length(genes <- DATA.HOLDER[[ifelse(length(input$taxa) > 1, 'human', input$taxa)]]@gene.meta[sample(1:.N, sample(1:10, 1))] %>%
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
                      selectInput('plot_type', 'Type', list('Boxplot', 'Violin plot', 'Heatmap', 'Scatterplot', 'Jitterplot')),
                      selectInput('plot_data', 'Data', list('Gene Expression')),
                      pickerInput('plot_genes', 'Genes', list('Loading...'), multiple = T),
                      pickerInput('plot_conditions', 'Contrasts', list('Loading...'), multiple = T),
                      pickerInput('plot_taxa', 'Taxa', list('Loading...'), multiple = T))
      )
    ), session)
    
    generatePlot <- function(value) {
      session$userData$plotData$exprInfo <- value
      
      updatePickerInput(session, 'plot_taxa',
                        choices = as.list(session$userData$plotData$options$taxa$value),
                        selected = session$userData$plotData$options$taxa$value[1])
      
      if(length(session$userData$plotData$options$taxa$value) > 1)
        mGenes <- unique(session$userData$plotData$genes[, gene.realName])
      else
        mGenes <- DATA.HOLDER[[session$userData$plotData$options$taxa$value]]@gene.meta[entrez.ID %in% unique(session$userData$plotData$genes), gene.Name]
      
      updatePickerInput(session, 'plot_genes',
                        choices = as.list(intersect(mGenes, rownames(value$expr))),
                        selected = intersect(mGenes, rownames(value$expr)))
      updatePickerInput(session, 'plot_conditions', # TODO Selection strategy
                        choices = as.list(session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)]),
                        selected = data.table::first(session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)]))
      
      output$plot <- renderPlotly({
        generateResultsPlot(session$userData$plotData$genes,
                            session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)],
                            value,
                            session$userData$plotData$options,
                            input$plot_taxa, input$plot_genes, input$plot_conditions, input$plot_type, input$plot_data)
      })
    }
    
    if(is.null(session$userData$plotData$exprInfo)) {
      # If we didn't already process we need to give the UI a tick to update and then fetch the data.
      if(!is.null(session$userData$plotData$conditions)) {
        shinyjs::delay(1, { 
          synchronise({
            mExpData <- rbindlist(lapply(session$userData$plotData$options$taxa$value,
                                         function(i) DATA.HOLDER[[i]]@experiment.meta[, .(rsc.ID, ee.ID, ee.NumSamples)]))
            
            getTags(session$userData$plotData$options$taxa$value) %>% # TODO Selection strategy
              .[cf.Cat %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.Cat)]) &
                  cf.BaseLongUri %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.BaseLongUri)]) &
                  cf.ValLongUri %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.ValLongUri)])] %>%
              unique %>%
              merge(mExpData, by = c('rsc.ID', 'ee.ID'), sort = F, allow.cartesian = T) %>% {
                if(nrow(.) == 0) return(NULL)
                
                # Select some samples to queue for gene expression visualization
                while((tmp <- .[sample(1:nrow(.))])[1, ee.NumSamples] > getOption('max.gemma')) {
                }
                
                tmp[cumsum(ee.NumSamples) <= getOption('max.gemma')]
              } %>% {
                geneExpression(.[, unique(ee.ID)],
                               .[, unique(rsc.ID)],
                               session$userData$plotData$options$taxa$value,
                               session$userData$plotData$genes)$then(generatePlot)
              }
          })
        })
      }
    } else {
      # If we already processed then we can just display it
      output$plot <- renderPlotly({
        generateResultsPlot(session$userData$plotData$genes,
                            session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)],
                            session$userData$plotData$exprInfo,
                            session$userData$plotData$options,
                            input$plot_taxa, input$plot_genes, input$plot_conditions, input$plot_type, input$plot_data)
      })
    }
  })
  
  #' Display a failure message
  #'
  #' Ends the search protocol with a failure message
  endFailure <- function() {
    setProgress()
    output$results_header <- renderUI({
      generateResultsHeader(HTML('<h2 style="margin-top: 0;">No results</h2>'))
    })
    output$results <- NULL
  }
  
  makeSuggestions <- function(original, taxon, suggestions) {
    if(!is.null(suggestions)) {
      if(length(taxon) == 1) {
        suggested <- lapply(getConfig('taxa')$core, function(t) {
          i <- 0
          oGenes <- NA
          tmp <- tidyGenes(original, t)
          if(!is.null(tmp)) {
            i <- length(tmp$genes)
            oGenes <- tmp$genes
          }
          
          list(n = i,
               genes = oGenes,
               taxon = t)
        }) %>% {
          suggestedTaxon <- lapply(., '[[', 'n') %>% which.max
          list(taxon = .[[suggestedTaxon]]$taxon, genes = .[[suggestedTaxon]]$genes)
        }
        
        if(suggested$taxon != taxon) {
          output$results_suggestions <- renderUI({
            fluidRow(class = 'info-text',
                     column(12, HTML(paste0('Did you mean to search ', paste0('<a search taxon=\'', suggested$taxon, '\'>', suggested$taxon, '</a>'),
                                            ' experiments? Otherwise, try using <a search genes=\'[', paste0(paste0('"', tidyGenes(suggested$genes, c(taxon, suggested$taxon))[taxon != suggested$taxon, gene.realName], '"'), collapse = ','), ']\'>', taxon, ' orthologs', '</a> of these ', suggested$taxon, ' genes.'))))
          })
          shinyjs::delay(100, { shinyjs::runjs('loadExamples();') })
          
          return(NULL)
        }
      }
      
      revisedSearch <- original
      for(i in names(suggestions)) {
        revisedSearch <- gsub(i, suggestions[i], revisedSearch, fixed = T)
      }
      
      output$results_suggestions <- renderUI({
        fluidRow(class = 'info-text', column(12, HTML(paste0('Did you mean: ', paste0('<a search genes=\'[', paste0(paste0('"', revisedSearch, '"'), collapse = ','), ']\'>', paste0(suggestions, collapse = ', '), '</a>'), '?'))))
      })
      shinyjs::delay(100, { shinyjs::runjs('loadExamples();') })
    }
  }
  
  #' Display an empty message
  #'
  #' Ends the search protocol with no results
  endEmpty <- function() {
    setProgress()
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
    if(exists('output')) {
      n_exp <- sapply(options$taxa$value, function(i) nrow(DATA.HOLDER[[i]]@experiment.meta)) %>% sum
      
      output$results_header <- renderUI({
        generateResultsHeader(HTML(paste0('<div style="margin-bottom: 10px"><h2 style="display: inline">Enriched ',
                                          format(n_exp, big.mark = ','), ' condition comparison', ifelse(n_exp > 1, 's', ''), ' for ',
                                          ifelse(nrow(genes) == 1, genes[, gene.Name],
                                                 paste0('<span data-toggle="tooltip" data-placement="top" title="',
                                                        paste0(unique(genes[, gene.Name]), collapse = ', '), '">',
                                                        nrow(genes), ' gene', ifelse(nrow(genes) > 1, 's', ''),  '</span>')),
                                          ifelse(length(options$taxa$value) > 1,
                                                 paste0(' across <span data-toggle="tooltip" data-placement="top" title="', paste0(options$taxa$value, collapse = ', '), '">',
                                                        length(options$taxa$value), ' taxa</span>'), ''),
                                          '</h2><span class="timestamp">in ',
                                          format(difftime(session$userData$endTime, session$userData$startTime), digits = 3), '.</span></span>')))
      })
    }
    
    conditions[, `Condition Comparison` := stringi::stri_c('<b>', cf.BaseLongUri, '</b> vs. <b>', cf.ValLongUri, '</b>')]
    conditions[, distance := round(distance, 3)]
    
    advanceProgress('Cross-linking')
    
    tmp <- rbindlist(lapply(options$taxa$value, function(i) DATA.HOLDER[[i]]@experiment.meta[, .(rsc.ID, ee.ID, ee.Name)]))
    
    # Associating experiments with tags
    tmp <- getTags(options$taxa$value) %>%
      .[cf.Cat %in% conditions[, unique(cf.Cat)] &
          cf.BaseLongUri %in% conditions[, unique(cf.BaseLongUri)] &
          cf.ValLongUri %in% conditions[, unique(cf.ValLongUri)]] %>%
      merge(tmp, by = c('rsc.ID', 'ee.ID'), sort = F, allow.cartesian = T) %>%
      .[, N := length(unique(ee.ID)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
    
    tmp[, Evidence := ifelse(N == 1,
                             stringi::stri_c('<span data-id="', ee.ID[1], '>1 Experiment</span>'),
                             {
                               dedup <- !duplicated(ee.ID)
                               stringi::stri_c('<span data-id="', stringi::stri_c(ee.ID[dedup], collapse = ','),
                                               '" data-ee="', stringi::stri_c(ee.Name[dedup], collapse = ','), '">', N, ' Experiments</span>')
                             }),
        .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
    
    tmp[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, N, Evidence)] %>%
      unique %>%
      merge(conditions, by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), sort = F) %>%
      setnames(c('stat', 'score', 'distance'), c('Effect Size', 'Test Statistic', 'Ontology Steps'))
  }
  
  #' Handle Search
  #' 
  #' Call the processing function to handle searching and update the display.
  #'
  #' @param genes A character vector of entrez IDs
  #' @param options The search options
  handleSearch <- function(genes, rawGenes, options) {
    advanceProgress('Ranking experiments')
    
    # For some reason genes becomes a Promise and this resolves it?
    genes
    
    future({
      if(!is.null(genes) && length(options$taxa$value) > 1 && 'data.table' %in% class(genes)) {
        lapply(options$taxa$value, function(t) {
          mOp <- options
          mOp$taxa$value <- t
          search(genes[taxon == t, entrez.ID], mOp)
        })
      } else if(length(options$taxa$value) == 1)
        search(genes$genes, options)
    }, globals = 'DATA.HOLDER') %...>%(function(experiments) {
      if(!is.null(genes) && !('data.table' %in% class(genes)))
        makeSuggestions(rawGenes, options$taxa$value, genes$suggestions)
      
      if(!is.null(session$userData$INTERRUPT))
        return(NULL)
      
      if(all(sapply(experiments, is.null)))
        endFailure()
      else {
        advanceProgress('Enriching')
        
        future({
          if(length(options$taxa$value) == 1)
            copy(enrich(experiments, options))
          else {
            lapply(1:length(options$taxa$value), function(t) {
              mOp <- options
              mOp$taxa$value <- options$taxa$value[t]
              mSearchable <- experiments[[t]]
              
              enrich(mSearchable, mOp) %>%
                setnames(genes[taxon == options$taxa$value[t], gene.realName],
                         genes[taxon == options$taxa$value[t], identifier], skip_absent = T) # TODO ?
            }) %>% rbindlist(fill = T) %>%
              reorderTags3 %>%
              .[, lapply(.SD, mean, na.rm = T), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
          }
        }, globals = c('CACHE.BACKGROUND', 'NULLS')) %...>%(function(conditions) {
          if(!is.null(session$userData$INTERRUPT))
            return(NULL)
          
          session$userData$endTime <- Sys.time()
          
          if(is.null(conditions))
            endEmpty()
          else {
            if(length(options$taxa$value) > 1) {
              geneInfo <- genes %>% copy %>% setnames('identifier', 'gene.Name')
              mGenes <- genes %>% copy
            } else {
              geneInfo <- DATA.HOLDER[[options$taxa$value]]@gene.meta %>%
                .[, .(.I, entrez.ID, gene.Name)] %>%
                .[entrez.ID %in% genes$genes, .(I, entrez.ID, gene.Name)]
              mGenes <- genes$genes
            }
            
            conditions <- endSuccess(geneInfo, experiments, conditions, options)
            
            if(!is.null(session$userData$INTERRUPT))
              return(NULL)
            
            output$results <- generateResults(conditions)
            
            # Prepare some plotting information.
            session$userData$plotData <- list(
              options = options,
              genes = mGenes,
              conditions = conditions %>% copy
            )
            
            output$dataDownload <- downloadHandler(function() {
              paste0('enrichment-', Sys.Date(), '.csv')
            }, function(file) {
              write.csv(session$userData$plotData$conditions %>%
                          .[, !c('Condition Comparison', 'Evidence')] %>%
                          setorder(-`Test Statistic`, `Ontology Steps`) %>%
                          setnames(c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), c('Category', 'Baseline', 'Value')), file)
            }, 'text/csv')
          }
        })
      }
    })
  }
  
  tidyGenes <- function(genes, taxa) {
    if(length(taxa) > 1) {
      taxIDs <- getConfig(key = 'taxa')$mapping[taxa]
      
      symbols <- lapply(taxa, function(x) {
        data.table(entrez.ID = tidyGenes(genes, x)$genes, taxon = x) %>% {
          if(is.null(.) || ncol(.) == 1) NULL
          else
            merge(., DATA.HOLDER[[x]]@gene.meta[, .(.I, entrez.ID, gene.Name)], by = 'entrez.ID')
        }
      }) %>% rbindlist
      
      if(!is.null(symbols) && nrow(symbols) > 0) {
        symbols <- symbols %>%
          merge(data.table(taxon.ID = taxIDs, taxon = names(taxIDs)), by = 'taxon') %>%
          .[, .SD[1], gene.Name]
      } else
        return(tidyGenes(genes, taxa[1]))
      
      orthologs <- symbols[, list(lapply(taxIDs[taxIDs != unique(taxon.ID)], function(t) {
        homologs <- homologene(entrez.ID, unique(taxon.ID), t)
        if(nrow(homologs) == 0) NULL
        else {
          homologs %>%
            as.data.table %>%
            melt(measure.vars = which(!grepl('_ID', colnames(.))), value.name = 'gene.realName', variable.name = 'taxon.ID') %>%
            .[, entrez.ID := unlist(lapply(1:nrow(.), function(I) {
              .SD[I, grepl(paste0(taxon.ID[I], '_ID'), colnames(.SD)), with = F]
            }))] %>%
            .[, !grepl('_ID', colnames(.)), with = F] %>%
            .[, taxon.ID := as.integer(levels(taxon.ID))[taxon.ID]] %>%
            cbind(identifier = unique(homologs[, as.character(unique(taxon.ID))])) %>%
            merge(DATA.HOLDER[[names(taxIDs)[taxIDs == t]]]@gene.meta[, .(.I, entrez.ID = as.integer(entrez.ID))], by = 'entrez.ID')
        }
      })), taxon] %>% .[, V1] %>% rbindlist %>%
        rbind(symbols[, .(entrez.ID, taxon.ID, gene.realName = gene.Name, identifier = gene.Name, I)]) %>%
        unique %>%
        merge(data.table(taxon.ID = taxIDs, taxon = names(taxIDs)), by = 'taxon.ID') %>%
        .[, entrez.ID := as.character(entrez.ID)] %>%
        .[, identifier := make.names(identifier, T), taxon.ID]
      
      return(orthologs)
    }
    
    oGenes <- genes
    # Clean numerics (interpreted as entrez IDs) and remove them from further processing.
    cleanGenes <- suppressWarnings(Filter(function(x) !is.na(as.integer(x)), genes))
    idMap <- sapply(as.character(cleanGenes), function(x) which(oGenes == x), USE.NAMES = F)
    
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
    
    # Try to match to gene names
    if(length(genes) > 0) {
      descriptors <- DATA.HOLDER[[taxa]]@gene.meta[gene.Name %in% genes, .(entrez.ID, gene.Name)]
      if(nrow(descriptors) != 0) {
        cleanGenes <- c(cleanGenes, descriptors[, entrez.ID])
        idMap <- c(idMap, sapply(descriptors[, gene.Name], function(x) which(oGenes == x), USE.NAMES = F))
        
        genes <- genes[!(genes %in% descriptors[, gene.Name])]
      }
    }
    
    if(F) {
      # If anything is left, try to match it to gene aliases.
      if(length(genes) > 0) { # TODO this can multimatch and shouldn't, also no idMap yet
        aliases <- DATA.HOLDER[[taxa]]@gene.meta[, parseListEntry(alias.Name), entrez.ID] %>%
          .[grepl(paste0(genes, collapse = '|'), V1)]
        if(nrow(aliases) > 0)
          cleanGenes <- c(cleanGenes, aliases[, unique(entrez.ID)])
      }
    }
    
    mSuggestions <- NULL
    if(length(genes) > 0) {
      mSuggestions <- sapply(genes, function(gene) {
        stringdist(tolower(gene), tolower(DATA.HOLDER[[taxa]]@gene.meta[, gene.Name])) %>% {
          DATA.HOLDER[[taxa]]@gene.meta[which.min(.), gene.Name]
        }
      }) %>% setNames(genes)
    }
    
    mGenes <- NULL
    if(length(cleanGenes) > 0)
      mGenes <- cleanGenes[order(unlist(idMap))]
    
    if(length(mGenes) > 0 || length(mSuggestions) > 0)
      list(genes = mGenes, suggestions = mSuggestions)
    else
      NULL
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
    output$results_suggestions <- NULL
    
    shinyjs::disable('search')
    shinyjs::disable('reset')
    
    session$userData$startTime <- as.POSIXct(Sys.time())
    setProgress(0, 'Validating input', n.steps = getOption('max.progress.steps'))
    
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
    
    # Done processing, handle the search.
    # TODO should prevent signature from working for length(taxa) > 1
    # TODO This { . } is a weird way to prevent lazy evaluation. Probably a better way
    tidyGenes(genes, taxa) %>% { . } %>% handleSearch(genes, options)
  }
}
