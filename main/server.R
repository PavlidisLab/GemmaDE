#' Server
#'
#' Query strings are supported as follows:
#' ?genes=[list of genes]&taxa=[human|mouse|rat]&scope=[list of ontologies]...
server <- function(input, output, session) {
  # On connect, observe the query string. Search if any genes are specified
  observeEvent(input$LOAD, {
    output$results_header <- renderUI(generateResultsHeader(
      HTML('<div style="margin-bottom: 10px"><h2 class="loading" style="display: inline">EnriChing</h2></div>')))
    
    query <- getQueryString(session)
    
    genes <- query$genes %>% jsonify
    scope <- switch(is.null(query$scope) + 1, jsonify(query$scope), getOption('app.ontology'))
    
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    session$sendCustomMessage('querySet', genes)
    updateSelectizeInput(session, 'taxa', selected = ifelse(is.null(query$taxa), getOption('app.taxa'), query$taxa))
    updateSelectizeInput(session, 'method', selected = ifelse(is.null(query$method), getOption('app.search_method'), query$method))
    updatePickerInput(session, 'scope', selected = scope)
    updateCheckboxInput(session, 'mfx', value = ifelse(is.null(query$mfx), getOption('app.mfx'), query$mfx))
    updateNumericInput(session, 'reqall', value = ifelse(is.null(query$required), getOption('app.req.all'), query$required))
    updateCheckboxInput(session, 'meanval', value = ifelse(is.null(query$expr), getOption('app.meanval'), query$expr))
    updateNumericInput(session, 'distance', value = ifelse(is.null(query$distance), getOption('app.distance_cutoff'), query$distance))
    updateNumericInput(session, 'min.tags', value = ifelse(is.null(query$min), getOption('app.min.tags'), query$min))
    updateNumericInput(session, 'max.rows', value = ifelse(is.null(query$rows), getOption('app.max.rows'), query$rows))
    updateNumericInput(session, 'pv', value = ifelse(is.null(query$pv), getOption('app.pv'), query$pv))
    updateSliderInput(session, 'fc', value = ifelse(is.null(query$fc), c(getOption('app.fc_lower'), getOption('app.fc_upper')), query$fc %>% jsonify %>% as.numeric))
    
    # Search if any genes are specified. Specify parameters to make sure no reactives.
    # Run it next ms to give the server time to propagate layout changes.
    if(!is.null(query$genes))
      shinyjs::delay(1, searchGenes(genes, query$taxa, scope, update = F))
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
    updateSelectizeInput(session, 'taxa', selected = NULL)
    updateSelectizeInput(session, 'method', selected = NULL)
    updatePickerInput(session, 'scope', selected = getOption('app.ontology'))
    updateCheckboxInput(session, 'mfx', value = getOption('app.mfx'))
    updateNumericInput(session, 'reqall', value = getOption('app.req.all'))
    updateCheckboxInput(session, 'meanval', value = getOption('app.meanval'))
    updateCheckboxInput(session, 'geeq', value = getOption('app.geeq'))
    updateNumericInput(session, 'min.tags', value = getOption('app.min.tags'))
    updateNumericInput(session, 'distance', value = getOption('app.distance_cutoff'))
    updateNumericInput(session, 'max.rows', value = getOption('app.max.rows'))
    updateNumericInput(session, 'pv', value = getOption('app.pv'))
    updateSliderInput(session, 'fc', value = c(getOption('app.fc_lower'), getOption('app.fc_upper')))
    session$sendCustomMessage('fileUpload', F)
  })
  
  # Search
  observeEvent(input$search, {
    genes <- input$genes
    if(!is.null(input$genes.csv))
      genes <- read.csv(input$genes.csv$datapath, header = F)$V1 %>% as.character
    
    searchGenes(genes)
  })
  
  # Open gene or GO tab
  observeEvent(input$tabs, {
    if(input$tabs == 'Gene Info' && !is.null(session$userData$plotData) && is.null(session$userData$genesRendered)) {
      session$userData$genesRendered <- T
      
      synchronise({
        geneEvidence(session$userData$plotData$gene.ID,
                     session$userData$plotData$taxa)$
          then(function(evidence) {
            output$results_genes <- generateGenePage(evidence)
          })
      })
    } else if(input$tabs == 'GO Enrichment' && !is.null(session$userData$plotData) && is.null(session$userData$goRendered)) {
      session$userData$goRendered <- T
      
      # We could precompute this but it seems to run rapidly and the other view is more convenient.
      #annList <- DATA.HOLDER[[session$userData$plotData$taxa]]@go[, list(list(as.character(id))), as.character(entrez.ID)] %>%
      #  lapply(c) %>% { `names<-`(.[[2]], .[[1]]) }
      
      #annotation <- makeAnnotation(
      #  annList,
      #  name = DATA.HOLDER[[session$userData$plotData$taxa]]@gene.meta[, .(entrez.ID, gene.Name)] %>%
      #    as.data.frame %>% `rownames<-`(.[, 'entrez.ID']) %>%
      #    .[names(annList), ]
      #)
      
      output$results_go <- generateGOPage(ora(hitlist = session$userData$plotData$gene.Name,
                                              annotation = paste0('Generic_', session$userData$plotData$taxa))) # annotation))
    }
  })
  
  # Force an updated of the gene view when we get new search results
  observeEvent(input$UPDATED, {
    session$userData$genesRendered <- NULL
    session$userData$goRendered <- NULL

    if(input$tabs == 'Gene Info') {
      shinyjs::delay(100, {
        if(input$tabs == 'Gene Info' && !is.null(session$userData$plotData) && is.null(session$userData$genesRendered)) {
          session$userData$genesRendered <- T
          
          synchronise({
            geneEvidence(session$userData$plotData$gene.ID,
                         session$userData$plotData$taxa)$
              then(function(evidence) {
                output$results_genes <- generateGenePage(evidence)
              })
          })
        }
      })
    } else if(input$tabs == 'GO Enrichment') {
      session$userData$goRendered <- T
      
      # We could precompute this but it seems to run rapidly and the other view is more convenient.
      #annList <- DATA.HOLDER[[session$userData$plotData$taxa]]@go[, list(list(as.character(id))), as.character(entrez.ID)] %>%
      #  lapply(c) %>% { `names<-`(.[[2]], .[[1]]) }
      
      #annotation <- makeAnnotation(
      #  annList,
      #  name = DATA.HOLDER[[session$userData$plotData$taxa]]@gene.meta[, .(entrez.ID, gene.Name)] %>%
      #    as.data.frame %>% `rownames<-`(.[, 'entrez.ID']) %>%
      #    .[names(annList), ]
      #)
      
      output$results_go <- generateGOPage(ora(hitlist = session$userData$plotData$gene.Name,
                                              annotation = paste0('Generic_', session$userData$plotData$taxa))) # annotation))
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
               column(10, plotlyOutput('plot') %>% withSpinner),
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
      #rownames(value$expr) <- merge(data.table(ID = rownames(value$expr)),
      #                          data.table(ID = session$userData$plotData$gene.ID,
      #                                     name = session$userData$plotData$gene.Name),
      #                          by = 'ID', sort = F, all.x = T)[, name]
      
      session$userData$plotData$processed <- value
      session$userData$plotData$queued <- NULL
      
      updatePickerInput(session, 'plot_genes',
                        choices = as.list(intersect(session$userData$plotData$gene.Name, rownames(value$expr))),
                        selected = intersect(session$userData$plotData$gene.Name, rownames(value$expr)))
      updatePickerInput(session, 'plot_conditions',
                        choices = as.list(session$userData$plotData$conditions),
                        selected = data.table::first(session$userData$plotData$conditions))
      
      output$plot <- renderPlotly({
        generateResultsPlot(session$userData$plotData$gene.Name,
                            session$userData$plotData$conditions,
                            value,
                            session$userData$plotData$options,
                            input$plot_genes, input$plot_conditions, input$plot_type, input$plot_data)
      })
    }
    
    if(is.null(session$userData$plotData$processed)) {
      # If we didn't already process we need to give the UI a tick to update and then
      # fetch the data.
      if(!is.null(session$userData$plotData$queued)) {
        shinyjs::delay(1, {
          synchronise(session$userData$plotData$queued %>% {
            geneExpression(.[, unique(ee.ID)],
                           .[, unique(rsc.ID)],
                           session$userData$plotData$taxa,
                           session$userData$plotData$gene.ID)$then(generatePlot)
          })
        })
      }
    } else {
      # If we already processed then we can just display
      output$plot <- renderPlotly({
        generateResultsPlot(session$userData$plotData$gene.Name,
                            session$userData$plotData$conditions,
                            session$userData$plotData$processed,
                            session$userData$plotData$options,
                            input$plot_genes, input$plot_type, input$plot_data)
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
  #' @param taxa The taxa scope that was queried
  #' @param scope The ontology scope that was queried
  #' @param options Any additional options that were passed
  endSuccess <- function(genes, experiments, conditions,
                         taxa = getOption('app.taxa'), scope = getOption('app.ontology'),
                         options = getOption('app.all_options')) {
    # Generate the results header
    if(exists('output'))
      output$results_header <- renderUI({
        generateResultsHeader(HTML(paste0('<div style="margin-bottom: 10px"><h2 style="display: inline">Found ',
                                          nrow(experiments), ' experiment', ifelse(nrow(experiments) > 1, 's', ''), ' differentially expressing ',
                                          ifelse(ncol(experiments) == 2, colnames(experiments)[1],
                                                 paste0('<span data-toggle="tooltip" data-placement="top" title="',
                                                        paste0(colnames(experiments)[1:(ncol(experiments) - 3)], collapse = ', '), '">',
                                                        ncol(experiments) - 3, ' gene', ifelse(ncol(experiments) > 4, 's', ''), '</span>')),
                                          '</h2><span class="timestamp">in ',
                                          format(difftime(Sys.time(), session$userData$startTime), digits = 3), '.</span></span>')))
      })
    
    advanceProgress('Fetching Gemma data')
    
    # Compute experiments in scope
    mCache <- getTags(taxa, scope, experiments$rn, options$distance) %>%
      .[cf.BaseLongUri %in% conditions[, unique(cf.BaseLongUri)] & cf.ValLongUri %in% conditions[, unique(cf.ValLongUri)]] %>%
      unique %>% merge(DATA.HOLDER[[taxa]]@experiment.meta[rsc.ID %in% experiments$rn, .(rsc.ID, ee.ID)], by = 'rsc.ID', sort = F, allow.cartesian = T) %>%
      .[, N := length(unique(ee.ID)), .(cf.BaseLongUri, cf.ValLongUri)] %>%
      merge(data.table(rsc.ID = experiments$rn, experiments[, !c('score', 'rn')]), by = 'rsc.ID')
    
    # @see[1] Currently only supporting visualizing conditions that are significant.
    geneExpr <- mCache[paste0(cf.BaseLongUri, cf.ValLongUri) %in% conditions[pv.fisher <= options$pv, paste0(cf.BaseLongUri, cf.ValLongUri)]]
    
    advanceProgress('Adding cross references')
    
    geneScores <- data.table(conditions[, .(cf.BaseLongUri, cf.ValLongUri)]) %>%
      merge(mCache[, !'rsc.ID'], by = c('cf.BaseLongUri', 'cf.ValLongUri'), sort = F, allow.cartesian = T) %>%
      .[, Evidence := paste0(unique(ee.ID), collapse = ','), .(cf.BaseLongUri, cf.ValLongUri)]
    
    geneScores <- geneScores[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, N, Evidence)] %>%
      merge(
        merge(geneScores[, !c('ee.ID', 'cf.Cat', 'Evidence', 'N', 'direction')] %>%
                .[, lapply(.SD, function(x) paste0(head(round(x, 3), 100), collapse = ',')),
                  .(cf.BaseLongUri, cf.ValLongUri)],
              geneScores[, .(cf.BaseLongUri, cf.ValLongUri, direction = ifelse(direction, 1, -1) * ifelse(reverse, -1, 1))] %>%
                .[, lapply(.SD, function(x) paste(sum(x == 1), sum(x == -1), sep = ',')),
                  .(cf.BaseLongUri, cf.ValLongUri)],
              by = c('cf.BaseLongUri', 'cf.ValLongUri'), sort = F, allow.cartesian = T),
        by = c('cf.BaseLongUri', 'cf.ValLongUri'), sort = F, allow.cartesian = T) %>% unique
    
    # Rename things and add the ES
    conditions <- merge(conditions, geneScores, by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), sort = F, allow.cartesian = T) %>% setorder(pv.fisher) %>%
      setnames(c('pv.fisher', 'C', 'direction'), c('P-value', 'Augmented Count', 'Direction')) %>%
      .[, `Augmented Count` := round(`Augmented Count`, 2)]
    
    advanceProgress('Finishing up')
    
    list(results = conditions, expression = geneExpr)
  }
  
  #' Handle Search
  #' 
  #' Call the processing function to handle searching and update the display.
  #'
  #' @param genes A character vector of entrez IDs
  #' @param taxa The taxon of interest.
  #' @param scope The ontology scope.
  #' @param options The search options
  handleSearch <- function(genes, taxa = getOption('app.taxa'), scope = getOption('app.ontology'), options = getOption('app.all_options')) {
    experiments <- search(genes, taxa, options)
    
    if(is.null(experiments))
      endFailure()
    else {
      conditions <- enrich(experiments, taxa, scope, options)
      
      if(is.null(conditions))
        endEmpty()
      else {
        results <- endSuccess(genes, experiments, head(conditions, options$max.rows), taxa, scope, options)
        
        output$results <- generateResults(results$results, taxa, options)
        
        # Prepare some plotting information.
        session$userData$plotData <- list(
          taxa = taxa,
          gene.ID = genes,
          gene.Name = DATA.HOLDER[[taxa]]@gene.meta[entrez.ID %in% genes, gene.Name],
          # @see[1] Make sure to update this if we change the strategy to select viz-able
          conditions = conditions[pv.fisher <= options$pv, paste(cf.BaseLongUri, 'vs.', cf.ValLongUri)],
          options = options,
          queued = results$expression %>%
            merge(DATA.HOLDER[[taxa]]@experiment.meta[, .(rsc.ID, ee.ID, ee.NumSamples)],
                  by = c('rsc.ID', 'ee.ID'), all.x = T, sort = F) %>% {
                    if(nrow(.) == 0) return(NULL)
                    
                    # Select some samples to queue for gene expression visualization
                    while((tmp <- .[sample(1:nrow(.))])[1, ee.NumSamples] > getOption('max.gemma')) {
                    }
                    tmp[cumsum(ee.NumSamples) <= getOption('max.gemma')]
                  }
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
  #' @param taxa The taxon ([human|mouse|rat])
  #' @param scope The ontology scope
  #' @param update Whether or not to update the browser's query string.
  searchGenes <- function(genes = input$genes,
                          taxa = input$taxa,
                          scope = input$scope,
                          update = T) {
    session$userData$startTime <- as.POSIXct(Sys.time())
    setProgress(environment(), 0, 'Validating input', n.steps = getOption('max.progress.steps'))
    
    options <- list(pv = input$pv,
                    fc.lower = input$fc[1], fc.upper = input$fc[2],
                    score.lower = input$score[1], score.upper = input$score[2],
                    mfx = input$mfx,
                    reqall = input$reqall,
                    meanval = input$meanval,
                    geeq = input$geeq,
                    distance = input$distance,
                    min.tags = input$min.tags,
                    max.rows = input$max.rows,
                    method = input$method)
    
    # Update the query string
    if(update) {
      query <- paste0('?genes=',
                      switch(min(2, length(genes)), genes, paste0('[', paste0(genes, collapse = ','), ']')),
                      switch((taxa == getOption('app.taxa')) + 1, paste0('&taxa=', taxa), ''),
                      switch(identical(scope, getOption('app.ontology')) + 1,
                             paste0('&scope=', switch(min(2, length(scope)), scope,
                                                      paste0('[', paste0(scope, collapse = ','), ']'), ''))),
                      switch((options$pv == getOption('app.pv')) + 1, paste0('&pv=', options$pv), ''),
                      switch(isTRUE(all.equal(input$fc, c(getOption('app.fc_lower'), getOption('app.fc_upper')))) + 1,
                             paste0('&fc=[', paste0(input$fc, collapse = ','), ']'), ''),
                      switch((options$mfx == getOption('app.mfx')) + 1, paste0('&mfx=', options$mfx), ''),
                      switch((options$reqall == getOption('app.req.all')) + 1, paste0('&required=', options$reqall), ''),
                      switch((options$meanval == getOption('app.meanval')) + 1, paste0('&expr=', options$meanval), ''),
                      switch((options$geeq == getOption('app.geeq')) + 1, paste0('&geeq=', options$geeq), ''),
                      switch((options$min.tags == getOption('app.min.tags')) + 1, paste0('&min=', options$min.tags), ''),
                      switch((options$max.rows == getOption('app.max.rows')) + 1, paste0('&rows=', options$max.rows), ''),
                      switch((options$distance == getOption('app.distance_cutoff')) + 1, paste0('&distance=', options$distance), ''),
                      switch((options$method == getOption('app.search_method')) + 1, paste0('&method=', options$method), ''))
      updateQueryString(query, 'push')
    }
    
    if(is.null(taxa)) taxa <- getOption('app.taxa')
    if(is.null(scope)) scope <- getOption('app.ontology')
    
    tidyGenes <- function(genes, taxa) {
      if(taxa == 'any') {
        taxOptions <- Filter(function(x) !(x %in% c('artificial', 'dope', 'any')), getOption('app.all_taxa'))
        taxIDs <- c(human = 9606, mouse = 10090, rat = 10116)
        
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
      
      # Clean numerics (interpreted as entrez IDs) and remove them from further processing.
      cleanGenes <- suppressWarnings(Filter(~!is.na(as.integer(.)), genes))
      genes <- genes[!(genes %in% cleanGenes)]
      
      # If it matches (ENSG|ENSMUS|ENSRNO)\d{11}, it's an Ensembl ID (for human, mouse or rat).
      if(length(genes) > 0) {
        ensembl <- grep('(ENSG|ENSMUS|ENSRNO)\\d{11}', genes, value = T)
        
        if(length(ensembl) != 0) {
          # Extract genes with a matching Ensembl ID and clean them too.
          ensembls <- DATA.HOLDER[[taxa]]@gene.meta[ensembl.ID %in% ensembl, .(entrez.ID, ensembl.ID)]
          cleanGenes <- c(cleanGenes, ensembls[, entrez.ID])
          genes <- genes[!(genes %in% ensembls[, ensembl.ID])]
        }
      }
      
      if(length(genes > 0)) {
        go <- grep('(GO:)\\d{7}', genes, value = T)
        
        if(length(go) != 0) {
          gos <- DATA.HOLDER[[taxa]]@go[id %in% go, entrez.ID]
          cleanGenes <- c(cleanGenes, gos)
          genes <- genes[!(genes %in% go)]
        }
      }
      
      # Try to match to gene names and descriptions.
      if(length(genes) > 0) {
        descriptors <- DATA.HOLDER[[taxa]]@gene.meta[gene.Name %in% genes | gene.Desc %in% genes,
                                                     .(entrez.ID, gene.Name, gene.Desc)]
        if(nrow(descriptors) != 0) {
          cleanGenes <- c(cleanGenes, descriptors[, entrez.ID])
          genes <- genes[!(genes %in% descriptors[, c(gene.Name, gene.Desc)])]
        }
      }
      
      # If anything is left, try to match it to gene aliases.
      if(length(genes) > 0) {
        aliases <- DATA.HOLDER[[taxa]]@gene.meta[, parseListEntry(alias.Name), entrez.ID] %>%
          .[grepl(genes, V1)]
        if(nrow(aliases) > 0)
          cleanGenes <- c(cleanGenes, aliases[, unique(entrez.ID)])
      }
      
      cleanGenes
    }
    
    cleanGenes <- tidyGenes(genes, taxa)
    
    # Done processing, handle the search.
    handleSearch(cleanGenes, taxa, scope, options)
  }
}
