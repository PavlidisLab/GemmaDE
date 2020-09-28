library(shiny)
library(shinyjs)
library(htmlwidgets)
library(data.table)
library(dplyr)
library(shinycssloaders)

fineProgress <- function(session, steps) {
  if(!is.null(session))
    session$userData$progress.bar.steps <- steps
}

advanceProgress <- function(session, detail, fine = F) {
  if(!is.null(session))
    setProgress(session, session$userData$progress + 1, detail, fine)
}

setProgress <- function(session, progress, detail = '', fine = F) {
  if(!fine)
    session$userData$progress <- progress
  
  if(fine)
    session$userData$progress.bar$inc(amount = 1 / session$userData$progress.bar.steps)
  else if(progress == 0) {
    session$userData$progress.bar <- shiny::Progress$new(min = 0, max = getOption('max.progress.steps'))
    session$userData$progress.bar$set(message = 'Searching...', detail = detail)
  } else {
    session$userData$progress.bar$set(value = progress, detail = detail)
    
    if(progress >= getOption('max.progress.steps'))
      session$userData$progress.bar$close()
  }
}

#' Server
#'
#' Query strings are supported as follows:
#' ?genes=[list of genes]&taxa=[human|mouse|rat]&scope=[list of ontologies]...
server <- function(input, output, session) {
  # On connect, observe the query string. Search if any genes are specified
  observeEvent(input$LOAD, {
    query <- getQueryString(session)
    
    genes <- query$genes %>% jsonify
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    session$sendCustomMessage('querySet', genes)
    updateSelectizeInput(session, 'taxa', selected = query$taxa)
    updateCheckboxInput(session, 'mfx', value = ifelse(is.null(query$mfx), getOption('app.mfx'), query$mfx))
    updateNumericInput(session, 'distance', value = ifelse(is.null(query$distance), getOption('app.distance_cutoff'), query$distance))
    updateNumericInput(session, 'pv', value = ifelse(is.null(query$pv), getOption('app.pv'), query$pv))
    updateSliderInput(session, 'fc', value = ifelse(is.null(query$fc), c(getOption('app.fc_lower'), getOption('app.fc_upper')), query$fc %>% jsonify %>% as.numeric))
    
    # Use DO as a default scope, if unspecified
    scope <- query$scope %>% jsonify
    if(is.null(scope)) scope <- 'DO'
    
    updateSelectizeInput(session, 'scope', selected = scope)
    
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
    updateSelectizeInput(session, 'scope', selected = getOption('app.ontology'))
    updateCheckboxInput(session, 'mfx', value = getOption('app.mfx'))
    updateCheckboxInput(session, 'geeq', value = getOption('app.geeq'))
    updateNumericInput(session, 'distance', value = getOption('app.distance_cutoff'))
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
  
  # Open plot
  observeEvent(input$plotData, {
    showModal(modalDialog(
      easyClose = T,
      size = 'l',
      footer = NULL,
      fluidRow(style = 'height: 85vh;',
               column(10, plotlyOutput('plot') %>% withSpinner),
               column(2, style = 'padding-top: 15px; padding-right: 30px;',
                      fluidRow(style = 'display: flex; flex-direction: row; justify-content: space-evenly; margin-bottom: 15px;',
                               downloadButton('plot_save_jpg', 'Save JPG'), downloadButton('plot_save_pdf', 'Save PDF')),
                      checkboxInput('plot_page', 'Only current page', T),
                      selectInput('plot_type', 'Type', list('Heatmap', 'Histogram', 'Bars', 'Lines')),
                      selectInput('plot_data', 'Data', list('Gene Score', 'P-value')),
                      selectInput('plot_by', 'By', list('Condition', 'Experiment')))
      )
    ), session)
  })
  
  #' Handle Search
  #' 
  #' Call the processing function to handle searching and update the display.
  #'
  #' @param genes A character vector of entrez IDs
  #' @param cache A subsetted version of CACHE.BACKGROUND for the taxa scope.
  #' @param data A subsetted version of DATA.HOLDER for the taxa scope.
  #' @param taxa The taxon ([human|mouse|rat])
  #' @param options The search options
  #' @param scope The ontology scope
  #' @param options The search options
  handleSearch <- function(genes, cache, data, taxa, scope, options = getOption('app.all_options')) {
    #future({
    experiments <- search(data, genes, taxa, options, session)
    
    if(is.null(experiments)) {
      setProgress(session, getOption('max.progress.steps'))
      output$results_header <- renderUI({
        generateResultsHeader(HTML('<h2 data-toggle="tooltip" data-placement="top" title="Relax thresholds or modify gene set.">Invalid search.</h2>'))
      })
      output$results <- NULL
    } else {
      conditions <- enrich(cache, experiments, taxa, scope, options, session)
      
      if(is.null(conditions)) {
        setProgress(session, getOption('max.progress.steps'))
        output$results_header <- renderUI({ generateResultsHeader('No conditions found in scope.') })
        output$results <- NULL
      } else {
        # Generate the results header
        output$results_header <- renderUI({
          generateResultsHeader(HTML(paste0('<div style="margin-bottom: 10px"><h2 style="display: inline">Found ',
                                            nrow(experiments), ' experiment', ifelse(nrow(experiments) > 1, 's', ''), ' differentially expressing ',
                                            ifelse(ncol(experiments) == 2, colnames(experiments)[1],
                                                   paste0('<span data-toggle="tooltip" data-placement="top" title="',
                                                          paste0(colnames(experiments)[1:(ncol(experiments) - 2)], collapse = ', '), '">',
                                                          ncol(experiments) - 2, ' gene', ifelse(ncol(experiments) > 3, 's', ''), '</span>')),
                                            '</h2><span class="timestamp">in ',
                                            format(difftime(Sys.time(), session$userData$startTime), digits = 3), '.</span></span>')))
        })
        
        advanceProgress(session, 'Fetching Gemma data')
        
        # Compute experiments in scope
        mCache <- getTags(cache, taxa, scope, rownames(experiments), options$distance) %>%
          .[cf.BaseLongUri %in% conditions[, unique(cf.BaseLongUri)] & cf.ValLongUri %in% conditions[, unique(cf.ValLongUri)]] %>%
          unique %>% merge(data@experiment.meta[rsc.ID %in% rownames(experiments), .(rsc.ID, ee.ID)], by = 'rsc.ID', sort = F, allow.cartesian = T) %>%
          .[, N := length(unique(ee.ID)), .(cf.BaseLongUri, cf.ValLongUri)] %>%
          merge(data.table(rsc.ID = rownames(experiments), experiments[, !'score']), by = 'rsc.ID')
        
        advanceProgress(session, 'Adding cross references')
        
        geneScores <- data.table(conditions[, .(cf.BaseLongUri, cf.ValLongUri)]) %>%
          merge(mCache[, !'rsc.ID'], by = c('cf.BaseLongUri', 'cf.ValLongUri'), sort = F, allow.cartesian = T) %>%
          .[, Evidence := paste0(head(unique(ee.ID), 20), collapse = ','), .(cf.BaseLongUri, cf.ValLongUri)]
        
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
        conditions <- merge(conditions, geneScores, by = c('cf.BaseLongUri', 'cf.ValLongUri'), sort = F, allow.cartesian = T) %>% setorder(pv.fisher) %>%
          setnames(c('pv.fisher', 'direction'), c('P-value', 'Direction'))
        
        advanceProgress(session, 'Finishing up')
        
        output$plot <- renderPlotly({ generateResultsPlot(taxa, scope, experiments, conditions, options, input) })
        output$results <- generateResults(data, taxa, scope, experiments, head(conditions, getOption('max.rows')), options)
      }
    }
    #}, globals = c(as.vector(lsf.str()), 'ONTOLOGIES', 'ONTOLOGIES.DEFS', 'ui',
    #               'genes', 'cache', 'data', 'taxa', 'scope', 'options'))
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
    setProgress(session, 0, 'Validating input')
    
    options <- list(n.display = input$`top-n`,
                    pv = input$pv,
                    fc.lower = input$fc[1], fc.upper = input$fc[2],
                    score.lower = input$score[1], score.upper = input$score[2],
                    mfx = input$mfx,
                    geeq = input$geeq,
                    distance = input$distance)
    
    # Update the query string
    if(update) {
      query <- paste0('?genes=',
                      switch(min(2, length(genes)), genes, paste0('[', paste0(genes, collapse = ','), ']')),
                      switch((taxa == getOption('app.taxa')) + 1, paste0('&taxa=', taxa), ''),
                      switch((length(scope) == 1 && scope == getOption('app.ontology')) + 1,
                             paste0('&scope=', switch(min(2, length(scope)), scope,
                                                      paste0('[', paste0(scope, collapse = ','), ']'), ''))),
                      switch((options$pv == getOption('app.pv')) + 1, paste0('&pv=', options$pv), ''),
                      switch(isTRUE(all.equal(input$fc, c(getOption('app.fc_lower'), getOption('app.fc_upper')))) + 1,
                             paste0('&fc=[', paste0(input$fc, collapse = ','), ']'), ''),
                      switch((options$mfx == getOption('app.mfx')) + 1, paste0('&mfx=', options$mfx), ''),
                      switch((options$geeq == getOption('app.geeq')) + 1, paste0('&geeq=', options$geeq), ''),
                      switch((options$distance == getOption('app.distance_cutoff')) + 1, paste0('&distance=', options$distance), ''))
      updateQueryString(query, 'push')
    }
    
    if(is.null(taxa)) taxa <- getOption('app.taxa')
    if(is.null(scope)) scope <- getOption('app.ontology')
    
    # Clean numerics (interpreted as entrez IDs) and remove them from further processing.
    cleanGenes <- suppressWarnings(Filter(function(x) !is.na(as.integer(x)), genes))
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
        .[V1 %in% genes] # TODO may have multiple aliases for one entry
      if(nrow(aliases) > 0)
        cleanGenes <- c(cleanGenes, aliases[, entrez.ID])
    }
    
    # Done processing, handle the search.
    handleSearch(cleanGenes, CACHE.BACKGROUND[[taxa]], DATA.HOLDER[[taxa]], taxa, scope, options)
  }
}
