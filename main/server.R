library(shiny)
library(shinyjs)
library(htmlwidgets)
library(data.table)
library(dplyr)
library(shinycssloaders)

MAX_PROGRESS_STEPS <- 5

fineProgress <- function(session, steps) {
  if(!is.null(session))
    session$userData$progress.bar.steps <- steps
}

advanceProgress <- function(session, detail, fine = F) {
  if(!is.null(session))
    setProgress(session, session$userData$progress + 1, detail, fine)
}

setProgress <- function(session, progress, detail, fine = F) {
  if(!fine)
    session$userData$progress <- progress
  
  if(fine)
    session$userData$progress.bar$inc(amount = 1 / session$userData$progress.bar.steps)
  else if(progress == 0) {
    session$userData$progress.bar <- shiny::Progress$new(min = 0, max = MAX_PROGRESS_STEPS)
    session$userData$progress.bar$set(message = 'Searching...', detail = detail)
  } else {
    session$userData$progress.bar$set(value = progress, detail = detail)
  
    if(progress >= MAX_PROGRESS_STEPS)
      session$userData$progress.bar$close()
  }
}

#' Server
#'
#' Query strings are supported as follows:
#' ?genes=[list of genes]&taxa=[human|mouse|rat]&scope=[list of ontologies]
server <- function(input, output, session) {
  # On connect, observe the query string. Search if any genes are specified
  observeEvent(input$LOAD, {
    query <- getQueryString(session)
    
    genes <- query$genes %>% jsonify
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    session$sendCustomMessage('querySet', genes)
    updateSelectizeInput(session, 'taxa', selected = query$taxa)
    updateCheckboxInput(session, 'mfx', value = ifelse(is.null(query$mfx), DEFAULT_OPTIONS$mfx, query$mfx))
    updateNumericInput(session, 'pv', value = ifelse(is.null(query$pv), DEFAULT_OPTIONS$pv, query$pv))
    updateSliderInput(session, 'fc', value = ifelse(is.null(query$fc), c(DEFAULT_OPTIONS$fc.lower, DEFAULT_OPTIONS$fc.upper), query$fc %>% jsonify %>% as.numeric))
    
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
    updateSelectizeInput(session, 'scope', selected = 'DO')
    updateCheckboxInput(session, 'mfx', value = DEFAULT_OPTIONS$mfx)
    updateCheckboxInput(session, 'filter', value = DEFAULT_OPTIONS$filterSame)
    updateNumericInput(session, 'pv', value = DEFAULT_OPTIONS$pv)
    updateSliderInput(session, 'fc', value = c(DEFAULT_OPTIONS$fc.lower, DEFAULT_OPTIONS$fc.upper))
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
  #' @param taxa The taxon ([human|mouse|rat])
  #' @param options The search options
  #' @param scope The ontology scope
  #' @param options The search options
  handleSearch <- function(genes, taxa, scope, options = DEFAULT_OPTIONS) {
    experiments <- search(genes, taxa, scope, options, session)
    
    if(is.null(experiments)) {
      setProgress(session, MAX_PROGRESS_STEPS, '')
      output$results_header <- renderUI({
        generateResultsHeader(HTML('<h2 data-toggle="tooltip" data-placement="top" title="Relax thresholds or modify gene set.">Invalid search.</h2>'))
      })
      output$results <- NULL
    } else {
      conditions <- enrich(experiments, taxa, scope, session)
      
      if(is.null(conditions)) {
        setProgress(session, MAX_PROGRESS_STEPS, '')
        output$results_header <- renderUI({ generateResultsHeader('No conditions found in scope.') })
        output$results <- NULL
      } else {
        # Filter equal scoring children
        filterResults <- function(conditions) {
          do.call(rbind, lapply(1:(nrow(conditions) - 1), function(indx) {
            if(conditions$N.ranked[indx] != conditions$N.ranked[indx + 1])
              conditions[indx, ]
          })) %>% as.data.table
        }
        if(options$filter && nrow(conditions) > 1)
          conditions <- filterResults(conditions)
        
        # Compute experiments in scope
        mCache <- do.call(rbind, lapply(scope, function(id) CACHE.BACKGROUND[[taxa]][[id]])) %>%
          .[rsc.ID %in% rownames(experiments) & Definition %in% conditions$Definition]
        eDF <- experiments %>% as.data.frame
        rownames(eDF) <- rownames(experiments)
        colnames(eDF) <- colnames(experiments)
        
        # Generate the results header
        output$results_header <- renderUI({
          generateResultsHeader(HTML(paste0('<h2>Found ',
                                            '<span data-toggle="tooltip" data-placement="top" title="',
                                            length(mCache[, unique(rsc.ID)]), ' in scope">',
                                            nrow(experiments), '</span> experiments differentially expressing ',
                                            '<span data-toggle="tooltip" data-placement="top" title="',
                                            paste0(colnames(experiments)[1:(ncol(experiments) - 1)], collapse = ', '), '">',
                                            (ncol(experiments) - 1), ' genes</span>')))
        })
        
        # Merge in gene scores
        geneScores <- data.table(def = unique(conditions$Definition))
        
        for(indx in 1:(ncol(experiments) - 1)) {
          geneScores[, Evidence := paste0(DATA.HOLDER[[taxa]]@experiment.meta[rsc.ID %in% mCache[Definition == def, rsc.ID], unique(ee.ID)], collapse = ','), def]
          geneScores[, colnames(experiments)[indx] := paste0(eDF[mCache[Definition == def, rsc.ID], indx], collapse = ','), def]
          # geneScores[, colnames(experiments)[indx] := mean(eDF[mCache[Definition == def, rsc.ID], indx], na.rm = T), def]
        }
        
        # Rename things and add the ES
        conditions <- merge(conditions, geneScores, by.x = 'Definition', by.y = 'def') %>% setorder(pv.chisq) %>%
          setnames(c('chisq', 'pv.chisq', 'pv.fisher'), c('χ2', 'P-value (χ2)', 'P-value (Fisher)')) %>%
          .[, `E` := (N.ranked * (N.bg - N.y)) / ((N.tags - N.ranked) * N.y)]
        
        output$plot <- renderPlotly({ generateResultsPlot(taxa, scope, experiments, conditions, options, input) })
        output$results <- generateResults(taxa, scope, experiments, conditions, options)
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
    setProgress(session, 0, 'Validating input')
    
    options <- list(n.display = input$`top-n`,
                    pv = input$pv,
                    fc.lower = input$fc[1], fc.upper = input$fc[2],
                    score.lower = input$score[1], score.upper = input$score[2],
                    mfx = input$mfx,
                    filter = input$filter)
    
    # Update the query string
    if(update) {
      query <- paste0('?genes=',
                      switch(min(2, length(genes)), genes, paste0('[', paste0(genes, collapse = ','), ']')),
                      switch((taxa == TAXA[[1]]) + 1, paste0('&taxa=', taxa), ''),
                      switch((length(scope) == 1 && scope == 'DO') + 1,
                             paste0('&scope=', switch(min(2, length(scope)), scope,
                                                      paste0('[', paste0(scope, collapse = ','), ']'), ''))),
                      switch((options$pv == DEFAULT_OPTIONS$pv) + 1, paste0('&pv=', options$pv), ''),
                      switch(isTRUE(all.equal(input$fc, c(DEFAULT_OPTIONS$fc.lower, DEFAULT_OPTIONS$fc.upper))) + 1,
                             paste0('&fc=[', paste0(input$fc, collapse = ','), ']'), ''),
                      switch((options$mfx == DEFAULT_OPTIONS$mfx) + 1, paste0('&mfx=', options$mfx), ''),
                      switch((options$filter == DEFAULT_OPTIONS$filterSame) + 1, paste0('&filter=', options$filter), ''))
      updateQueryString(query, 'push')
    }
    
    if(is.null(taxa)) taxa <- 'human'
    if(is.null(scope)) scope <- 'DO'
    
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
        .[V1 %in% genes]
      if(nrow(aliases) > 0)
        cleanGenes <- c(cleanGenes, aliases[, entrez.ID])
    }
    
    # Done processing, handle the search.
    handleSearch(cleanGenes, taxa, scope, options)
  }
}
