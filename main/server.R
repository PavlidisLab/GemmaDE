library(shiny)
library(shinyjs)
library(htmlwidgets)
library(data.table)
library(dplyr)

MAX_PROGRESS_STEPS <- 7

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
  } else 
    session$userData$progress.bar$set(value = progress, detail = detail)
  
  if(progress >= MAX_PROGRESS_STEPS)
    session$userData$progress.bar$close()
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
    updateCheckboxInput(session, 'top-n', value = ifelse(is.null(query$N), DEFAULT_OPTIONS$n.display, query$N))
    updateCheckboxInput(session, 'mfx', value = ifelse(is.null(query$mfx), DEFAULT_OPTIONS$mfx, query$mfx))
    updateNumericInput(session, 'pv', value = ifelse(is.null(query$pv), DEFAULT_OPTIONS$pv, query$pv))
    updateSliderInput(session, 'fc', value = ifelse(is.null(query$fc), c(DEFAULT_OPTIONS$fc.lower, DEFAULT_OPTIONS$fc.upper), query$fc %>% jsonify %>% as.numeric))
    updateSliderInput(session, 'score', value = ifelse(is.null(query$score), c(DEFAULT_OPTIONS$score.lower, DEFAULT_OPTIONS$score.upper), query$score %>% jsonify %>% as.numeric))
    
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
    updateNumericInput(session, 'top-n', value = DEFAULT_OPTIONS$n.display)
    updateNumericInput(session, 'pv', value = DEFAULT_OPTIONS$pv)
    updateSliderInput(session, 'fc', value = c(DEFAULT_OPTIONS$fc.lower, DEFAULT_OPTIONS$fc.upper))
    updateSliderInput(session, 'score', value = c(DEFAULT_OPTIONS$score.lower, DEFAULT_OPTIONS$score.upper))
    session$sendCustomMessage('fileUpload', F)
  })
  
  # Search
  observeEvent(input$search, {
    genes <- input$genes
    if(!is.null(input$genes.csv))
      genes <- read.csv(input$genes.csv$datapath, header = F)$V1 %>% as.character
    
    searchGenes(genes)
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
    experiments <- search(genes, taxa, options, session)
    
    if(is.null(experiments)) {
      setProgress(session, MAX_PROGRESS_STEPS, '')
      output$results <- renderUI({ generateResultsHeader('Invalid search.') })
    } else {
      conditions <- enrich(experiments, taxa, scope, session)
      output$results <- generateResults(experiments, conditions, options)
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
                    mfx = input$mfx)
    
    # Update the query string
    if(update) {
      query <- paste0('?genes=',
                      switch(min(2, length(genes)), genes, paste0('[', paste0(genes, collapse = ','), ']')),
                      switch((taxa == TAXA[[1]]) + 1, paste0('&taxa=', taxa), ''),
                      switch((length(scope) == 1 && scope == 'DO') + 1,
                             paste0('&scope=', switch(min(2, length(scope)), scope,
                                                      paste0('[', paste0(scope, collapse = ','), ']'), ''))),
                      switch((options$n.display == DEFAULT_OPTIONS$n.display) + 1, paste0('&N=', options$n.display), ''),
                      switch((options$pv == DEFAULT_OPTIONS$pv) + 1, paste0('&pv=', options$pv), ''),
                      switch(isTRUE(all.equal(input$fc, c(DEFAULT_OPTIONS$fc.lower, DEFAULT_OPTIONS$fc.upper))) + 1,
                             paste0('&fc=[', paste0(input$fc, collapse = ','), ']'), ''),
                      switch(isTRUE(all.equal(input$score, c(DEFAULT_OPTIONS$score.lower, DEFAULT_OPTIONS$score.upper))) + 1,
                             paste0('&score=[', paste0(input$score, collapse = ','), ']'), ''),
                      switch((options$mfx == DEFAULT_OPTIONS$mfx) + 1, paste0('&mfx=', options$mfx), ''))
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
