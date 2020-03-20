library(shiny)
library(htmlwidgets)
library(data.table)
library(dplyr)

#' Server
#'
#' Query strings are supported as follows:
#' ?genes=[list of genes]&taxa=[human|mouse|rat]&scope=all|[list of ontologies]
server <- function(input, output, session) {
  
  # On connect, observe the query string. Search if any genes are specified
  observeEvent(input$LOAD, {
    query <- getQueryString(session)
    
    genes <- query$genes %>% jsonify
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    session$sendCustomMessage('querySet', genes)
    updateSelectizeInput(session, 'taxa', selected = query$taxa)
    
    # Use DO as a default scope, if unspecified
    scope <- query$scope %>% jsonify
    if(is.null(scope)) scope <- 'DO'
    
    updateSelectizeInput(session, 'scope', selected = scope)
    
    # Search if any genes are specified. Specify parameters to make sure no reactives.
    if(!is.null(query$genes))
      searchGenes(query$genes, query$taxa, scope, update = F)
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
    session$sendCustomMessage('fileUpload', F)
  })
  
  # Search
  observeEvent(input$search, {
    genes <- input$genes
    if(!is.null(input$genes.csv))
      genes <- read.csv(input$genes.csv$datapath, header = F)$V1 %>% as.character
    
    searchGenes(genes)
  })
  
  searchGenes <- function(genes = input$genes,
                          taxa = input$taxa,
                          scope = input$scope,
                          update = T) {
    if(update) {
      query <- paste0('?genes=',
                      switch(min(2, length(genes)), genes, paste0('[', paste0(genes, collapse = ','), ']')),
                      switch((taxa == TAXA[[1]]) + 1, paste0('&taxa=', taxa), ''),
                      switch((length(scope) == 1 && scope == 'DO') + 1,
                             paste0('&scope=', switch(min(2, length(scope)), scope,
                                                      paste0('[', paste0(scope, collapse = ','), ']'), ''))))
      updateQueryString(query, 'push')
    }
  }
}
