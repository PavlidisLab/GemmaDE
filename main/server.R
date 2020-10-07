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
    output$results_header <- renderUI(generateResultsHeader(
      HTML('<div style="margin-bottom: 10px"><h2 class="loading" style="display: inline">EnriChing</h2></div>')))
    
    query <- getQueryString(session)
    
    genes <- query$genes %>% jsonify
    scope <- switch(is.null(query$scope) + 1, jsonify(query$scope), getOption('app.ontology'))
    
    updateSelectizeInput(session, 'genes', options = list(persist = F, create = T, createOnBlur = T))
    session$sendCustomMessage('querySet', genes)
    updateSelectizeInput(session, 'taxa', selected = ifelse(is.null(query$taxa), getOption('app.taxa'), query$taxa))
    updatePickerInput(session, 'scope', selected = scope)
    updateCheckboxInput(session, 'mfx', value = ifelse(is.null(query$mfx), getOption('app.mfx'), query$mfx))
    updateNumericInput(session, 'distance', value = ifelse(is.null(query$distance), getOption('app.distance_cutoff'), query$distance))
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
    updatePickerInput(session, 'scope', selected = getOption('app.ontology'))
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
                      selectInput('plot_type', 'Type', list('Heatmap', 'Scatterplot')),
                      selectInput('plot_data', 'Data', list('Gene Expression')),
                      pickerInput('plot_genes', 'Genes', list('Loading...'), multiple = T),
                      pickerInput('plot_conditions', 'Contrasts', list('Loading...'), multiple = T))
      )
    ), session)
    
    if(is.null(session$userData$plotData$processed)) {
      shinyjs::delay(1, {
        synchronise(session$userData$plotData$queued %>% {
          geneExpression(DATA.HOLDER[[session$userData$plotData$taxa]],
                         .[, unique(ee.ID)],
                         .[, unique(rsc.ID)],
                         session$userData$plotData$gene.ID)$then(function(x) {
                           rownames(x$expr) <- merge(data.table(ID = rownames(x$expr)),
                                                     data.table(ID = session$userData$plotData$gene.ID,
                                                                name = session$userData$plotData$gene.Name),
                                                     by = 'ID', sort = F, all.x = T)[, name]
                           
                           session$userData$plotData$processed <- x
                           session$userData$plotData$queued <- NULL
                           
                           output$plot <- renderPlotly({
                             generateResultsPlot(session$userData$plotData$gene.Name,
                                                 session$userData$plotData$conditions,
                                                 x, session$userData$plotData$options, input, session)
                           })
                         })
        })
      })
    } else {
      output$plot <- renderPlotly({
        generateResultsPlot(session$userData$plotData$gene.Name,
                            session$userData$plotData$conditions,
                            session$userData$plotData$processed,
                            session$userData$plotData$options, input, session)
      })
    }
  })
  
  geneEvidence <- async(function(genes, taxa = getOption('app.taxa')) {
    parse <- function(json) {
      lapply(json$data[[1]]$evidence, function(evidence) {
        list(cf.ValLongUri = evidence$phenotypes[[1]]$valueUri,
             cf.CatLongUri = evidence$phenotypes[[1]]$categoryUri,
             relationship = evidence$relationship,
             evidence.Code = evidence$evidenceCode,
             score.Name = evidence$scoreValueObject$scoreName,
             score.Value = evidence$scoreValueObject$scoreValue,
             score.Strength = evidence$scoreValueObject$strength,
             citation.Url = evidence$phenotypeAssPubVO[[1]]$citationValueObject$pubmedURL,
             citation.Name = evidence$phenotypeAssPubVO[[1]]$citationValueObject$citation)
      })
    }
    
    lapply(genes, function(gene) {
      http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/taxa/', taxa, '/genes/', gene, '/evidence'))$then(function(response)
        parse(parse_json(rawToChar(response$content))))
    }) %>% { when_all(.list = .)$then(function(x) x %>% `names<-`(genes)) }
  })
  
  #' Gene Expression
  #' 
  #' Gets gene expression data from Gemma. This is kind of gross and hacky. For each ee.ID provided, it will
  #' 1. Request the gene expression data for the genes of interest
  #' 2. Request sample information for each sample
  #'   --> Subset only samples that have the baseline/contrasting factor
  #'
  #' @param data 
  #' @param ee.IDs 
  #' @param rsc.IDs 
  #' @param genes 
  #' @param keepNonSpecific 
  #' @param consolidate 
  #'
  #' @return
  #' @export
  #'
  #' @examples
  geneExpression <- async(function(data, ee.IDs, rsc.IDs, genes, keepNonSpecific = F, consolidate = 'average') {
    parse <- function(json, meta) {
      http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/datasets/',
                      data.table::first(meta[, ee.Name]), '/samples'))$
        then(function(response) {
          list(ee.Name = as.character(data.table::first(meta[, ee.Name])),
               geneData = lapply(json$data[[1]]$geneExpressionLevels, function(gene) {
                 if(length(gene$vectors) == 0)
                   mData <- NULL
                 else
                  mData <- unlist(gene$vectors[[1]]$bioAssayExpressionLevels)
                 
                 if(is.null(mData))
                   NULL
                 else {
                   mJson <- parse_json(rawToChar(response$content))
                   mMeta <- lapply(mJson$data, function(sample) {
                     vals <- sapply(sample$sample$factorValueObjects, function(fv) fv$fvSummary)
                     if((baseline <- any(meta[, cf.Baseline] %in% vals)) || any(meta[, cf.Val] %in% vals))
                       list(name = sample$name,
                            accession = ifelse(is.null(sample$accession$accession), 'N/A', sample$accession$accession),
                            baseline = ifelse(baseline, meta[, data.table::first(cf.Baseline)], meta[, data.table::first(cf.Val)]))
                   }) %>% rbindlist %>% setorder(baseline)
                   
                   list(
                     metadata = mMeta,
                     expr = data.table(gene = gene$geneNcbiId,
                                       t(mData[mMeta[, name]] %>% `class<-`('numeric')))
                   )
                 }
               }) %>% { Filter(Negate(is.null), .) }
          ) %>% {
            if(length(.$geneData) == 0) NULL
            else
              list(ee.Name = .$ee.Name,
                   geneData = list(
                     metadata = .$geneData[[1]]$metadata,
                     expr = rbindlist(lapply(.$geneData, '[[', 'expr'))
                   ))
          }
        })
    }
    
    lapply(ee.IDs, function(dataset) {
      meta <- data@experiment.meta[ee.ID == dataset & rsc.ID %in% rsc.IDs,
                                   .(rsc.ID, ee.Name, cf.Baseline, cf.Val)] %>% {
          merge(.[, .(rsc.ID, ee.Name)],
                merge(
                  .[, .(cf.Baseline = parseListEntry(cf.Baseline) %>% {
                    apply(permutation(1:length(.)), 1, function(perm) paste0(.[perm], collapse = ', '))
                  }), rsc.ID],
                  .[, .(cf.Val = parseListEntry(cf.Val) %>% {
                    apply(permutation(1:length(.)), 1, function(perm) paste0(.[perm], collapse = ', '))
                  }), rsc.ID],
                  by = 'rsc.ID', sort = F, allow.cartesian = T
                ),
                by = 'rsc.ID', sort = F, allow.cartesian = T
          )
        }
      
      http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/datasets/', dataset,
                       '/expressions/genes/', paste0(genes, collapse = '%2C'),
                       '?keepNonSpecific=', ifelse(keepNonSpecific, 'true', 'false'),
                       '&consolidate=', consolidate))$then(function(response) {
                         parse(parse_json(rawToChar(response$content)), meta) })
    }) %>% { when_all(.list = .)$then(function(x) {
      Filter(Negate(is.null), x) %>% {
        list(metadata = rbindlist(lapply(., '[[', 'geneData') %>% lapply('[[', 'metadata')),
             expr = lapply(., '[[', 'geneData') %>% lapply('[[', 'expr') %>% {
               Reduce(function(...) merge.data.table(..., all = T, by = 'gene', sort = F), .)
             } %>% as.data.frame %>% `rownames<-`(.[, 'gene']) %>% .[, -1] %>% as.matrix
        )
      }
    })
    }
  })
  
  endFailure <- function(session) {
    setProgress(session, getOption('max.progress.steps'))
    output$results_header <- renderUI({
      generateResultsHeader(HTML('<h2 data-toggle="tooltip" data-placement="top" title="Relax thresholds or modify gene set.">Invalid search.</h2>'))
    })
    output$results <- NULL
  }
  
  endEmpty <- function(session) {
    setProgress(session, getOption('max.progress.steps'))
    output$results_header <- renderUI({ generateResultsHeader('No conditions found in scope.') })
    output$results <- NULL
  }
  
  endSuccess <- function(genes, experiments, conditions, cache, data,
                         scope = getOption('app.ontology'), options = getOption('app.all_options'), session = NULL) {
    # Generate the results header
    if(exists('output'))
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
    mCache <- getTags(cache, scope, rownames(experiments), options$distance) %>%
      .[cf.BaseLongUri %in% conditions[, unique(cf.BaseLongUri)] & cf.ValLongUri %in% conditions[, unique(cf.ValLongUri)]] %>%
      unique %>% merge(data@experiment.meta[rsc.ID %in% rownames(experiments), .(rsc.ID, ee.ID)], by = 'rsc.ID', sort = F, allow.cartesian = T) %>%
      .[, N := length(unique(ee.ID)), .(cf.BaseLongUri, cf.ValLongUri)] %>%
      merge(data.table(rsc.ID = rownames(experiments), experiments[, !'score']), by = 'rsc.ID')
    
    geneExpr <- mCache[paste0(cf.BaseLongUri, cf.ValLongUri) %in% conditions[pv.fisher < options$pv, paste0(cf.BaseLongUri, cf.ValLongUri)]]
    
    advanceProgress(session, 'Adding cross references')
    
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
    conditions <- merge(conditions, geneScores, by = c('cf.BaseLongUri', 'cf.ValLongUri'), sort = F, allow.cartesian = T) %>% setorder(pv.fisher) %>%
      setnames(c('pv.fisher', 'direction'), c('P-value', 'Direction'))
    
    advanceProgress(session, 'Finishing up')
    
    list(results = conditions, expression = geneExpr)
  }
  
  #' Handle Search
  #' 
  #' Call the processing function to handle searching and update the display.
  #'
  #' @param genes A character vector of entrez IDs
  #' @param cache A subsetted version of CACHE.BACKGROUND for the taxa scope.
  #' @param data A subsetted version of DATA.HOLDER for the taxa scope.
  #' @param taxa The taxon of interest.
  #' @param scope The ontology scope.
  #' @param options The search options
  handleSearch <- function(genes, cache, data, taxa = getOption('app.taxa'), scope = getOption('app.ontology'), options = getOption('app.all_options')) {
    synchronise({
      gene.evidence <- geneEvidence(genes, taxa)
      
      experiments <- search(data, genes, options, session)
      
      if(is.null(experiments))
        endFailure(session)
      else {
        conditions <- enrich(cache, experiments, scope, options, session)
        
        if(is.null(conditions))
          endEmpty(session)
        else {
          results <- endSuccess(genes, experiments, head(conditions, getOption('max.rows')), cache, data,
                                scope, options, session)
          
          output$results <- generateResults(data, experiments, results$results, options)
          
          session$userData$plotData <- list(
            taxa = taxa,
            gene.ID = genes,
            gene.Name = data@gene.meta[entrez.ID %in% genes, gene.Name],
            conditions = conditions[pv.fisher < options$pv, paste(cf.BaseLongUri, 'vs.', cf.ValLongUri)],
            options = options,
            queued = results$expression %>%
              merge(data@experiment.meta[, .(rsc.ID, ee.ID, ee.NumSamples)],
                    by = c('rsc.ID', 'ee.ID'), all.x = T, sort = F) %>% {
                      while((tmp <- .[sample(1:nrow(.))])[1, ee.NumSamples] > getOption('max.gemma')) {
                      }
                      tmp[cumsum(ee.NumSamples) <= getOption('max.gemma')]
                    }
          )
        }
      }
    })
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
    
    options <- list(pv = input$pv,
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
                      switch(identical(scope, getOption('app.ontology')) + 1,
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
