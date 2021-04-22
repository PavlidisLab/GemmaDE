#' Gene Evidence
#'
#' Gets a set of gene evidence from Gemma
#' 
#' @param genes The genes to query
#' @param taxa The taxa to query
geneEvidence <- async(function(genes, taxa = getConfig('taxa')$value) {
  prettyPrint <- function(evidence) {
    list(cf.ValLongUri = evidence$phenotypes[[1]]$valueUri,
         cf.CatLongUri = evidence$phenotypes[[1]]$categoryUri,
         symbol = evidence$geneOfficialSymbol,
         relationship = evidence$relationship,
         evidence.Code = evidence$evidenceCode,
         score.Name = evidence$scoreValueObject$scoreName,
         score.Value = evidence$scoreValueObject$scoreValue,
         score.Strength = evidence$scoreValueObject$strength,
         citation.Url = evidence$phenotypeAssPubVO[[1]]$citationValueObject$pubmedURL,
         citation.Name = evidence$phenotypeAssPubVO[[1]]$citationValueObject$citation)
  }
  
  parse <- function(json) {
    if(length(json$data) == 0) return(NULL)
    lapply(json$data[[1]]$evidence, prettyPrint)
  }
  
  if(length(taxa) > 1) {
    lapply(unique(genes[, taxon]), function(tax) {
      print(tax)
      lapply(genes[tax == taxon, entrez.ID], function(gene) {
        print(gene)
        http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/taxa/', tax, '/genes/', gene, '/evidence'))$then(function(response) parse(parse_json(rawToChar(response$content))))
      }) %>% { when_all(.list = .)$then(function(x) x %>% `names<-`(paste0(genes[tax == taxon, gene.realName], ' (', tax, ')'))) }
    }) %>% { when_all(.list = .)$then(function(x) unique(unlist(x, F))) }
  } else {
    lapply(genes, function(gene) {
      http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/taxa/', taxa, '/genes/', gene, '/evidence'))$then(function(response) parse(parse_json(rawToChar(response$content))))
    }) %>% { when_all(.list = .)$then(function(x) x %>% `names<-`(genes)) }
  }
})

#' Gene Expression
#' 
#' Gets gene expression data from Gemma. This is kind of gross and hacky. For each ee.ID provided, it will
#' 1. Request the gene expression data for the genes of interest
#' 2. Request sample information for each sample
#'   --> Subset only samples that have the baseline/contrasting factor
#'
#' @param ee.IDs Experiment IDs to search in
#' @param rsc.IDs The sub IDs to search for too
#' @param taxa The taxa scope
#' @param genes Genes to search for
#' @param keepNonSpecific, consolidate Options passed on to Gemma
geneExpression <- async(function(ee.IDs, rsc.IDs, taxa = getConfig('taxa')$value, genes, keepNonSpecific = T, consolidate = 'average') {
  parse <- function(content, json, meta) {
    mJson <- parse_json(rawToChar(content))
    
    # Tidy gene expression output
    prettyPrint <- function(gene) {
      if(length(gene$vectors) == 0) mData <- NULL
      else mData <- unlist(gene$vectors[[1]]$bioAssayExpressionLevels)
      
      extractSampleInfo <- function(sample) {
        # Factor value objects seems to be sufficient, should we merge with characteristic values?
        vals <- sapply(sample$sample$factorValueObjects, '[[', 'fvSummary') %>%
          `c`(sample$sample$characteristicValues %>% unname) %>% unique
        
        if((baseline <- any(meta[, cf.Baseline] %in% vals)) || any(meta[, cf.Val] %in% vals))
          list(name = sample$name,
               accession = ifelse(is.null(sample$accession$accession), 'N/A', sample$accession$accession),
               baseline = ifelse(baseline, meta[, data.table::first(cf.Baseline)], meta[, data.table::first(cf.Val)]))
      }
      
      if(is.null(mData)) NULL
      else {
        mMeta <- lapply(mJson$data, extractSampleInfo) %>% rbindlist %>% setorder(baseline)
        
        list(
          metadata = mMeta,
          expr = data.table(gene = gene$geneOfficialSymbol, # geneNcbiId
                            t(mData[mMeta[, name]] %>% `class<-`('numeric')))
        )
      }
    }
    
    list(ee.Name = as.character(data.table::first(meta[, ee.Name])),
         ee.Scale = as.character(data.table::first(meta[, ee.Scale])),
         geneData = lapply(json$data[[1]]$geneExpressionLevels, prettyPrint) %>% { Filter(Negate(is.null), .) }
    ) %>% {
      if(length(.$geneData) == 0) NULL
      else {
        # TODO this gross block could be made nicer when I'm less lazy
        expr <- rbindlist(lapply(.$geneData, '[[', 'expr'))
        expr.rn <- expr[, gene]
        expr <- expr[, !'gene']
        
        # Try to put everything on a log2 scale. Not sure about this.
        expr <- switch(.$ee.Scale,
                       LOG2 = expr,
                       LOG10 = expr / log10(2),
                       LN = expr / ln(2),
                       log2(expr + 1)
        )
        expr <- data.table(gene = expr.rn, expr)
        
        list(ee.Name = .$ee.Name,
             geneData = list(
               metadata = .$geneData[[1]]$metadata,
               expr = expr
             ))
      }
    }
  }
  
  # Get expression information for ee.ID by sending `length(ee.IDs)` async requests
  mLongData <- rbindlist(lapply(taxa, function(i) DATA.HOLDER[[i]]@experiment.meta[, .(ee.ID, rsc.ID, ee.Name = as.character(ee.Name), cf.Baseline, cf.Val, ee.Scale)]))
  
  lapply(ee.IDs, function(dataset) {
    # No guarantee that the contrast ordering we have is the same as in Gemma :(
    meta <- mLongData[ee.ID == dataset & rsc.ID %in% rsc.IDs] %>% {
      merge(.[, .(rsc.ID, ee.Name, ee.Scale)],
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
    
    if(length(taxa) > 1)
      mGenes <- unique(genes[, entrez.ID])
    else
      mGenes <- unique(genes)
    
    # TODO Look into these
    http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/datasets/', dataset,
                    '/expressions/genes/', paste0(mGenes, collapse = '%2C'),
                    '?keepNonSpecific=', ifelse(keepNonSpecific, 'true', 'false'),
                    '&consolidate=', consolidate))$then(function(response) {
                      content <- response$content
                      # Get the sample information for every expression set
                      http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/datasets/',
                                      data.table::first(meta[, ee.Name]), '/samples'))$
                        then(function(response2) parse(response2$content, parse_json(rawToChar(content)), meta))
                    })
  }) %>% { when_all(.list = .)$then(function(x) {
    # Once everything is done, pretty it up for return
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
