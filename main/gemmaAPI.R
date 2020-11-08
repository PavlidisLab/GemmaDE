#' Gene Evidence
#'
#' Gets a set of gene evidence from Gemma
#' 
#' @param genes The genes to query
#' @param taxa The taxa to query
geneEvidence <- async(function(genes, taxa = getOption('app.taxa')) {
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
  
  parseInner <- function(response) {
    parse(parse_json(rawToChar(response$content)))
  }
  
  lapply(genes, function(gene) {
    http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/taxa/', taxa, '/genes/', gene, '/evidence'))$then(parseInner)
  }) %>% { when_all(.list = .)$then(function(x) x %>% `names<-`(genes)) }
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
geneExpression <- async(function(ee.IDs, rsc.IDs, taxa = getOption('app.taxa'), genes, keepNonSpecific = F, consolidate = 'average') {
  extractSampleInfo <- function(sample, meta) {
    vals <- sapply(sample$sample$factorValueObjects, function(fv) fv$fvSummary)
    if((baseline <- any(meta[, cf.Baseline] %in% vals)) || any(meta[, cf.Val] %in% vals))
      list(name = sample$name,
           accession = ifelse(is.null(sample$accession$accession), 'N/A', sample$accession$accession),
           baseline = ifelse(baseline, meta[, data.table::first(cf.Baseline)], meta[, data.table::first(cf.Val)]))
  }
  
  prettyPrintGene <- function(gene, mJson, meta) {
    if(length(gene$vectors) == 0)
      mData <- NULL
    else
      mData <- unlist(gene$vectors[[1]]$bioAssayExpressionLevels)
    
    if(is.null(mData))
      NULL
    else {
      mMeta <- lapply(mJson$data, extractSampleInfo, meta) %>% rbindlist %>% setorder(baseline)
      
      list(
        metadata = mMeta,
        expr = data.table(gene = gene$geneNcbiId,
                          t(mData[mMeta[, name]] %>% `class<-`('numeric')))
      )
    }
  }
  
  parseInner <- function(response, json, meta) {
    mJson <- parse_json(rawToChar(response$content))
    list(ee.Name = as.character(data.table::first(meta[, ee.Name])),
         geneData = lapply(json$data[[1]]$geneExpressionLevels, prettyPrintGene, mJson, meta) %>% { Filter(Negate(is.null), .) }
    ) %>% {
      if(length(.$geneData) == 0) NULL
      else
        list(ee.Name = .$ee.Name,
             geneData = list(
               metadata = .$geneData[[1]]$metadata,
               expr = rbindlist(lapply(.$geneData, '[[', 'expr'))
             ))
    }
  }
  
  parse <- function(json, meta) {
    http_get(paste0('https://gemma.msl.ubc.ca/rest/v2/datasets/',
                    data.table::first(meta[, ee.Name]), '/samples'))$
      then(function(response) parseInner(response, json, meta))
  }
  
  lapply(ee.IDs, function(dataset) {
    meta <- DATA.HOLDER[[taxa]]@experiment.meta[ee.ID == dataset & rsc.ID %in% rsc.IDs,
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
