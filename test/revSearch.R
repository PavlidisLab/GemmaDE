library(pbapply)
source('dependencies.R')

revSearch <- function(data, cache, contrasts,
                      taxa = getOption('app.taxa'), scope = getOption('app.ontology'),
                      options = getOption('app.all_options'), session = NULL) {
  lapply(contrasts, function(contrast) {
    candidates <<- data@gene.meta$entrez.ID[
      matrixStats::rowMins(data@data$adj.pv[, cache[
        cf.BaseLongUri %in% contrast | cf.ValLongUri %in% contrast, rsc.ID]
      ], na.rm = T) < options$pv
    ]
    
    searched <<- pblapply(candidates, function(gene) search(data, gene, taxa, options, session))
    enriched <<- pblapply(searched, function(search) enrich(cache, search, taxa, scope, options, session))
    
    list(searched = searched, enriched = enriched)
  })
}

DATA.HOLDER$artificial <- NULL
DATA.HOLDER$mouse <- NULL
DATA.HOLDER$rat <- NULL

CACHE.BACKGROUND$artificial <- NULL
CACHE.BACKGROUND$mouse <- NULL
CACHE.BACKGROUND$rat <- NULL

saveRDS(revSearch(DATA.HOLDER$human, CACHE.BACKGROUND$human, list(first = c('male', 'female'))), 'revSearch.rds')
