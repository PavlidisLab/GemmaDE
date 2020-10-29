library(pbapply)
source('dependencies.R')

#' Reverse Search
#'
#' @param contrasts 
#' @param taxa 
#' @param scope 
#' @param options 
#'
#' @return
#' @export
#'
#' @examples
revSearch <- function(contrasts,
                      taxa = getOption('app.taxa'), scope = getOption('app.ontology'),
                      options = getOption('app.all_options')) {
  lapply(contrasts, function(contrast) {
    candidates <<- DATA.HOLDER$human@gene.meta$entrez.ID[
      matrixStats::rowMins(DATA.HOLDER$human@data$adj.pv[, CACHE.BACKGROUND$human[
        cf.BaseLongUri %in% contrast | cf.ValLongUri %in% contrast, rsc.ID]
      ], na.rm = T) < options$pv
    ]
    
    searched <<- pblapply(candidates, function(gene) search(gene, taxa, options, verbose = F))
    enriched <<- pblapply(searched, function(search) enrich(search, taxa, scope, options, verbose = F))
    
    list(searched = searched, enriched = enriched)
  })
}

DATA.HOLDER$artificial <- NULL
DATA.HOLDER$mouse <- NULL
DATA.HOLDER$rat <- NULL

CACHE.BACKGROUND$artificial <- NULL
CACHE.BACKGROUND$mouse <- NULL
CACHE.BACKGROUND$rat <- NULL

saveRDS(revSearch(list(first = c('female', 'male'))), 'revSearch.rds')
