library(pbapply)
source('../main/process.R')
source('../main/server.R')
source('../main/load.R')
source('../main/ui.R')

options(app.name = 'EnriCh')
options(app.description = 'Pending')
options(app.tags = 'genomics,bioinformatics,genetics,transcriptomes,rnaseq,microarrays,biotechnology,medicine,biomedical,meta-analysis,statistics,search,open source,database,software')
options(app.author = 'Jordan Sicherman (jordan.sicherman@msl.ubc.ca)')
options(spinner.color = '#002145')
options(spinner.type = 6)
options('future.globals.maxSize' = 5 * 1024^3)

options(max.progress.steps = 8)
options(max.rows = 500)
options(app.pv = 0.05)
options(app.fc_lower = 0, app.fc_upper = 10)
options(app.distance_cutoff = 2.25)
options(app.mfx = T, app.geeq = T)
options(app.taxa = 'human', app.ontology = 'DO')
options(app.all_taxa = list(`H. sapiens` = 'human', `M. musculus` = 'mouse', `R. norvegicus` = 'rat', artificial = 'artificial'))

options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                               fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                               geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff')))

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
