library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)
library(parallel)

options(app.name = 'EnriCh')
options(app.description = 'Pending')
options(app.tags = 'genomics,bioinformatics,genetics,transcriptomes,rnaseq,microarrays,biotechnology,medicine,biomedical,meta-analysis,statistics,search,open source,database,software')
options(app.author = 'Jordan Sicherman (jordan.sicherman@msl.ubc.ca)')
options(spinner.color = '#002145')
options(spinner.type = 6)

options(max.progress.steps = 7)
options(max.gemma = 500)
options(max.rows = 100)
options(app.pv = 0.05)
options(app.fc_lower = 0, app.fc_upper = 10)
options(app.distance_cutoff = 2.25)
options(app.mfx = T, app.geeq = T)
options(app.taxa = 'human',
        app.ontology = c('GO', 'CLO', 'CL', 'CHEBI', 'DO', 'EFO', 'TGEMO', 'HP', 'MP', 'OBI', 'UBERON'))
options(app.all_taxa = list(`H. sapiens` = 'human', `M. musculus` = 'mouse', `R. norvegicus` = 'rat', artificial = 'artificial'))

options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                               fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                               geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff'),
                               max.rows = getOption('max.rows')))

source('/home/jsicherman/Thesis Work/main/process.R')
source('/home/jsicherman/Thesis Work/main/load.R')

N_ITERS <- 20000
SCOPE <- getOption('app.ontology')

options(mc.cores = 2)

lapply(names(DATA.HOLDER), function(taxa) {
  N_GENES <- rnorm(N_ITERS, 20, 4) %>% round
  mclapply(1:N_ITERS, function(iteration) {
    message(paste0(taxa, ': ', round(100*iteration/N_ITERS, 2), '%'))
    s <- search(DATA.HOLDER[[taxa]]@gene.meta[, sample(entrez.ID, N_GENES[iteration])], taxa, verbose = F)
    list(search = s,
         enrich = enrich(s, taxa, SCOPE, verbose = F))
  }) -> bootstrap
  
  saveRDS(bootstrap, paste0('/space/scratch/jsicherman/Thesis Work/data/bootstrap_', taxa, '.rds'))
  rm(bootstrap)
})

q('no')
