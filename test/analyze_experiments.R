library(parallel)
library(data.table)
library(dplyr)

set.seed(18232)
options('mc.cores' = 5)

N.replicates <- 1000
N.sample <- 10

lapply(c('OBI', 'UBERON', 'DO', 'HP', 'MP', 'CL'), function(ontology) {
  mclapply(1:N.replicates, function(n) {
    experiments <- search(sample(DATA.HOLDER$human@gene.meta$entrez.ID, N.sample))
    
    if(!is.null(experiments)) {
      mExp <- experiments$score
      names(mExp) <- rownames(experiments)
      mExp
    }
  }) %>% saveRDS(paste0('data/bootstrap/experiment.', ontology, '.rds'))
})
