source('/home/jsicherman/Thesis Work/requirements.R')

library(parallel)
options(mc.cores = 16)

OPTIONS <- c('human', 'artificial', 'mouse', 'rat')

for(x in OPTIONS) {
  if(!exists('DATA.HOLDER') || names(DATA.HOLDER) != x) {
    rm(DATA.HOLDER, NULLS.EXP)
    source('/home/jsicherman/Thesis Work/dependencies.R')
    DATA.HOLDER[OPTIONS[-which(OPTIONS == x)]] <- NULL
    NULLS.EXP[OPTIONS[-which(OPTIONS == x)]] <- NULL
  }
  
  mRange <- 1:length(DATA.HOLDER[[x]]@gene.meta$entrez.ID) %>% split(ceiling(seq_along(.[]) / getOption('chunk.size', 200)))
  mclapply(names(mRange), function(mID) {
    lapply(1:length(mRange[[mID]]), function(indx) {
      i <- mRange[[mID]][indx]
      if(indx %% 25 == 0)
        message(paste0(x, ' - ', mID, ': ', (100 * indx / length(mRange[[mID]])), '%'))
      
      tmp <- search(DATA.HOLDER[[x]]@gene.meta$entrez.ID[i], getConfig(taxa = x), F)
      if(is.null(tmp)) data.table(I = i)
      else {
        tmp %>% enrich(getConfig(taxa = x)) %>%
        .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, distance, stat, I = i)]
      }
    }) %>% rbindlist(fill = T) %>%
      saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/singlegene/', x, '_', mID, '.rds'))
    rm(mBlock)
    gc()
    NULL
  })
}
