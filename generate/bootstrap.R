source('/home/jsicherman/Thesis Work/requirements.R')

library(parallel)
options(mc.cores = 14)

ITERS <- 1000
BLOCK <- 500

OPTIONS <- c('artificial', 'human', 'mouse', 'rat')

for(x in OPTIONS) {
  if(!exists('DATA.HOLDER') || names(DATA.HOLDER) != x) {
    rm(DATA.HOLDER)
    source('/home/jsicherman/Thesis Work/dependencies.R')
    rm(NULLS)
    DATA.HOLDER[OPTIONS[-which(OPTIONS == x)]] <- NULL
  }
  
  message(paste0(Sys.time(), ' ... Starting ', x))
  if(F && file.exists(paste0('/space/scratch/jsicherman/Thesis Work/data/updated_nulls2/', x, '.rds'))) {
    message(paste0('File for ', x, ' already exists... Skipping.'))
  } else {
    mclapply(1:ITERS, function(j) {
      if(j %% (ITERS / 20) == 0)
        message(paste0(Sys.time(), ' ... ', x, ' ... ', round(100 * j / ITERS, 2), '%'))
      
      DATA.HOLDER[[x]]@gene.meta[sample(1:nrow(DATA.HOLDER[[x]]@gene.meta), BLOCK), entrez.ID] %>%
        search(getConfig(taxa = x)) %>%
        .[, c(1:500, 505)] %>% # TODO BLOCK
        melt(id.vars = 'rn') %>%
        .[, !'variable'] %>%
        .[, .(score.sum = sum(value, na.rm = T),
              score.sqsum = sum(value^2, na.rm = T)), rn]
    }) %>% rbindlist %>% {
      message('Saving...')
      
      .[, .(score.mean = sum(score.sum) / (BLOCK * ITERS),
            score.sd = sqrt(sum(score.sqsum) / (BLOCK * ITERS) - (sum(score.sum) / (BLOCK * ITERS))^2)), rn] %>%
        saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/updated_nulls2/', x, '.rds'))
      
      gc()
      NULL
    }
  }
}
