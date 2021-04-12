source('/home/jsicherman/Thesis Work/requirements.R')

library(parallel)
options(mc.cores = 14)

ITERS <- 10000

OPTIONS <- c('human', 'artificial', 'mouse', 'rat')

for(x in OPTIONS) {
  if(!exists('DATA.HOLDER') || names(DATA.HOLDER) != x) {
    rm(DATA.HOLDER)
    source('/home/jsicherman/Thesis Work/dependencies.R')
    rm(NULLS.EXP)
    DATA.HOLDER[OPTIONS[-which(OPTIONS == x)]] <- NULL
  }
  
  i <- 1
  message(paste0(Sys.time(), ' ... Starting ', x, ' ', i))
  if(file.exists(paste0('/space/scratch/jsicherman/Thesis Work/data/nulls2/', x, '_', i, '.rds'))) {
    message(paste0('File for ', x, '_', i, ' already exists... Skipping.'))
  } else {
    mclapply(1:ITERS, function(j) {
      if(j %% (ITERS / 20) == 0)
        message(paste0(Sys.time(), ' ... ', x, ' ', i, ' ... ', round(100 * j / ITERS, 2), '%'))
      
      tmp <- DATA.HOLDER[[x]]@gene.meta[sample(1:nrow(DATA.HOLDER[[x]]@gene.meta), i), entrez.ID] %>%
        search(getConfig(taxa = x), inprod = F)
      
      if(is.null(tmp)) {
        message('No rankings on genes')
        NULL
      } else {
        tmp[, .(score, rn)]
      }
    }) %>% rbindlist %>% {
      message('Saving...')
      
      .[, .(score.mean = mean(score, na.rm = T),
            score.median = median(score, na.rm = T),
            score.sd = sd(score, na.rm = T)), rn] %>%
        saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/nulls2/', x, '_', i, '.rds'))
      
      gc()
      NULL
    }
  }
}
