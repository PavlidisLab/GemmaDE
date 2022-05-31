source('/home/jsicherman/Thesis Work/requirements.R')
source('/home/jsicherman/Thesis Work/dependencies.R')

library(parallel)

if(!exists('contrasts'))

  contrasts <- readRDS(paste(DATADIR, 'artificial/contrast_aff.rds', sep='/'))
uContrasts <- contrasts[, unique(contrast)]

DATA.HOLDER[c('human', 'mouse', 'rat')] <- NULL
NULLS[c('human', 'mouse', 'rat')] <- NULL

options(mc.cores = 30)
mclapply(1:250, function(n) {
  lapply(c(T, F), function(k) {
    lapply(c(1:10, 20, 50), function(i) {
      system(paste0('echo ', n, k, i))
      if(k)
        GOIs <- sample(1:nrow(DATA.HOLDER$artificial@gene.meta), i)
      else
        GOIs <- contrasts[contrast == sample(uContrasts, 1)] %>% setorder(-probability) %>%
          .[, head(entrez.ID, i)]
      
      contrasts[entrez.ID %in% GOIs] %>%
        .[, .(probability = sum(probability)), contrast] %>%
        merge(
          merge(DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
                CACHE.BACKGROUND$artificial[, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                by = 'rsc.ID') %>% .[, !'rsc.ID'] %>% unique,
          by.x = 'contrast', by.y = 'cf.ID', all = T) %>%
        unique %>%
        setorder(-probability) %>%
        .[, rank := rank(-probability, ties.method = 'min')] %>%
        .[, rank := as.integer(as.character(factor(rank, labels = 1:length(unique(rank)))))] -> tmp
      
      mSearch <- search(GOIs, getConfig(taxa = 'artificial'))
      if(is.null(mSearch))
        return(NULL)
      
      mSearch %>% enrich(getConfig(taxa = 'artificial', categories = getConfig('categories')$choices %>% unlist %>% unname)) %>%
        .[, I := .I] %>%
        merge(tmp[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, rank, probability)],
              by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all = T) %>%
        .[, .(stat = max(stat), I = min(I),
              rank = min(rank),
              probability = max(probability),
              score = max(score)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
        .[is.na(score), score := 0] %>%
        .[is.na(stat), stat := 0] %>%
        setorder(I) %>% {
          data.table(genes = i,
                     cor.stat = cor(.$stat, .$probability, method = 'spearman'),
                     cor.score = cor(.$score, .$probability, method = 'spearman'),
                     n = .[rank <= 5, sum(I <= 50, na.rm = T)],
                     random = k)
        }
    }) %>% rbindlist
  }) %>% rbindlist
}) %>% rbindlist -> new.rankings2

saveRDS(new.rankings2, paste(DATADIR, 'artificial/search_grid2.rds', sep='/'))
