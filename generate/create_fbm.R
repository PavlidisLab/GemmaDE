print('file based storage')

devtools::load_all()
library(gemma.R)
library(magrittr)
library(dplyr)

data.holder <- readRDS(file.path(RAWDIR,'DATA.HOLDER_fixed_terms.rds'))

if('matrix' %in% class(data.holder[[1]]@data$adj.pv)) {
  for(taxon in names(data.holder)) {
    message(paste0('Converting in-memory matrices for ', taxon, ' to file-backed...'))
    
    file.remove(list.files(file.path(DATADIR, 'fbm', taxon), full.names = T))
    
    dir.create(file.path(DATADIR,'fbm',taxon),recursive = TRUE, showWarnings =FALSE)
    
    dimnames(data.holder[[taxon]]@data$zscore) %>% 
      rev() %>%
      saveRDS(file.path(DATADIR, 'fbm', taxon, 'z.dimnames.rds'))
    
    data.holder[[taxon]]@data$zscore <- bigstatsr::as_FBM(data.holder[[taxon]]@data$zscore %>% t,
                                                          backingfile = file.path(DATADIR, 'fbm', taxon, 'zscores'),
                                                          is_read_only = T)$save()
    
    dimnames(data.holder[[taxon]]@data$adj.pv) %>% 
      rev() %>%
      saveRDS(file.path(DATADIR, 'fbm', taxon, 'p.dimnames.rds'))
    
    data.holder[[taxon]]@data$adj.pv <- bigstatsr::as_FBM(data.holder[[taxon]]@data$adj.pv %>% t,
                                                          backingfile = file.path(DATADIR, 'fbm', taxon, 'adjpvs'),
                                                          is_read_only = T)$save()
  }
}

dimnames.FBM <- function(object, ...) {
  attr(object, '.dimnames')
}

for (taxon in names(data.holder)) {
  data.holder[[taxon]]@data$zscore <- NULL
  data.holder[[taxon]]@data$adj.pv <- NULL
}
saveRDS(data.holder, paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/'))
