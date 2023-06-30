print('file based storage')

devtools::load_all()
library(gemma.R)
library(magrittr)
library(dplyr)

data.holder <- readRDS(file.path(RAWDIR,'DATA.HOLDER_fixed_terms.rds'))

fbm_path = file.path(DATADIR,'data_fbm')

for(taxon in names(data.holder)) {
  create_fbm(data.holder[[taxon]]@data$fc,
             file.path(fbm_path,taxon,'fc'))
  
  create_fbm(data.holder[[taxon]]@data$zscore,
             file.path(fbm_path,taxon,'zscore'))
  
  
  create_fbm(data.holder[[taxon]]@data$adj.pv,
             file.path(fbm_path,taxon,'adj.pv'))
}


for (taxon in names(data.holder)) {
  data.holder[[taxon]]@data$zscore <- NULL
  data.holder[[taxon]]@data$adj.pv <- NULL
  data.holder[[taxon]]@data$fc <- NULL
}
saveRDS(data.holder, paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/'))
