# this is my re-organization of things that are exclusively required for
# Jordan's thesis.
# - ogan

library(here)
library(gemma.R)
library(magrittr)
options('gemma.cache' = 'cache_in_memory')
# options('gemma.memoised' = FALSE)


gemma_user <- readLines('generate/auth')
gemma.R::set_gemma_user(gemma_user[1],gemma_user[2])
rm(gemma_user)

experiments <- lapply(c('human', 'mouse', 'rat'), function(taxon){
  print(taxon)
  poke_call <- gemma.R::get_taxon_datasets(taxon,limit = 1)
  
  seq(0,attributes(poke_call)$totalElements,100) %>% lapply(function(x){
    cat('=')
    out <- gemma.R::get_taxon_datasets(taxon,offset = x,limit = 100)
  }) %>% do.call(rbind,.) -> all_datasets
  
}) %>% do.call(rbind,.)


platforms <- experiments$experiment.ID %>% lapply(function(id){
  cat('=')
  p = gemma.R::get_dataset_platforms(id,memoised = TRUE)
  p$experiment.ID = id
  return(p)
}) %>% do.call(rbind,.)

platforms %>% dplyr::select(experiment.ID,platform.Type) %>% 
  unique %>% .[!.$experiment.ID %in% .$experiment.ID[duplicated(.$experiment.ID)],] %>% dim

experiments$platform.Type = platforms$platform.Type[match(experiments$experiment.ID,platforms$experiment.ID)]

saveRDS(experiments,file = 'report/experiments.rds')
