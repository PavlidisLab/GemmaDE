DATADIR <- '/cosmos/data/project-data/GemmaDE'
FREEZEDIR <- '/cosmos/data/project-data/GemmaDE/gemma_freeze'

source(here::here("main/requirements.R"))
source(here::here("main/dependencies.R"))


#* DE Search
#* 
#* @param genes
#* @param taxon
#* @param max_dist
#* @get /de_search
de_search = function(genes,
                     taxon,
                     max_dist = 1.5){
  print(genes)
  print(taxon)
  print(max_dist)
}
