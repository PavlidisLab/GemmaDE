library(data.table)
library(jsonlite)
library(dplyr)

isDataLoaded <- function() {
  exists('ONTOLOGIES') & exists('ONTOLOGIES.DEFS') & exists('DATA.HOLDER') & exists('BLACKLISST')
}

# Load the data into the global environment
if(!isDataLoaded()) {
  setClass('EData', representation(taxon = 'character', data = 'list',
                                   experiment.meta = 'data.table', gene.meta = 'data.table'))
  
  TAXA <<- list(`H. sapiens` = 'human')#, `M. musculus` = 'mouse', `R. norvegicus` = 'rat')
  BLACKLIST <<- read.csv('data/blacklist.csv', header = F)$V1
  
  ONTOLOGIES <<- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV')
  ONTOLOGIES.DEFS <<- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV')
  
  DATA.HOLDER <<- lapply(TAXA, function(taxon) {
    load(paste0('/home/nlim/MDE/RScripts/DataFreeze/Packaged/Current/', taxon, '.RDAT.XZ'))
    
    dataHolder$adj.pv <- apply(dataHolder$pv, 2, function(pv)
      p.adjust(pv, method = 'BH', n = length(Filter(Negate(is.nan), pv))))
    
    new('EData', taxon = taxon, data = dataHolder,
        experiment.meta = metaData %>% as.data.table, gene.meta = metaGene %>% as.data.table)
  })
  
  names(DATA.HOLDER) <<- TAXA
}


#' JSONify
#' 
#' Assuming the input string is entirely unquoted, make it suitable for JSON.
jsonify <- function(str) {
  if(is.null(str) || length(grep('(\\[|\\])', str, value = F)) == 0) return(str)
  
  gsub('\\[', '["', gsub('\\]', '"]', gsub(',', '","', gsub(', ', ',', str)))) %>% parse_json(simplifyVector = T)
}
