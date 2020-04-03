library(data.table)
library(jsonlite)
library(dplyr)

isDataLoaded <- function() {
  exists('ONTOLOGIES') & exists('ONTOLOGIES.DEFS') & exists('DATA.HOLDER') &
    exists('BLACKLIST') & exists('CACHE.BACKGROUND') & exists('SUMMARY')
}

# Load the data into the global environment
if(!isDataLoaded()) {
  setClass('EData', representation(taxon = 'character', data = 'list',
                                   experiment.meta = 'data.table', gene.meta = 'data.table'))
  
  TAXA <- list(`H. sapiens` = 'human')#, `M. musculus` = 'mouse', `R. norvegicus` = 'rat')
  BLACKLIST <- read.csv('data/blacklist.csv', header = F)$V1
  SUMMARY <- readRDS('data/bootstrap_summary.rds')
  
  CACHE.BACKGROUND <- list()
  
  DEFAULT_OPTIONS <- list(n.display = 50, pv = 0.05, fc.lower = 0, fc.upper = 10, score.lower = 0, score.upper = 500, mfx = T)
  
  # TODO
  CATEGORIES <- list(DO = c('disease', 'treatment', 'dose', 'disease staging', 'clinical history', 'medical procedure'),
                      CLO = c('cell type', 'cell line', 'genotype', 'strain'),
                      CL = c('cell type', 'genotype', 'strain', 'biological process'),
                      CHEBI = c('treatment', 'clinical history', 'disease staging', 'growth condition', 'molecular entity', 'dose', 'medical procedure'),
                      GO = c('genotype', 'phenotype', 'biological process'),
                      EFO = c('treatment', 'timepoint', 'cell type', 'age', 'developmental stage', 'generation', 'block', 'temperature', 'diet', 'environmental stress', 'collection of material'),
                      HP = 'phenotype',
                      MP = 'phenotype',
                      OBI = c('treatment', 'disease staging', 'growth condition', 'collection of material', 'dose'),
                      UBERON = c('cell type', 'phenotype', 'collection of material', 'molecular entity', 'medical procedure', 'cell line'))
  
  ONTOLOGIES <- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV')
  ONTOLOGIES.DEFS <- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV')
  
  DATA.HOLDER <- lapply(TAXA, function(taxon) {
    load(paste0('/home/nlim/MDE/RScripts/DataFreeze/Packaged/Current/', taxon, '.RDAT.XZ'))
    
    dataHolder$adj.pv <- apply(dataHolder$pv, 2, function(pv)
      p.adjust(pv, method = 'BH', n = length(Filter(Negate(is.nan), pv))))
    
    metaData$n.DE <- apply(dataHolder$adj.pv, 2, function(pv) sum(pv < 0.05, na.rm = T))
    
    new('EData', taxon = taxon, data = dataHolder,
        experiment.meta = metaData %>% as.data.table, gene.meta = metaGene %>% as.data.table)
  })
  
  names(DATA.HOLDER) <- unlist(TAXA)
}

#' JSONify
#' 
#' Assuming the input string is entirely unquoted, make it suitable for JSON.
jsonify <- function(str) {
  if(is.null(str) || length(grep('(\\[|\\])', str, value = F)) == 0) return(str)
  
  gsub('\\[', '["', gsub('\\]', '"]', gsub(',', '","', gsub(', ', ',', str)))) %>% parse_json(simplifyVector = T)
}

#' Parse List Entry
#' 
#' Given an entry (that may or may not be NA) of semicolon delimited (; ) list (in which, any may be NA),
#' return a vector of those values, coercing NA strings into true NAs.
#'
#' @param entry NA or a character of semicolon delimited (; ) values.
#' @param withKey A key to associate with this
#'
#' @return A character vector of values in the input character. If @param withKey is non-null, a data table with
#' two columns, one for each key and entry.
parseListEntry <- function(entry, withKey = NULL) {
  if(is.na(entry) || !grepl('; ', entry, fixed = T)) return(entry)
  
  cleaned <- strsplit(entry, '; ', fixed = T) %>% unlist
  ifelse('NA' == cleaned, NA, cleaned)
}
