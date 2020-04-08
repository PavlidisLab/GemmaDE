library(data.table)
library(jsonlite)
library(dplyr)
library(matrixStats)

isDataLoaded <- function() {
  exists('ONTOLOGIES') & exists('ONTOLOGIES.DEFS') & exists('DATA.HOLDER') &
    exists('BLACKLIST') & exists('CACHE.BACKGROUND')
}

# Load the data into the global environment
if(!isDataLoaded()) {
  setClass('EData', representation(taxon = 'character', data = 'list',
                                   experiment.meta = 'data.table', gene.meta = 'data.table'))
  
  TAXA <- list(`H. sapiens` = 'human')
  BLACKLIST <- c()
  
  CACHE.BACKGROUND <- list()
  
  DEFAULT_OPTIONS <- list(pv = 0.05,
                          fc.lower = 0, fc.upper = 10000,
                          mfx = T, filterSame = T)

  ONTOLOGIES <- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV')
  ONTOLOGIES.DEFS <- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV')
  
  ONTOLOGIES$ChildNode <- NULL
  ONTOLOGIES$ParentNode <- NULL
  ONTOLOGIES$ChildNode_Long <- ONTOLOGIES$ChildNode_Long %>% as.factor
  ONTOLOGIES$ParentNode_Long <- ONTOLOGIES$ParentNode_Long %>% as.factor
  ONTOLOGIES$RelationType <- ONTOLOGIES$RelationType %>% as.factor
  ONTOLOGIES$OntologyScope <- ONTOLOGIES$OntologyScope %>% as.factor
  
  ONTOLOGIES.DEFS$Node <- NULL
  ONTOLOGIES.DEFS$Node_Long <- ONTOLOGIES.DEFS$Node_Long %>% as.factor
  ONTOLOGIES.DEFS$OntologyScope <- ONTOLOGIES.DEFS$OntologyScope %>% as.factor
  
  DATA.HOLDER <- readRDS('data/human-lite.rds')
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
  if(length(entry) == 0) return(NA)
  
  if(is.na(entry) || !grepl('; ', entry, fixed = T)) return(entry)
  
  cleaned <- strsplit(entry %>% as.character, '; ', fixed = T) %>% unlist
  ifelse('NA' == cleaned, NA, cleaned)
}
