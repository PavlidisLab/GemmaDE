library(data.table)
library(jsonlite)
library(dplyr)
library(matrixStats)
library(igraph)

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

isDataLoaded <- function() {
  exists('ONTOLOGIES') & exists('ONTOLOGIES.DEFS') &
    exists('BLACKLIST') & exists('CACHE.BACKGROUND')
}

# Load the data into the global environment
if(!isDataLoaded()) {
  setClass('EData', representation(taxon = 'character', data = 'list',
                                   experiment.meta = 'data.table', gene.meta = 'data.table'))
  
  TAXA <- list(`H. sapiens` = 'human', `M. musculus` = 'mouse', `R. norvegicus` = 'rat', TEST = 'artificial')
  BLACKLIST <- read.csv('data/blacklist.csv', header = F)$V1
  
  CACHE.BACKGROUND <- list()
  
  DEFAULT_OPTIONS <- list(pv = 0.05,
                          fc.lower = 0, fc.upper = 10,
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
}

# Load the lite versions if they're already created.
if(!exists('DATA.HOLDER')) {
  if(file.exists('data/DATA.HOLDER.rds'))
    DATA.HOLDER <- readRDS('data/DATA.HOLDER.rds')
  else {
    DATA.HOLDER <- lapply(TAXA, function(taxon) {
      load(paste0('/home/nlim/MDE/RScripts/DataFreeze/Packaged/Current/', taxon, '.RDAT.XZ'))
      
      dataHolder$adj.pv <- apply(dataHolder$pv, 2, function(pv)
        p.adjust(pv, method = 'BH', n = length(Filter(Negate(is.nan), pv))))
      
      metaData$n.DE <- colSums2(dataHolder$adj.pv < 0.05, na.rm = T)
      metaData$mean.fc <- colMeans2(dataHolder$fc, na.rm = T)
      
      dataHolder$ts <- NULL
      dataHolder$pv <- NULL
      dataHolder$meanrank <- NULL
      dataHolder$meanval <- NULL
      dataHolder$baserank <- NULL
      dataHolder$baseval <- NULL
      
      metaData <- metaData %>% as.data.table %>% .[, .(rsc.ID, ee.ID, ee.Name, ee.Source, ee.NumSamples,
                                                       ee.TagLongUri, ad.Name, ad.Company, ad.Sequencing,
                                                       sf.Subset, sf.Cat, sf.CatLongUri, sf.ValLongUri,
                                                       cf.Cat, cf.CatLongUri, cf.ValLongUri, cf.BaseLongUri,
                                                       n.DE, mean.fc)]
      
      metaGene <- metaGene %>% as.data.table %>% .[, .(entrez.ID, gene.ID, ensembl.ID, gene.Name, alias.Name,
                                                       gene.Desc, mfx.Rank)]
      
      metaData$ee.ID <- metaData$ee.ID %>% as.factor
      metaData$ee.Name <- metaData$ee.Name %>% as.factor
      metaData$ee.Source <- metaData$ee.Source %>% as.factor
      metaData$ee.TagLongUri <- metaData$ee.TagLongUri %>% as.factor
      metaData$ad.Name <- metaData$ad.Name %>% as.factor
      metaData$ad.Company <- metaData$ad.Company %>% as.factor
      metaData$ad.Sequencing <- metaData$ad.Sequencing %>% as.factor
      metaData$sf.Cat <- metaData$sf.Cat %>% as.factor
      metaData$sf.CatLongUri <- metaData$sf.CatLongUri %>% as.factor
      metaData$sf.ValLongUri <- metaData$sf.ValLongUri %>% as.factor
      metaData$cf.Cat <- metaData$cf.Cat %>% as.factor
      metaData$cf.CatLongUri <- metaData$cf.CatLongUri %>% as.factor
      metaData$cf.ValLongUri <- metaData$cf.ValLongUri %>% as.factor
      metaData$cf.BaseLongUri <- metaData$cf.BaseLongUri %>% as.factor
      
      new('EData', taxon = taxon, data = dataHolder,
          experiment.meta = metaData, gene.meta = metaGene)
    })
    
    names(DATA.HOLDER) <- unlist(TAXA)
    
    saveRDS(DATA.HOLDER, 'data/DATA.HOLDER.rds')
  }
  
  # Pre-load all ontology expansions
  lapply(unique(ONTOLOGIES$OntologyScope), function(scope) {
    lapply(TAXA, function(taxa) {
      getTags(taxa = taxa, scope = scope)
    })
  })
}
