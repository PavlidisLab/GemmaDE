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
  
  TAXA <- list(`H. sapiens` = 'human')#, `M. musculus` = 'mouse', `R. norvegicus` = 'rat')
  BLACKLIST <- read.csv('data/blacklist.csv', header = F)$V1
  SUMMARY <- readRDS('data/bootstrap/bootstrap.rds')
  
  CACHE.BACKGROUND <- list()
  
  DEFAULT_OPTIONS <- list(pv = 0.05,
                          fc.lower = 0, fc.upper = 10,
                          mfx = T, filterSame = T)
  
  SUMMARY <- lapply(c('CL', 'DO', 'HP', 'MP', 'OBI', 'UBERON'), function(ontology) {
    do.call(rbind, lapply(readRDS(paste0('data/bootstrap/experiment.', ontology, '.rds')), function(entry) {
      data.table(rsc.ID = names(entry), score = entry)
    }))
  })
  names(SUMMARY) <- c('CL', 'DO', 'HP', 'MP', 'OBI', 'UBERON')
  
  # TODO
  #CATEGORIES <- list(DO = c('disease', 'treatment', 'dose', 'disease staging', 'clinical history', 'medical procedure'),
  #                    CLO = c('cell type', 'cell line', 'genotype', 'strain'),
  #                    CL = c('cell type', 'genotype', 'strain', 'biological process'),
  #                    CHEBI = c('treatment', 'clinical history', 'disease staging', 'growth condition', 'molecular entity', 'dose', 'medical procedure'),
  #                    GO = c('genotype', 'phenotype', 'biological process'),
  #                    EFO = c('treatment', 'timepoint', 'cell type', 'age', 'developmental stage', 'generation', 'block', 'temperature', 'diet', 'environmental stress', 'collection of material'),
  #                    HP = 'phenotype',
  #                    MP = 'phenotype',
  #                    OBI = c('treatment', 'disease staging', 'growth condition', 'collection of material', 'dose'),
  #                    UBERON = c('cell type', 'phenotype', 'collection of material', 'molecular entity', 'medical procedure', 'cell line'))
  
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
