
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
#'
#' @return A character vector of values in the input character.
parseListEntry <- function(entry) {
  if(length(entry) == 0) return(NA)
  
  if(is.na(entry) || !grepl('; ', entry, fixed = T)) return(entry)
  
  cleaned <- strsplit(entry %>% as.character, '; ', fixed = T) %>% unlist
  ifelse('NA' == cleaned, NA, cleaned)
}

fixOntoGenes <- function() {
  lapply(getConfig(key = 'taxa')$core, function(taxa) {
    matches <- do.call(rbind,
                       str_match_all(DATA.HOLDER[[taxa]]@experiment.meta[, c(as.character(cf.BaseLongUri),
                                                                             as.character(cf.ValLongUri))],
                                     'http://purl.org/commons/record/ncbi_gene/(\\d*)')) %>% unique
    
    matches[, 2] <- lapply(getConfig(key = 'taxa')$core, function(taxa) {
      DATA.HOLDER[[taxa]]@gene.meta[entrez.ID %in% matches[, 2],
                                    .(entrez.ID, gene.Name = paste0(gene.Name, ' [', taxa, ']'))]
    }) %>% rbindlist %>% unique(by = 'entrez.ID') %>% .[match(matches[, 2], entrez.ID), gene.Name]
    
    data.table(Node_Long = matches[, 1], Definition = matches[, 2], OntologyScope = 'TGEMO')
  }) %>% rbindlist %>% {
    rbind(ONTOLOGIES.DEFS[, .(Node_Long = as.character(Node_Long),
                              Definition = as.character(Definition),
                              OntologyScope = as.character(OntologyScope))], .) %>%
      .[, c('Node_Long', 'Definition', 'OntologyScope') := list(as.factor(Node_Long),
                                                                as.factor(Definition),
                                                                as.factor(OntologyScope))]
  }
}

loadDrugbank <- function() {
  if(file.exists('/space/scratch/jcastillo/Thesis Work/data/drugbank/drugbank.rds'))
    return(readRDS('/space/scratch/jcastillo/Thesis Work/data/drugbank/drugbank.rds'))
  
  dbank <- xmlParse('/space/scratch/jcastillo/Thesis Work/data/drugbank/full database.xml')
  droot <- xmlRoot(dbank)
  dsize <- xmlSize(droot)
  
  tmp <- lapply(1:dsize, function(i) {
    droot[[i]] %>% xmlToList %>%
      .[c('name', 'synonyms', 'categories', 'targets')] %>%
      rbind %>% as.data.table
  }) %>% rbindlist(fill = T)
  
  tmp[, I := .I]
  
  tmp[, name := unlist(name)]
  tmp <- tmp[, .(name,
                 synonym = unlist(synonyms, F) %>% sapply('[[', 'text') %>% unname %>% list,
                 category = unlist(categories, F) %>% sapply('[[', 'category') %>% unname %>% list,
                 target = unlist(targets, F) %>% sapply('[[', 'name') %>% unname %>% list), I] %>%
    .[, !'I']
  
  saveRDS(tmp, '/space/scratch/jcastillo/Thesis Work/data/drugbank/drugbank.rds')
  
  tmp
}

setClass('EData', representation(taxon = 'character', data = 'list',
                                 experiment.meta = 'data.table', gene.meta = 'data.table',
                                 go = 'data.table'))

# Load the data into the global environment
if(!exists('ONTOLOGIES') || !exists('ONTOLOGIES.DEFS')) {
  ONTOLOGIES <- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED.TSV')
  ONTOLOGIES.DEFS <- fread('/space/grp/nlim/CronGemmaDump/Ontology/Ontology_Dump_MERGED_DEF.TSV')
  
  ONTOLOGIES[, c('ChildNode', 'ParentNode')] <- NULL
  ONTOLOGIES.DEFS$Node <- NULL
  
  ONTOLOGIES$ChildNode_Long <- ONTOLOGIES$ChildNode_Long %>% as.factor
  ONTOLOGIES$ParentNode_Long <- ONTOLOGIES$ParentNode_Long %>% as.factor
  ONTOLOGIES$RelationType <- ONTOLOGIES$RelationType %>% as.factor
  ONTOLOGIES$OntologyScope <- ONTOLOGIES$OntologyScope %>% as.factor
  
  ONTOLOGIES.DEFS$Node_Long <- ONTOLOGIES.DEFS$Node_Long %>% as.factor
  ONTOLOGIES.DEFS$OntologyScope <- ONTOLOGIES.DEFS$OntologyScope %>% as.factor
}

.DATA_PATH <- '/space/scratch/jcastillo/Thesis Work/data/DATA.HOLDER.light.rds'

# Load the lite versions if they're already created.
if(!exists('DATA.HOLDER')) {
  if(file.exists(.DATA_PATH))
    DATA.HOLDER <- readRDS(.DATA_PATH)
  else {
    stop("Couldn't find DATA.HOLDER, run generate/fbm.R first.")
  }
  
  DATA.HOLDER$artificial <- NULL
  
  ONTOLOGIES.DEFS <- fixOntoGenes()
}

# Read existing FBMs
for(i in names(DATA.HOLDER)) {
  DATA.HOLDER[[i]]@data$zscore <- big_attach(paste0('/space/scratch/jcastillo/Thesis Work/data/fbm/', i, '/zscores.rds'))
  attr(DATA.HOLDER[[i]]@data$zscore, '.dimnames') <- readRDS(paste0('/space/scratch/jcastillo/Thesis Work/data/fbm/', i, '/z.dimnames.rds'))
  
  DATA.HOLDER[[i]]@data$adj.pv <- big_attach(paste0('/space/scratch/jcastillo/Thesis Work/data/fbm/', i, '/adjpvs.rds'))
  attr(DATA.HOLDER[[i]]@data$adj.pv, '.dimnames') <- readRDS(paste0('/space/scratch/jcastillo/Thesis Work/data/fbm/', i, '/p.dimnames.rds'))
}
rm(i)

gc()

dimnames.FBM <- function(object, ...) {
  attr(object, '.dimnames')
}

if(!exists('CACHE.BACKGROUND')) {
  # Pre-load all ontology expansions
  if(file.exists('/space/scratch/jcastillo/Thesis Work/data/CACHE.BACKGROUND.rds'))
    CACHE.BACKGROUND <- readRDS('/space/scratch/jcastillo/Thesis Work/data/CACHE.BACKGROUND.rds')
  else {
    CACHE.BACKGROUND <- lapply(names(DATA.HOLDER), precomputeTags) %>%
      setNames(names(DATA.HOLDER))
    
    saveRDS(CACHE.BACKGROUND, '/space/scratch/jcastillo/Thesis Work/data/CACHE.BACKGROUND.rds')
  }
}

if(!exists('NULLS')) {
  NULLS <- lapply(names(DATA.HOLDER), function(taxa) {
    tryCatch(readRDS(paste0('/space/scratch/jcastillo/Thesis Work/data/updated_nulls2/', taxa, '.rds')) %>%
      .[, .(rn, score.mean, score.sd)], error = function(e) NULL)
  }) %>% `names<-`(names(DATA.HOLDER)) %>% {
    Filter(Negate(is.null), .)
  }
}

if(!exists('DRUGBANK') && Sys.getenv('RSTUDIO') == '1') {
  DRUGBANK <- loadDrugbank()
}
