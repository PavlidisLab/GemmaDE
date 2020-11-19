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

isDataLoaded <- function() {
  exists('ONTOLOGIES') & exists('ONTOLOGIES.DEFS')
}

setClass('EData', representation(taxon = 'character', data = 'list',
                                 experiment.meta = 'data.table', gene.meta = 'data.table'))

# Load the data into the global environment
if(!isDataLoaded()) {
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
  if(file.exists('/space/scratch/jsicherman/Thesis Work/data/DATA.HOLDER.rds'))
    DATA.HOLDER <- readRDS('/space/scratch/jsicherman/Thesis Work/data/DATA.HOLDER.rds')
  else {
    DATA.HOLDER <- lapply(Filter(function(x) x != 'artificial', getOption('app.all_taxa')), function(taxon) {
      load(paste0('/home/nlim/MDE/RScripts/DataFreeze/Packaged/Current/', taxon, '.RDAT.XZ'))
      
      dataHolder$ts <- NULL
      dataHolder$meanrank <- NULL
      dataHolder$meanval <- NULL
      dataHolder$baserank <- NULL
      dataHolder$baseval <- NULL
      
      metaData <- metaData %>% as.data.table %>% .[, .(rsc.ID, ee.Troubled, ee.Public, ee.ID, ee.Name, ee.Source,
                                                       ee.NumSamples, ee.TagLongUri, ad.Name, ad.Company,
                                                       ad.Sequencing, sf.Subset, sf.Cat, sf.CatLongUri, sf.ValLongUri,
                                                       cf.Cat, cf.CatLongUri, cf.Val, cf.ValLongUri, cf.Baseline, cf.BaseLongUri)]
      
      # Clean up the data.
      
      # Add structured text entries that will be dealt with specially
      metaData[is.na(cf.BaseLongUri), cf.BaseLongUri := cf.Baseline]
      metaData[is.na(cf.ValLongUri), cf.ValLongUri := cf.Val]
      
      clean <- function(baselines, baseUris) {
        unlist(lapply(1:length(baselines), function(i) {
          tmp <- parseListEntry(baseUris[i])
          tmp[is.na(tmp)] <- parseListEntry(baselines[i])[is.na(tmp)]
          paste0(tmp, collapse = '; ')
        }))
      }
      
      metaData[grepl('NA', cf.BaseLongUri, fixed = T), cf.BaseLongUri := clean(cf.Baseline, cf.BaseLongUri)]
      metaData[grepl('NA', cf.ValLongUri, fixed = T), cf.ValLongUri := clean(cf.Val, cf.ValLongUri)]
      
      # If it was a free-text entry, it becomes NA here.
      metaData[grepl('NA', cf.BaseLongUri, fixed = T), cf.BaseLongUri := gsub('; $', '', gsub('NA(; )?', '', cf.BaseLongUri))]
      metaData[grepl('NA', cf.ValLongUri, fixed = T), cf.ValLongUri := gsub('; $', '', gsub('NA(; )?', '', cf.ValLongUri))]
      
      # After filtering NAs, we should make them real NAs
      metaData[cf.BaseLongUri == '', cf.BaseLongUri := NA]
      metaData[cf.ValLongUri == '', cf.ValLongUri := NA]
      
      # We don't want:
      # Troubled or private experiments
      # Experiments where one/both contrasts is/are unknown
      # Experiments when the contrasts are identical (seemingly dose-dependent?)
      bad.rscs <- metaData[ee.Troubled | !ee.Public |
                             is.na(cf.BaseLongUri) |
                             is.na(cf.ValLongUri) |
                             cf.BaseLongUri == cf.ValLongUri, rsc.ID]
      
      # Drop experiment samples that don't meet our needs
      dataHolder$fc <- dataHolder$fc[, !(colnames(dataHolder$fc) %in% bad.rscs)]
      dataHolder$pv <- dataHolder$pv[, !(colnames(dataHolder$pv) %in% bad.rscs)]
      metaData <- metaData[!(rsc.ID %in% bad.rscs)]
      
      dataHolder$adj.pv <- apply(dataHolder$pv, 2, function(pv)
        p.adjust(pv, method = 'BH', n = length(pv))) # Assuming NaN p-values should be 1.
      dataHolder$pv <- NULL
      
      metaData[, n.DE := colSums2(dataHolder$adj.pv < 0.05, na.rm = T)]
      metaData[, mean.fc := colMeans2(dataHolder$fc, na.rm = T)]
      
      metaData$ee.ID <- metaData$ee.ID %>% as.integer
      
      metaGene <- metaGene %>% as.data.table %>% .[, .(entrez.ID, gene.ID, ensembl.ID, gene.Name, alias.Name,
                                                       gene.Desc, mfx.Rank)]
      
      metaGene[, n.DE := rowSums2(dataHolder$adj.pv < 0.05, na.rm = T)]
      metaGene[, dist.Mean := rowMeans2(dataHolder$fc, na.rm = T)]
      metaGene[, dist.SD := Rfast::rowVars(dataHolder$fc, na.rm = T, std = T)]
      
      # Precompute z-scores and p-weighted z-scores. Need to maintain logFC and adj.pv for user-selected filtering
      dataHolder$zscore <- (dataHolder$fc - metaGene$dist.Mean) / metaGene$dist.SD
      
      dataHolder$pvz <- dataHolder$zscore %>% {
        tmp <- dataHolder$adj.pv
        tmp[is.na(tmp)] <- 1
        tmp[tmp < 1e-20] <- 1e-20
        abs(.) * -log(tmp, 100)
      }
      
      # Split into 500 ee.ID chunks (so the URI doesn't get too long) and fetch quality scores 
      # from Gemma for all experiments.
      metaData <- rbindlist(lapply(lapply(metaData$ee.ID %>% unique %>% split(ceiling(seq_along(.[]) / 500)),
                                               datasetInfo, memoized = T) %>% unlist(recursive = F), function(ee.ID) {
                                                 data.table(ee.ID = ee.ID$id,
                                                            ee.qScore = ee.ID$geeq$publicQualityScore,
                                                            ee.sScore = ee.ID$geeq$publicSuitabilityScore)
                                                 })) %>% merge(metaData, by = 'ee.ID', all.y = T)
      forgetGemmaMemoised()
      
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
    
    names(DATA.HOLDER) <- unlist(Filter(function(x) x != 'artificial', getOption('app.all_taxa')))
    
    saveRDS(DATA.HOLDER, '/space/scratch/jsicherman/Thesis Work/data/DATA.HOLDER.rds')
  }
  
  # Pre-load all ontology expansions
  if(file.exists('/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND.rds'))
    CACHE.BACKGROUND <- readRDS('/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND.rds')
  else {
    CACHE.BACKGROUND <- lapply(Filter(function(x) x %in% names(DATA.HOLDER), getOption('app.all_taxa')), precomputeTags)
    
    names(CACHE.BACKGROUND) <- Filter(function(x) x %in% names(DATA.HOLDER), getOption('app.all_taxa'))
    
    saveRDS(CACHE.BACKGROUND, '/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND.rds')
  }
}

rm(isDataLoaded)
