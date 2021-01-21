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
  if(file.exists('/space/scratch/jsicherman/Thesis Work/data/drugbank/drugbank.rds'))
    return(readRDS('/space/scratch/jsicherman/Thesis Work/data/drugbank/drugbank.rds'))
  
  dbank <- xmlParse('/space/scratch/jsicherman/Thesis Work/data/drugbank/full database.xml')
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
  
  saveRDS(tmp, '/space/scratch/jsicherman/Thesis Work/data/drugbank/drugbank.rds')
  
  tmp
}

setClass('EData', representation(taxon = 'character', data = 'list',
                                 experiment.meta = 'data.table', gene.meta = 'data.table',
                                 go = 'data.table'))

# Load the data into the global environment
if(!exists('ONTOLOGIES') || !exists('ONTOLOGIES.DEFS')) {
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
    DATA.HOLDER <- lapply(getConfig(key = 'taxa')$core, function(taxon) {
      load(paste0('/home/nlim/MDE/RScripts/DataFreeze/Packaged/Current/', taxon, '.RDAT.XZ'))
      
      dataHolder$ts <- NULL
      dataHolder$meanrank <- NULL
      dataHolder$baserank <- NULL
      dataHolder$meanval <- NULL
      dataHolder$baseval <- NULL
      
      metaData <- metaData %>% as.data.table %>% .[, .(rsc.ID, ee.Troubled, ee.Public, ee.ID, ee.Name, ee.Source, ee.Scale,
                                                       ee.NumSamples, ee.TagLongUri, ad.Name, ad.Company,
                                                       ad.Sequencing, sf.Subset, sf.Cat, sf.CatLongUri, sf.Val, sf.ValLongUri,
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
      
      # Experiments who have any rsc that are DE_Exclude or Include should be ignored
      bad.ees <- metaData[sf.ValLongUri == 'http://purl.obolibrary.org/obo/TGEMO_00014' |
                            sf.Val == 'DE_Exclude' |
                            sf.ValLongUri == 'http://purl.obolibrary.org/obo/TGEMO_00013' |
                            sf.Val == 'DE_Include', ee.ID]
      bad.rscs <- metaData[ee.ID %in% bad.ees, rsc.ID]
      
      # We don't want:
      # Troubled or private experiments
      # Experiments where one/both contrasts is/are unknown
      # Experiments when the contrasts are identical (seemingly dose-dependent?)
      bad.rscs <- c(bad.rscs, metaData[ee.Troubled | !ee.Public |
                                         is.na(cf.BaseLongUri) |
                                         is.na(cf.ValLongUri) |
                                         cf.BaseLongUri == cf.ValLongUri, rsc.ID]) %>% unique
      
      # Drop experiment samples that don't meet our needs
      dataHolder$fc <- dataHolder$fc[, !(colnames(dataHolder$fc) %in% bad.rscs)]
      dataHolder$pv <- dataHolder$pv[, !(colnames(dataHolder$pv) %in% bad.rscs)]
      metaData <- metaData[!(rsc.ID %in% bad.rscs)]
      
      dataHolder$fc[is.nan(dataHolder$fc)] <- NA
      dataHolder$pv[is.nan(dataHolder$pv)] <- NA
      
      dataHolder$adj.pv <- apply(dataHolder$pv, 2, function(pv)
        p.adjust(pv, method = 'BH', n = length(pv))) # TODO assuming NaN p-values should be 1, they may also be missing genes.
      dataHolder$pv <- NULL
      
      metaData[, n.DE := colSums2(dataHolder$adj.pv <= 0.05, na.rm = T)]
      metaData[, mean.fc := colMeans2(dataHolder$fc, na.rm = T)]
      metaData[, n.detect := nrow(dataHolder$fc) - colSums2(dataHolder$fc %>% is.na)]
      
      metaData$ee.ID <- metaData$ee.ID %>% as.integer
      
      metaGene <- metaGene %>% as.data.table %>% .[, .(entrez.ID, gene.ID, ensembl.ID, gene.Name,
                                                       alias.Name, gene.Desc, mfx.Rank)]
      
      metaGene[, n.DE := rowSums2(dataHolder$adj.pv <= 0.05, na.rm = T)]
      metaGene[, dist.Mean := rowMeans2(dataHolder$fc, na.rm = T)]
      metaGene[, dist.SD := Rfast::rowVars(dataHolder$fc, na.rm = T, std = T)]
      
      # Precompute z-scores. Need to maintain logFC and adj.pv for user-selected filtering
      dataHolder$zscore <- (dataHolder$fc - metaGene$dist.Mean) / metaGene$dist.SD
      
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
      metaData$sf.Val <- metaData$sf.Val %>% as.factor
      metaData$sf.ValLongUri <- metaData$sf.ValLongUri %>% as.factor
      metaData$cf.Cat <- metaData$cf.Cat %>% as.factor
      metaData$cf.CatLongUri <- metaData$cf.CatLongUri %>% as.factor
      metaData$cf.ValLongUri <- metaData$cf.ValLongUri %>% as.factor
      metaData$cf.BaseLongUri <- metaData$cf.BaseLongUri %>% as.factor
      
      goTerms <- queryMany(metaGene$entrez.ID, scopes = 'entrezgene', fields = 'go', species = taxon)
      
      goTerms <- rbindlist(lapply(c('CC', 'BP', 'MF'), function(cat) {
        gocat <- paste0('go.', cat)
        rbindlist(lapply(1:nrow(goTerms), function(indx) {
          row <- goTerms@listData[[gocat]][[indx]]
          if(!is.null(row))
            data.frame(entrez.ID = goTerms@listData$query[indx], category = cat, id = row$id, term = row$term)
        }))
      }))
      
      new('EData', taxon = taxon, data = dataHolder,
          experiment.meta = metaData, gene.meta = metaGene, go = unique(goTerms))
    })
    
    names(DATA.HOLDER) <- getConfig(key = 'taxa')$core
    
    saveRDS(DATA.HOLDER, '/space/scratch/jsicherman/Thesis Work/data/DATA.HOLDER.rds')
  }
  
  ONTOLOGIES.DEFS <- fixOntoGenes()
}

if(!exists('CACHE.BACKGROUND')) {
  # Pre-load all ontology expansions
  if(file.exists('/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND.rds'))
    CACHE.BACKGROUND <- readRDS('/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND.rds')
  else {
    CACHE.BACKGROUND <- lapply(names(DATA.HOLDER), precomputeTags)
    names(CACHE.BACKGROUND) <- names(DATA.HOLDER)
    
    saveRDS(CACHE.BACKGROUND, '/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND.rds')
  }
}

if(!exists('TAGS')) {
  TAGS <- lapply(names(CACHE.BACKGROUND), getTags) %>% `names<-`(names(CACHE.BACKGROUND))
}

if(!exists('NULLS')) {
  mFiles <- list.files('/space/scratch/jsicherman/Thesis Work/data/nulls')
  NULLS <- lapply(c('any', names(DATA.HOLDER)), function(taxa) {
    mDT <- NULL
    for(f in mFiles[startsWith(mFiles, taxa)]) {
      message(paste('Loading', f))
      N <- gsub(paste0(taxa, '_(\\d+)\\.rds'), '\\1', f)
      dt <- readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/nulls/', f)) %>%
        .[, .(st.mean, st.sd), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
        setnames(4:5, paste0(c('M', 'S'), N))
      if(is.null(mDT))
        mDT <- dt
      else
        mDT <- merge(mDT, dt, by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all = T, sort = F)
    }
    
    mDT
  }) %>% `names<-`(c('any', names(DATA.HOLDER)))
  rm(mFiles)
}

if(!exists('DRUGBANK')) {
  DRUGBANK <- loadDrugbank()
}
