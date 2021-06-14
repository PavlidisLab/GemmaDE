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
  
  ONTOLOGIES[, c('ChildNode', 'ParentNode')] <- NULL
  ONTOLOGIES.DEFS$Node <- NULL
  
  ONTOLOGIES$ChildNode_Long <- ONTOLOGIES$ChildNode_Long %>% as.factor
  ONTOLOGIES$ParentNode_Long <- ONTOLOGIES$ParentNode_Long %>% as.factor
  ONTOLOGIES$RelationType <- ONTOLOGIES$RelationType %>% as.factor
  ONTOLOGIES$OntologyScope <- ONTOLOGIES$OntologyScope %>% as.factor
  
  ONTOLOGIES.DEFS$Node_Long <- ONTOLOGIES.DEFS$Node_Long %>% as.factor
  ONTOLOGIES.DEFS$OntologyScope <- ONTOLOGIES.DEFS$OntologyScope %>% as.factor
}

# Load the lite versions if they're already created.
if(!exists('DATA.HOLDER')) {
  if(file.exists('/space/scratch/jsicherman/Thesis Work/data/DATA.HOLDER2.rds'))
    DATA.HOLDER <- readRDS('/space/scratch/jsicherman/Thesis Work/data/DATA.HOLDER2.rds')
  else {
    DATA.HOLDER <- lapply(getConfig(key = 'taxa')$core, function(taxon) {
      message(paste('Loading', taxon, 'metadata'))
      
      load(paste0('/space/grp/nlim/MDE/RDataRepo/Packaged/Current/', taxon, '/metadata.RDAT'))
      
      meta.platformCoverage <- melt(meta.platformCoverage, id.vars = 'gene.ID') %>%
        .[, variable := gsub('AD_(.*)', '\\1', variable)]
      
      # TODO Yes I know this is no longer optimal. In the old data dump, values were semicolon delimited...
      # It's too deeply integrated for me to change everything else at this point
      metaData <- meta.fvAnnot %>%
        merge(meta.full, by = c('rsc.ID', 'ee.ID'), sort = F) %>%
        merge(meta.platformCoverage[, .(n.detect = sum(value)), variable],
              by.x = 'ad.ID', by.y = 'variable', sort = F) %>%
        .[, .(ee.ID, rsc.ID, ee.Name, ee.Source, ee.Scale = ee.DataScale, ee.Reprocessed = ee.IsReprocessed,
              ef.IsBatchConfounded, ad.ID, ad.Type, ad.NumGenes = n.detect, ee.NumSample, sf.NumSample,
              ef.Cat, ef.CatLongUri = ef.CatUri, at.Cat, at.CatLongUri,
              at.Val, at.ValLongUri, at.Type = at.Type_1, at.Subtype = at.Type_2)] %>%
        .[, .(ee.ID, ee.Name, ee.Source, ee.Scale, ee.Reprocessed, ef.IsBatchConfounded,
              ad.ID, ad.Type, ad.NumGenes, ee.NumSample, sf.NumSample, cf.Cat = ef.Cat, cf.CatLongUri = ef.CatLongUri,
              cf.Baseline = paste0(.SD[at.Subtype == 'Baseline', at.Val], collapse = '; '),
              cf.BaseLongUri = paste0(.SD[at.Subtype == 'Baseline', at.ValLongUri], collapse = '; '),
              cf.Val = paste0(.SD[at.Subtype == 'Contrast', at.Val], collapse = '; '),
              cf.ValLongUri = paste0(.SD[at.Subtype == 'Contrast', at.ValLongUri], collapse = '; '),
              sf.Val = paste0(.SD[at.Type == 'SubsetFactor', at.Val], collapse = '; '),
              sf.ValLongUri = paste0(.SD[at.Type == 'SubsetFactor', at.ValLongUri], collapse = '; '),
              ee.Tag = paste0(.SD[at.Type %in% c('ExperimentTag', 'BioMaterial'), at.Val], collapse = '; '),
              ee.TagLongUri = paste0(.SD[at.Type %in% c('ExperimentTag', 'BioMaterial'), at.ValLongUri], collapse = '; ')), rsc.ID] %>%
        unique
        #dcast(... ~ at.Type + at.Subtype, list, value.var = c('at.Cat', 'at.CatLongUri', 'at.Val', 'at.ValLongUri')) %>%
        #setnames(gsub('at\\.(Cat|CatLongUri|Val|ValLongUri)_(BioMaterial|ExperimentTag|Factor|SubsetFactor)_(Baseline|Contrast)',
        #              'at.\\3.\\1', names(.))) %>%
        #setnames(gsub('at\\.(Cat|CatLongUri|Val|ValLongUri)_(BioMaterial|ExperimentTag|Factor|SubsetFactor)_Common',
        #              'at.\\2.\\1', names(.)))
      
      metaData[is.na(cf.BaseLongUri) | cf.BaseLongUri == 'NA', cf.BaseLongUri := cf.Baseline]
      metaData[is.na(cf.ValLongUri) | cf.ValLongUri == 'NA', cf.ValLongUri := cf.Val]
      
      clean <- function(baselines, baseUris) {
        unlist(lapply(1:length(baselines), function(i) {
          tmp <- parseListEntry(baseUris[i])
          tmp[is.na(tmp)] <- parseListEntry(baselines[i])[is.na(tmp)]
          paste0(tmp, collapse = '; ')
        }))
      }
      
      metaData[grepl('(^|; )(NA($|; ))+', cf.BaseLongUri), cf.BaseLongUri := clean(cf.Baseline, cf.BaseLongUri)]
      metaData[grepl('(^|; )(NA($|; ))+', cf.ValLongUri), cf.ValLongUri := clean(cf.Val, cf.ValLongUri)]
      
      # If it was a free-text entry, it becomes NA here.
      metaData[grepl('(^|; )(NA($|; ))+', cf.BaseLongUri), cf.BaseLongUri := gsub('; $', '', gsub('(^|; )(NA($|; ))+', '', cf.BaseLongUri))]
      metaData[grepl('(^|; )(NA($|; ))+', cf.ValLongUri), cf.ValLongUri := gsub('; $', '', gsub('(^|; )(NA($|; ))+', '', cf.ValLongUri))]
      
      # After filtering NAs, we should make them real NAs
      metaData[cf.BaseLongUri == '' | cf.BaseLongUri == 'NA', cf.BaseLongUri := NA]
      metaData[cf.ValLongUri == '' | cf.BaseLongUri == 'NA', cf.ValLongUri := NA]
      
      message('Loading data')
      
      dataHolder <- list(
        pv = readRDS(paste0('/space/grp/nlim/MDE/RDataRepo/Packaged/Current/', taxon, '/pv.RDS')) %>%
          `rownames<-`(gsub('GENE_(.*)', '\\1', rownames(.))),
        fc = readRDS(paste0('/space/grp/nlim/MDE/RDataRepo/Packaged/Current/', taxon, '/fc.RDS')) %>%
          `rownames<-`(gsub('GENE_(.*)', '\\1', rownames(.)))
      )
      
      message('Preprocessing')
      
      metaData <- metaData[order(match(rsc.ID, colnames(dataHolder$pv)))]
      
      dataHolder$fc[is.nan(dataHolder$fc)] <- NA
      dataHolder$pv[is.nan(dataHolder$pv)] <- NA
      
      # TODO validity of FDR correction based on number of genes each platform is capable of detecting?
      dataHolder$adj.pv <- sapply(1:ncol(dataHolder$pv), function(col)
        p.adjust(dataHolder$pv[, col], method = 'BH', n = metaData[col, ad.NumGenes])
      ) %>% `colnames<-`(colnames(dataHolder$pv))
      
      dataHolder$pv <- NULL
      
      metaData[, n.DE := colSums2(dataHolder$adj.pv <= 0.05, na.rm = T)]
      metaData[, mean.fc := colMeans2(abs(dataHolder$fc), na.rm = T)]
      metaData[, mean.up := colMeans2(dataHolder$fc * ifelse(dataHolder$fc > 0, 1, NA), na.rm = T)]
      metaData[, mean.down := colMeans2(dataHolder$fc * ifelse(dataHolder$fc < 0, 1, NA), na.rm = T)]
      
      metaData$ee.ID <- metaData$ee.ID %>% as.integer
      
      metaGene <- meta.gene[, .(entrez.ID, gene.ID = gemmaGene.ID,
                                ensembl.ID, gene.Name, alias.Name,
                                gene.Desc, gene.Type, gene.Chromosome,
                                mfx.Rank = ((mfx.Score - min(mfx.Score)) / (max(mfx.Score) - min(mfx.Score))) * (1 - 0) + 0)] # Rescale between 0 and 1
      
      metaGene[, n.DE := rowSums2(dataHolder$adj.pv <= 0.05, na.rm = T)]
      metaGene[, dist.Mean := rowMeans(dataHolder$fc[, metaData$ee.Reprocessed], na.rm = T)]
      metaGene[, dist.SD := Rfast::rowVars(dataHolder$fc[, metaData$ee.Reprocessed], na.rm = T, std = T)]
      
      # Precompute z-scores
      dataHolder$zscore <- (dataHolder$fc - metaGene$dist.Mean) / metaGene$dist.SD
      # TODO Maintained for GeneChaser # dataHolder$fc <- NULL
      
      message('Adding Gemma info')
      
      # Split into 500 ee.ID chunks (so the URI doesn't get too long) and fetch quality scores 
      # from Gemma for all experiments.
      metaData <- rbindlist(lapply(lapply(metaData$ee.ID %>% unique %>% split(ceiling(seq_along(.[]) / 500)),
                                          datasetInfo) %>% unlist(recursive = F), function(ee.ID) {
                                            data.table(ee.ID = ee.ID$id,
                                                       ee.qScore = ee.ID$geeq$publicQualityScore,
                                                       ee.sScore = ee.ID$geeq$publicSuitabilityScore)
                                          }), fill = T) %>% merge(metaData, by = 'ee.ID', all.y = T)
      
      goTerms <- queryMany(metaGene$entrez.ID, scopes = 'entrezgene', fields = 'go', species = taxon)
      
      goTerms <- rbindlist(lapply(c('CC', 'BP', 'MF'), function(cat) {
        gocat <- paste0('go.', cat)
        rbindlist(lapply(1:nrow(goTerms), function(indx) {
          row <- goTerms@listData[[gocat]][[indx]]
          if(!is.null(row))
            data.frame(entrez.ID = goTerms@listData$query[indx], category = cat, id = row$id, term = row$term)
        }), fill = T)
      }), fill = T)
      
      metaData$ee.Name <- metaData$ee.Name %>% as.factor
      metaData$ee.Source <- metaData$ee.Source %>% as.factor
      metaData$ee.Scale <- metaData$ee.Scale %>% as.factor
      metaData$ee.Tag <- metaData$ee.Tag %>% as.factor
      metaData$ee.TagLongUri <- metaData$ee.TagLongUri %>% as.factor
      metaData$ad.Type <- metaData$ad.Type %>% as.factor
      metaData$ad.ID <- metaData$ad.ID %>% as.factor
      metaData$sf.Val <- metaData$sf.Val %>% as.factor
      metaData$sf.ValLongUri <- metaData$sf.ValLongUri %>% as.factor
      metaData$cf.Cat <- metaData$cf.Cat %>% as.factor
      metaData$cf.CatLongUri <- metaData$cf.CatLongUri %>% as.factor
      metaData$cf.Baseline <- metaData$cf.Baseline %>% as.factor
      metaData$cf.Val <- metaData$cf.Val %>% as.factor
      metaData$cf.BaseLongUri <- metaData$cf.BaseLongUri %>% as.factor
      metaData$cf.ValLongUri <- metaData$cf.ValLongUri %>% as.factor
      
      new('EData', taxon = taxon, data = dataHolder,
          experiment.meta = metaData, gene.meta = metaGene, go = unique(goTerms))
    }) %>%
      setNames(getConfig(key = 'taxa')$core)
    
    saveRDS(DATA.HOLDER, '/space/scratch/jsicherman/Thesis Work/data/DATA.HOLDER2.rds')
  }
  
  ONTOLOGIES.DEFS <- fixOntoGenes()
}

if(!exists('CACHE.BACKGROUND')) {
  # Pre-load all ontology expansions
  if(file.exists('/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND2.rds'))
    CACHE.BACKGROUND <- readRDS('/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND2.rds')
  else {
    CACHE.BACKGROUND <- lapply(names(DATA.HOLDER), precomputeTags) %>%
      setNames(names(DATA.HOLDER))
    
    saveRDS(CACHE.BACKGROUND, '/space/scratch/jsicherman/Thesis Work/data/CACHE.BACKGROUND2.rds')
  }
}

if(!exists('NULLS')) {
  NULLS <- lapply(names(DATA.HOLDER), function(taxa) {
    readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/updated_nulls2/', taxa, '.rds')) %>%
      .[, .(rn, score.mean, score.sd)]
  }) %>% `names<-`(names(DATA.HOLDER))
}

if(!exists('DRUGBANK') && Sys.getenv('RSTUDIO') == '1') {
  DRUGBANK <- loadDrugbank()
}
