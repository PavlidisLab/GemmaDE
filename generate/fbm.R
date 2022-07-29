PROJDIR <- here::here()
DATADIR <- '/cosmos/data/project-data/GemmaDE'
FREEZEDIR <- '/cosmos/data/project-data/GemmaDE/gemma_freeze'
devtools::load_all()
dir.create(DATADIR,showWarnings = FALSE)

# source(file.path(PROJDIR, 'main', 'requirements.R'))

setClass('EData', representation(taxon = 'character', data = 'list',
                                 experiment.meta = 'data.table', gene.meta = 'data.table',
                                 go = 'data.table'))


.DATA_PATH <- file.path(DATADIR, 'DATA.HOLDER.rds')

# Load the lite versions if they're already created.
if(!exists('DATA.HOLDER')) {
  if(file.exists(.DATA_PATH))
    DATA.HOLDER <- readRDS(.DATA_PATH)
  else {
    DATA.HOLDER <- lapply(c('human', 'mouse', 'rat'), function(taxon) {
      message(paste('Loading', taxon, 'metadata'))
      
      # nathaniel's data files. need documentation on these - ogan
      # backed up everything in cosmos with the rest of the data. 
      # original path is kept in case it's needed later
      # objs = load(file.path('/space/grp/nlim/MDE/RDataRepo/Packaged/Current', taxon, 'metadata.RDAT'))
      objs = load(file.path(FREEZEDIR,'Current', taxon, 'metadata.RDAT'))
      
      
      meta.platformCoverage <- data.table::melt(meta.platformCoverage, id.vars = 'gene.ID') %>%
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
      
      # original paths to nathaniel's freeze is kept in case needed later
      # dataHolder <- list(
      #   pv = readRDS(paste0('/space/grp/nlim/MDE/RDataRepo/Packaged/Current/', taxon, '/pv.RDS')) %>%
      #     `rownames<-`(gsub('GENE_(.*)', '\\1', rownames(.))),
      #   fc = readRDS(paste0('/space/grp/nlim/MDE/RDataRepo/Packaged/Current/', taxon, '/fc.RDS')) %>%
      #     `rownames<-`(gsub('GENE_(.*)', '\\1', rownames(.)))
      # )
      
      dataHolder <- list(
        pv = readRDS(file.path(FREEZEDIR,'Current/', taxon, '/pv.RDS')) %>%
          `rownames<-`(gsub('GENE_(.*)', '\\1', rownames(.))),
        fc = readRDS(file.path(FREEZEDIR, 'Current/', taxon, '/fc.RDS')) %>%
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
      
      metaData[, n.DE := matrixStats::colSums2(dataHolder$adj.pv <= 0.05, na.rm = T)]
      metaData[, mean.fc := matrixStats::colMeans2(abs(dataHolder$fc), na.rm = T)]
      metaData[, mean.up := matrixStats::colMeans2(dataHolder$fc * ifelse(dataHolder$fc > 0, 1, NA), na.rm = T)]
      metaData[, mean.down := matrixStats::colMeans2(dataHolder$fc * ifelse(dataHolder$fc < 0, 1, NA), na.rm = T)]
      
      metaData$ee.ID <- metaData$ee.ID %>% as.integer
      
      metaGene <- meta.gene[, .(entrez.ID, gene.ID = gemmaGene.ID,
                                ensembl.ID, gene.Name, alias.Name,
                                gene.Desc, gene.Type, gene.Chromosome,
                                mfx.Rank = ((mfx.Score - min(mfx.Score)) / (max(mfx.Score) - min(mfx.Score))) * (1 - 0) + 0)] # Rescale between 0 and 1
      
      metaGene[, n.DE := matrixStats::rowSums2(dataHolder$adj.pv <= 0.05, na.rm = T)]
      metaGene[, dist.Mean := rowMeans(dataHolder$fc[, metaData$ee.Reprocessed], na.rm = T)]
      metaGene[, dist.SD := Rfast::rowVars(dataHolder$fc[, metaData$ee.Reprocessed], na.rm = T, std = T)]
      
      # Precompute z-scores
      dataHolder$zscore <- (dataHolder$fc - metaGene$dist.Mean) / metaGene$dist.SD
      # TODO Maintained for GeneChaser # dataHolder$fc <- NULL
      
      message('Adding Gemma info')
      
      # Split into 500 ee.ID chunks (so the URI doesn't get too long) and fetch quality scores 
      # from Gemma for all experiments.
      # TODO this is broken with the new Gemma API. Needs to be updated to get this info back
      metaData <- data.table::rbindlist(lapply(lapply(metaData$ee.ID %>% unique %>% split(ceiling(seq_along(.[]) / 500)),
                                          gemma.R::getDatasetsInfo) %>% unlist(recursive = F), function(ee.ID) {
                                            data.table::data.table(ee.ID = ee.ID$ee.ID,
                                                       ee.DescriptiveName = ee.ID$name,
                                                       ee.qScore = ee.ID$geeq.qScore,
                                                       ee.sScore = ee.ID$geeq.sScore)
                                          }), fill = T) %>% merge(metaData, by = 'ee.ID', all.y = T)
      
      goTerms <- mygene::queryMany(metaGene$entrez.ID, scopes = 'entrezgene', fields = 'go', species = taxon)
      
      goTerms <- data.table::rbindlist(lapply(c('CC', 'BP', 'MF'), function(cat) {
        gocat <- paste0('go.', cat)
        data.table::rbindlist(lapply(1:nrow(goTerms), function(indx) {
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
      setNames(c('human', 'mouse', 'rat'))
    
    
    saveRDS(DATA.HOLDER, file.path(DATADIR, 'DATA.HOLDER.rds'))
  }
}



# File-backed, gene-major matrices area HUGE upgrade compared to in-memory experiment-major matrices because 
# 1) way faster to access gene-slices (which are always longer than experiment slices), and 
# 2) data in memory can be reduced to < 1 GB (!)
if(class(DATA.HOLDER[[1]]@data$adj.pv) == 'matrix') {
  for(i in names(DATA.HOLDER)) {
    message(paste0('Converting in-memory matrices for ', i, ' to file-backed...'))
    
    file.remove(list.files(file.path(DATADIR, 'fbm', i), full.names = T))
    
    dir.create(file.path(DATADIR,'fbm',i),recursive = TRUE, showWarnings =FALSE)
    
    dimnames(DATA.HOLDER[[i]]@data$zscore) %>% 
      rev() %>%
      saveRDS(file.path(DATADIR, 'fbm', i, 'z.dimnames.rds'))
    
    DATA.HOLDER[[i]]@data$zscore <- bigstatsr::as_FBM(DATA.HOLDER[[i]]@data$zscore %>% t,
                                           backingfile = file.path(DATADIR, 'fbm', i, 'zscores'),
                                           is_read_only = T)$save()
    
    dimnames(DATA.HOLDER[[i]]@data$adj.pv) %>% 
      rev() %>%
      
      saveRDS(file.path(DATADIR, 'fbm', i, 'p.dimnames.rds'))
    DATA.HOLDER[[i]]@data$adj.pv <- bigstatsr::as_FBM(DATA.HOLDER[[i]]@data$adj.pv %>% t,
                                           backingfile = file.path(DATADIR, 'fbm', i, 'adjpvs'),
                                           is_read_only = T)$save()
  }
}

dimnames.FBM <- function(object, ...) {
  attr(object, '.dimnames')
}

# Remove z-scores and t-scores from DATA.HOLDER to allow faster loading with FBMs in load.R
for (taxon in names(DATA.HOLDER)) {
  DATA.HOLDER[[taxon]]@data$zscore <- NULL
  DATA.HOLDER[[taxon]]@data$adj.pv <- NULL
}


saveRDS(DATA.HOLDER, paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/'))
