# Rscript generate/compile.R >gen.log 2>gen_err.log
# cat gen.log;echo ""

print('data compilation')
devtools::load_all()
library(gemma.R)
library(magrittr)
library(dplyr)

load_from_backups = FALSE
trust_cache = TRUE


mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% dplyr::mutate(...)
  .data
}
# some fields require a logged in user to access
# gemma_user <- readLines('generate/auth')
# gemma.R::set_gemma_user(gemma_user[1],gemma_user[2])
# rm(gemma_user)

dir.create(RAWDIR,showWarnings = FALSE,recursive = TRUE)

# set up cache and memoise so repeated calls to platforms are faster
# alternative is to process everything before that but lets stick to
# the easy way for now
options('gemma.cache' = CACHEDIR)

if(!trust_cache){
  # just in case...
  warning("You are trying to delete the cache. It's currently disabled for your own good")
  # gemma.R::forget_gemma_memoised()
}

options('gemma.memoised' = TRUE)

dir.create(file.path(RAWDIR,'midway_backups'),showWarnings = FALSE,recursive = TRUE)

# objects to create for each species
# taxon,dataHolder,metaData,metaGene,goTerms
lapply(c('human', 'mouse', 'rat'), function(taxon){
  dir.create(file.path(RAWDIR,'midway_backups',taxon),showWarnings =FALSE,recursive = TRUE)
  print(taxon)
  # poke call to get the number of elements
  poke_call = gemma.R::get_taxon_datasets(taxon,limit = 1)
  
  # temporary for debuggin
  # attributes(poke_call)$totalElements = 99
  
  print('/nGetting all datasets')
  bc_path = file.path(RAWDIR,'midway_backups',taxon,'all_datasets.rds')
  if(load_from_backups && file.exists(bc_path)){
    all_datasets = readRDS(bc_path)
  } else {
    seq(0,attributes(poke_call)$totalElements,100) %>% lapply(function(x){
      cat('=')
      out = gemma.R::get_taxon_datasets(taxon,offset = x,limit = 100)
    }) %>% do.call(rbind,.) -> all_datasets
    saveRDS(all_datasets,bc_path)
  }
  
  print('\nGetting platforms')
  bc_path = file.path(RAWDIR,'midway_backups',taxon,'platform_ids.rds')
  if(load_from_backups && file.exists(bc_path)){
    all_platform_ids <- readRDS(bc_path)
  } else {
    all_datasets$experiment.ID %>% lapply(function(x){
      cat('=')
      gemma.R::get_dataset_platforms(x)$platform.ID
    }) %>% unlist %>% unique -> all_platform_ids
    saveRDS(all_platform_ids, bc_path)
  }
  
  
  
  print('/nGetting platform annotations')
  bc_path = file.path(RAWDIR,'midway_backups',taxon,'all_platform_annotations.rds')
  if(load_from_backups && file.exists(bc_path)){
    all_platform_annotations <- readRDS(bc_path)
  } else {
    all_platform_annotations <- all_platform_ids %>% lapply(function(x){
      cat('=')
      gemma.R::get_platform_annotations(x)
    })
    names(all_platform_annotations) <- all_platform_ids
    saveRDS(all_platform_annotations,bc_path)
  }
  
  
  
  null_platforms <- all_platform_annotations %>% purrr::map_lgl(is.null)
  all_platform_annotations <- all_platform_annotations[!null_platforms]
  
  print('/nCompiling the datasets table')
  
  bc_path = file.path(RAWDIR,'midway_backups',taxon,'metadata_first_pass.rds')
  if(load_from_backups && file.exists(bc_path)){
    contrast_metaData <- readRDS(bc_path)
  } else {
    seq_len(nrow(all_datasets)) %>% lapply(function(i){
      cat("=")
      dataset <- all_datasets[i,]
      
      dataset_annotations <- gemma.R::get_dataset_annotations(dataset$experiment.ID)
      differential <- gemma.R::get_dataset_differential_expression_analyses(dataset$experiment.ID)
      
      if(length(differential)==0){
        return(NULL)
      }
      
      samples <- tryCatch(gemma.R::get_dataset_samples(dataset$experiment.ID),
                          error = function(e){
                            if(e$message == "500: Internal server error."){
                              return(list())
                            }else{
                              stop(e$message)
                            }
                          })
      
      # count sf.NumSample
      if (any(c('DE_Include','DE_Exclude') %in% samples$sample.FactorValues[[1]]$name)){
        samples$sample.FactorValues %>% purrr::map('name') %>% purrr::map_lgl(function(x){
          'DE_Include' %in% x
        }) %>% sum -> sf.NumSample
      } else{
        sf.NumSample <- dataset$experiment.SampleCount
      }
      
      platform <- gemma.R::get_dataset_platforms(dataset$experiment.ID)
      platform_annots <- platform$platform.ID %>% lapply(function(x){
        all_platform_annotations[[as.character(x)]]
      }) %>% do.call(rbind,.)
      # data frame is manually re-created from the components to ensure
      # compatibility with downstream use by replicating nathaniel's structure
      data.frame(ee.ID = dataset$experiment.ID,
                 ee.qScore = dataset$geeq.qScore,
                 ee.sScore = dataset$geeq.sScore,
                 rsc.ID = paste0('RSCID.',differential$result.ID,'.',differential$contrast.id),
                 result.id = differential$result.ID,
                 contrast.id = differential$contrast.id,
                 ee.Name = dataset$experiment.ShortName,
                 ee.Source = dataset$experiment.Database,
                 ee.Scale = NA, # not accessible with the current API, not used in the dataset either
                 ee.Reprocessed = dataset$geeq.rawData==1, # not sure if this field matches 1-1
                 ef.IsBatchConfounded = dataset$geeq.batchConfound == -1,
                 ad.ID = paste0(platform$platform.ID,collapse = '; '),
                 ad.Type = paste0(platform$technology.Type,collapse = '; '),
                 ad.NumGenes = platform_annots$NCBIids %>% {.[.!='' | grepl('|',.,fixed = TRUE)]} %>% unique %>% length, # this counts genes a bit differently than what gemma displays as I think it should be closer to the intended purpose. gemma counts by splitting genes that are aligned to the same probeset which presumably shouldn't be used as differential expression results. need to reconsider
                 ee.NumSample = dataset$experiment.SampleCount,
                 sf.NumSample =  sf.NumSample,
                 cf.Cat = differential$baseline.category,
                 cf.CatLongUri = differential$baseline.categoryURI,
                 cf.Baseline = differential$baseline.factorValue,
                 cf.BaseLongUri = differential$baseline.factorValueURI,
                 cf.Val = differential$experimental.factorValue,
                 cf.ValLongUri = differential$experimental.factorValueURI,
                 sf.Val = differential$subsetFactor.factorValue,
                 sf.ValLongUri = differential$subsetFactor.factorValueURI,
                 ee.Tag = dataset_annotations$term.Name %>% paste(collapse = '; '),
                 ee.TagLongUri = dataset_annotations$term.URI%>% paste(collapse = '; '))
      
    })  %>% do.call(rbind,.) -> contrast_metaData
    
    saveRDS(contrast_metaData, bc_path)
  }
  contrast_metaData %<>% mutate_cond(is.na(cf.BaseLongUri), cf.BaseLongUri = cf.Baseline)
  contrast_metaData %<>% mutate_cond(is.na(cf.ValLongUri), cf.ValLongUri = cf.Val)
  contrast_metaData <- contrast_metaData[!is.na(contrast_metaData$cf.BaseLongUri),]
  
  
  differentials <- contrast_metaData$result.id %>% unique
  
  print("/nGetting differential expression data")
  bc_path =  file.path(RAWDIR,'midway_backups',taxon,'differentials.rds')
  bc_path2 = file.path(RAWDIR,'midway_backups',taxon,'differentials_filtered.rds')
  if(load_from_backups && file.exists(bc_path) && !file.exists(bc_path2)){
    all_differential_values <- readRDS(bc_path)
  } else if(!(file.exists(bc_path2) && load_from_backups)){
    all_differential_values <- differentials %>% lapply(function(x){
      cat('=')
      tryCatch(
        gemma.R::get_differential_expression_values(resultSet = x)[[1]],
        error = function(e){
          if(grepl('502',e$message)){
            return(data.frame())
          } else{
            stop(e$message)
          }
        }
      )
      
    })
    names(all_differential_values)= differentials
    saveRDS(all_differential_values, bc_path)
  }
  
  print('/nProcessing differential expression data')
  if(load_from_backups && file.exists(bc_path2)){
    all_differential_values <- readRDS(bc_path2)
  } else{
    bad_differentials = all_differential_values %>% sapply(nrow) %>% {.==0}
    all_differential_values = all_differential_values[!bad_differentials]
    
    # remove duplicates, calculate additional statistics
    # looking at nathaniel's old freeze, it appears that 
    # probeset with the lowest p value was selected to
    # to represent a gene. need a deeper look to make sure
    # but assuming that is the case for now -Ogan
    all_differential_values %<>% lapply(function(x){
      cat('=')
      x %<>% dplyr::filter(NCBIid!='' & !grepl('|',NCBIid,fixed = TRUE))
      p_value_cols = names(x)[grepl('contrast.*?pvalue',names(x))]
      logfc_cols = names(x)[grepl('contrast.*?log2fc',names(x))]
      stat_cols =names(x)[grepl('contrast.*?tstat',names(x))]
      dup_genes = x$NCBIid[duplicated(x$NCBIid)]
      out = x %>% dplyr::filter(!NCBIid %in% dup_genes) %>% dplyr::select(-Probe,-pvalue,-corrected_pvalue,-rank)
      # add dups with lowest p value
      if (length(p_value_cols)==0){
        return(data.frame())
      }
      
      lapply(seq_along(p_value_cols),function(i){
        x %>% dplyr::filter(NCBIid %in% dup_genes) %>% dplyr::arrange(!!sym(p_value_cols[i])) %>% 
          dplyr::filter(!duplicated(NCBIid)) %>%
          dplyr::select(NCBIid,GeneSymbol,GeneName,!!sym(logfc_cols[i]),!!sym(stat_cols[i]),!!sym(p_value_cols[i]))
      }) %>% {
        if(length(.)>1){
          out = .[[1]]
          for (j in seq_len(length(.)-1)){
            out = merge(out,.[[j+1]])
          }
          out
        } else if (length(.) == 1){
          .[[1]]
        }
      } -> to_append
      
      out = rbind(out,to_append)
      
      # calculate additional stats
      for(i in seq_along(p_value_cols)){
        out[[paste0(p_value_cols[[i]],'_adjusted')]] = p.adjust(out[[p_value_cols[i]]],'BH')
      }
      return(out)
    })
    bad_differentials = all_differential_values %>% sapply(nrow) %>% {.==0}
    all_differential_values = all_differential_values[!bad_differentials]
    
    saveRDS(all_differential_values, bc_path2)
  }
  
  # construct the dataHolder object-----------
  dh = list()
  
  ncbi_ids = all_differential_values %>% lapply(function(x){
    as.integer(x$NCBIid)
  }) %>% unlist %>% unique %>% sort
  
  print('calculating adj.p values')
  bc_path = file.path(RAWDIR,'midway_backups',taxon,'adj.pv.rds')
  if(load_from_backups && file.exists(bc_path)){
    adj.pv = readRDS(bc_path)
  } else{
    lapply(seq_len(nrow(contrast_metaData)),function(i){
      result_id = contrast_metaData$result.id[i]
      contrast_id = contrast_metaData$contrast.id[i]
      differential = all_differential_values[[as.character(result_id)]]
      differential[[
        glue::glue('contrast_{contrast_id}_pvalue_adjusted')
      ]][match(ncbi_ids,differential$NCBIid)]
    }) %>% {names(.) = contrast_metaData$rsc.ID;.}%>% do.call(cbind,.) -> adj.pv
    rownames(adj.pv) = ncbi_ids
    full_na_genes =  adj.pv %>% apply(1,function(x){all(is.na(x))})
    full_na_samples = adj.pv %>% apply(2,function(x){all(is.na(x))})
    adj.pv = adj.pv[!full_na_genes,!full_na_samples]
    
    saveRDS(adj.pv, file.path(bc_path))
  }
  
  # remove any filtered gene from the ncbi_ids
  ncbi_ids = rownames(adj.pv)

  dh$adj.pv = adj.pv
  
  
  # remove differentials without differential expression from the metadata
  contrast_metaData %<>% dplyr::filter(rsc.ID %in% colnames(dh$adj.pv))
  
  dh$fc = lapply(seq_len(nrow(contrast_metaData)),function(i){
    result_id = contrast_metaData$result.id[i]
    contrast_id = contrast_metaData$contrast.id[i]
    differential = all_differential_values[[as.character(result_id)]]
    differential[[
      glue::glue('contrast_{contrast_id}_log2fc')
    ]][match(ncbi_ids,differential$NCBIid)]
  }) %>% {names(.) = contrast_metaData$rsc.ID;.}%>% do.call(cbind,.) %>% 
    {rownames(.) = ncbi_ids;.}
  
  
  assertthat::assert_that(all(contrast_metaData$rsc.ID == colnames(adj.pv)))
  contrast_metaData$n.DE = matrixStats::colSums2(adj.pv <= 0.05, na.rm = T)
  contrast_metaData$mean.fc = matrixStats::colMeans2(abs(dh$fc), na.rm = T)
  contrast_metaData$mean.up = matrixStats::colMeans2(dh$fc * ifelse(dh$fc > 0, 1, NA), na.rm = T)
  contrast_metaData$mean.down = matrixStats::colMeans2(dh$fc * ifelse(dh$fc < 0, 1, NA), na.rm = T)
  
  saveRDS(contrast_metaData, file.path(RAWDIR,'midway_backups',taxon,'metadata_final.rds'))
  
  print('\ncompiling gene metadata')
  # the output order isn't guaranteed to be identical to the input order for get_genes
  seq(1,length(ncbi_ids),50) %>% lapply(function(i){
    cat('=')
    gemma.R::get_genes(ncbi_ids[seq(i,min(i+49,length(ncbi_ids)))],raw= TRUE, memoised = FALSE)
  }) %>% do.call(c,.) -> all_genes
  
  # not used, can add later if needed
  # ncbi_ids %>% lapply(function(id){
  #   print(id)
  #   gemma.R::get_gene_locations(id)
  # }) %>% do.call(rbind,.)
  
  gene_metaData = data.frame(
    entrez.ID = all_genes %>% purrr::map_int('ncbiId'),
    gene.ID = all_genes %>% purrr::map_int('id'),
    ensembl.ID =  all_genes %>% purrr::map_chr(function(x){if(is.null(x$ensemblId)){NA_character_}else{x$ensemblId}}),
    gene.Name = all_genes %>% purrr::map_chr('officialSymbol'),
    alias.Name = "", # not populated by gemma api
    gene.Desc = all_genes %>% purrr::map_chr('officialName'),
    gene.Type = NA, # not populated by gemma api
    gene.Chromosome = NA, # populated by gemmaAPI but takes quite long and we don't need it
    mfx.Rank = all_genes %>% purrr::map_dbl('multifunctionalityRank'),
    species = all_genes %>% purrr::map('taxon') %>% purrr::map_chr('commonName'))
  gene_metaData = gene_metaData[match(ncbi_ids,gene_metaData$entrez.ID),]
  gene_metaData$n.DE = matrixStats::rowSums2(dh$adj.pv <= 0.05, na.rm = T)
  gene_metaData$dist.Mean =  rowMeans(dh$fc[, contrast_metaData$ee.Reprocessed], na.rm = T)
  gene_metaData$dist.SD = Rfast::rowVars(dh$fc[, contrast_metaData$ee.Reprocessed], na.rm = T, std = T)
  
  
  # identify genes from the wrong species that snuck inside the metadata
  bad_genes = !gene_metaData$species %in% taxon
  assertthat::assert_that(all(gene_metaData$entrez.ID == rownames(dh$adj.pv)))
  assertthat::assert_that(all(gene_metaData$entrez.ID == rownames(dh$fc)))
  
  # identify contrasts that includes those genes
  bad_contrasts = dh$adj.pv[bad_genes,] %>% apply(2,function(x){!all(is.na(x))})
  
  # remove bad genes and bad contrasts from the results. this should be moved
  # up, before dealing with inflated data for better efficiency later but for this
  # corrective run, it's easier to keep it here
  gene_metaData = gene_metaData[!bad_genes,]
  dh$adj.pv = dh$adj.pv[!bad_genes,]
  dh$fc = dh$fc[!bad_genes,]
  dh$adj.pv = dh$adj.pv[,!bad_contrasts]
  dh$fc = dh$fc[,!bad_contrasts]
  contrast_metaData = contrast_metaData[!bad_contrasts,]
  assertthat::assert_that(all(contrast_metaData$rsc.ID == colnames(dh$adj.pv)))
  
  
  
  print('calculating z scores')
  dh$zscore <- (dh$fc - gene_metaData$dist.Mean) / gene_metaData$dist.SD
  
  print('getting go terms')
  
  go_terms <- mygene::queryMany(gene_metaData$entrez.ID, scopes = 'entrezgene', fields = 'go', species = taxon)
  
  go_terms <- data.table::rbindlist(lapply(c('CC', 'BP', 'MF'), function(cat) {
    gocat <- paste0('go.', cat)
    data.table::rbindlist(lapply(1:nrow(go_terms), function(indx) {
      row <- go_terms@listData[[gocat]][[indx]]
      if(!is.null(row))
        data.frame(entrez.ID = go_terms@listData$query[indx], category = cat, id = row$id, term = row$term)
    }), fill = T)
  }), fill = T)
  
  # metaData is contrast_metaData : complete
  # dataHolder is dh: finished
  # metaGene is gene_metaData : complete
  # goTerms is go_terms: finished
  
  print('final touches')
  contrast_metaData$ee.Name %<>% as.factor
  contrast_metaData$ee.Source %<>% as.factor
  contrast_metaData$ee.Scale %<>% as.factor
  contrast_metaData$ee.Tag %<>% as.factor
  contrast_metaData$ee.TagLongUri %<>% as.factor
  contrast_metaData$ad.Type %<>% as.factor
  contrast_metaData$ad.ID %<>% as.factor
  contrast_metaData$sf.Val %<>% as.factor
  contrast_metaData$sf.ValLongUri %<>% as.factor
  contrast_metaData$cf.Cat %<>% as.factor
  contrast_metaData$cf.CatLongUri %<>% as.factor
  contrast_metaData$cf.Baseline %<>% as.factor
  contrast_metaData$cf.Val %<>% as.factor
  contrast_metaData$cf.BaseLongUri %<>% as.factor
  contrast_metaData$cf.ValLongUri %<>% as.factor
  
  contrast_metaData %<>% data.table()
  gene_metaData %<>% data.table()
  
  
  out <- new('EData', taxon = taxon, data = dh,
             experiment.meta = contrast_metaData, gene.meta = gene_metaData, go = unique(go_terms))
  
  
  out$experiment.meta$nonNa.numGenes = out@data$fc %>% apply(2,\(x){
    sum(!is.na(x))
  })
  
  
  
  dir.create(file.path(RAWDIR,taxon),showWarnings = FALSE)
  
  
  
  
  saveRDS(out, file.path(RAWDIR,taxon,'Edata.rds'))
  return(out)
}) ->  data.holder
names(data.holder) = c('human','mouse','rat')


# this part is no longer necesarry and should be removed in the next run
# continuing to use this for now to avoid deleting the cache 
# here we load the multifunctionality scores. these are generated by Paul at the time
# of data update but I believe they'd require access to re-generate. below is the mysql
# comman used but not adding it into the workflow since it can't be replicated by
# arbitrary people
# mysql -ugemmaadmin -p -habe -B -e ‘select g.OFFICIAL_SYMBOL, g.NCBI_GENE_ID, g.TAXON_FK, m.SCORE, m.RANK, m.NUM_GO_TERMS from CHROMOSOME_FEATURE g inner join MULTIFUNCTIONALITY m on m.ID=g.MULTIFUNCTIONALITY_FK where g.TAXON_FK in (1,2,3)’ gemd > multifunctionality.dump.txt
# file.copy('/home/paul/gemma-work/multifunctionality.dump.txt',file.path(RAWDIR,'mfx.txt'))
# mfx_dump <- readr::read_tsv(file.path(RAWDIR,'mfx.txt'))
# 
# for (taxon in c('human','mouse','rat')){
#   mfx <- mfx_dump$RANK[match(data.holder[[taxon]]@gene.meta$entrez.ID,mfx_dump$NCBI_GENE_ID)]
#   data.holder[[taxon]]@gene.meta$mfx.Rank <- mfx
# }

saveRDS(data.holder, file.path(RAWDIR, 'DATA.HOLDER.rds'))


# conversion to file based
# matrixes are transposed before file based storage since
# data is stored column-wise and we access subsets of genes
# we also save the data in non transposed form as _contrast
# and load it to DATA.HOLDER to allow for accessing contrast
# level data easier if needed. One use case is calculation of
# null sets
fbm_path = file.path(DATADIR,'data_fbm')

for(taxon in names(data.holder)) {
  create_fbm(data.holder$human@data$fc,
             file.path(fbm_path,taxon,'fc'))
  
  create_fbm(data.holder$human@data$zscore,
             file.path(fbm_path,taxon,'zscore'))
  
  
  create_fbm(data.holder$human@data$adj.pv,
             file.path(fbm_path,taxon,'adj.pv'))
}


for (taxon in names(data.holder)) {
  data.holder[[taxon]]@data$zscore <- NULL
  data.holder[[taxon]]@data$adj.pv <- NULL
  data.holder[[taxon]]@data$fc <- NULL
}
saveRDS(data.holder, paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/'))


