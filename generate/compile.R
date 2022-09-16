devtools::load_all()
library(gemma.R)
library(magrittr)


mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
# some fields require a logged in user to access
gemma_user <- readLines('generate/auth')
gemma.R::set_gemma_user(gemma_user[1],gemma_user[2])
rm(gemma_user)

dir.create(RAWDIR,showWarnings = FALSE,recursive = TRUE)

# set up cache and memoise so repeated calls to platforms are faster
# alternative is to process everything before that but lets stick to
# the easy way for now
options('gemma.cache' = file.path(DATADIR,'gemma_cache'))
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
  attributes(poke_call)$totalElements = 99
  
  print('Getting all datasets')
  seq(0,attributes(poke_call)$totalElements,100) %>% lapply(function(x){
    out = gemma.R::get_taxon_datasets(taxon,offset = x,limit = 100)
  }) %>% do.call(rbind,.) -> all_datasets
  
  saveRDS(all_datasets,file.path(RAWDIR,'midway_backups',taxon,'all_datasets.rds'))
  
  print('Getting dataset platforms')
  all_datasets$experiment.ID %>% lapply(function(x){
    gemma.R::get_dataset_platforms(x)$platform.ID
  }) %>% unlist %>% unique -> all_platform_ids
  
  print('Getting platform annotations')
  all_platform_annotations <- all_platform_ids %>% lapply(function(x){
    gemma.R::get_platform_annotations(x)
  })
  names(all_platform_annotations) <- all_platform_ids
  saveRDS(all_platform_annotations,file.path(RAWDIR,'midway_backups',taxon,'all_platform_annotations.rds'))
  
  null_platforms <- all_platform_annotations %>% purrr::map_lgl(is.null)
  all_platform_annotations <- all_platform_annotations[!null_platforms]
  
  
  # seq_len(10) %>% lapply(function(i){
  seq_len(nrow(all_datasets)) %>% lapply(function(i){
    print(i)
    dataset <- all_datasets[i,]
    dataset_raw <- all_datasets_raw[[i]]
    
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
               ee.Reprocessed = NA, # not accessible with the current API, it's used in processing but we'll see...
               ef.IsBatchConfounded = dataset$geeq.batchConfound == -1,
               ad.ID = paste0(platform$platform.ID,collapse = '; '),
               ad.Type = paste0(platform$technology.Type,collapse = '; '),
               ad.NumGenes = platform_annots$GeneSymbols %>% unique %>% length(), # this counts genes a bit differently than what gemma displays as I think it should be closer to the intended purpose. gemma counts by splitting genes that are aligned to the same probeset which presumably shouldn't be used as differential expression results. need to reconsider
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
    
  })  %>% do.call(rbind,.) -> new_metaData
  
  saveRDS(new_metaData, file.path(RAWDIR,'midway_backups',taxon,'metadata_first_pass.rds'))
  
  new_metaData %<>% mutate_cond(is.na(cf.BaseLongUri), cf.BaseLongUri = cf.Baseline)
  new_metaData %<>% mutate_cond(is.na(cf.ValLongUri), cf.ValLongUri = cf.Val)
  new_metaData <- new_metaData[!is.na(new_metaData$cf.BaseLongUri),]
  
  
  differentials <- new_metaData$result.id %>% unique
  
  print("Getting differential expression data")
  all_differential_values <- differentials %>% lapply(function(x){
    print(x)
    gemma.R::get_differential_expression_values(resultSet = x)[[1]]
  })
  names(all_differential_values)= differentials
  saveRDS(all_differential_values, file.path(RAWDIR,'midway_backups',taxon,'differentials.rds'))
  bad_differentials = all_differential_values %>% sapply(nrow) %>% {.==0}
  all_differential_values = all_differential_values[!bad_differentials]
  # probeset selection
  # comparing with the old results, this doesn't appear to be the method used
  # print('Selecting probesets by using expression data')
  # for (x in unique(new_metaData$ee.ID)){
  #   expr <- gemma.R::get_dataset_expression(x)
  #   expr %<>%dplyr::filter(NCBIid!='' & !grepl('|',NCBIid,fixed = TRUE))
  #   dup_genes <- expr$NCBIid[duplicated(expr$NCBIid)]
  #   expr %<>% dplyr::filter(NCBIid %in% dup_genes)
  #   gene_var <- expr %>% ogbox::sepExpr() %>% {.[[2]]} %>% apply(1,function(x){var(na.omit(x))})
  #   probes_to_remove <- expr[order(gene_var,decreasing = TRUE),] %>% filter(duplicated(NCBIid)) %$%  Probe
  #   relevant_differentials = new_metaData %>% filter(ee.ID %in% x) %$% result.id
  #   all_differential_values[as.character(relevant_differentials)] %<>% lapply(function(diff){
  #     diff %>% dplyr::filter(!Probe %in% probes_to_remove)
  #   })
  # }
  

  # remove duplicates, calculate additional statistics
  # looking at nathaniel's old freeze, it appears that 
  # probeset with the lowest p value was selected to
  # to represent a gene. need a deeper look to make sure
  # but assuming that is the case for now -Ogan
  z = 1
  all_differential_values %<>% lapply(function(x){
    print(z)
    z<<- z + 1
    x %<>% dplyr::filter(NCBIid!='' & !grepl('|',NCBIid,fixed = TRUE))
    p_value_cols = names(x)[grepl('contrast.*?pvalue',names(x))]
    logfc_cols = names(x)[grepl('contrast.*?log2fc',names(x))]
    stat_cols =names(x)[grepl('contrast.*?tstat',names(x))]
    dup_genes = x$NCBIid[duplicated(x$NCBIid)]
    out = x %>% filter(!NCBIid %in% dup_genes) %>% select(-Probe,-pvalue,-corrected_pvalue,-rank)
    
    # add dups with lowest p value
    lapply(seq_along(p_value_cols),function(i){
      x %>% filter(NCBIid %in% dup_genes) %>% arrange(!!sym(p_value_cols[i])) %>% 
        filter(!duplicated(NCBIid)) %>%
        select(NCBIid,GeneSymbol,GeneName,!!sym(logfc_cols[i]),!!sym(stat_cols[i]),!!sym(p_value_cols[i]))
    }) %>% {
      if(length(.)>1){
        out = .[[1]]
        for (j in seq_len(length(.)-1)){
          out = merge(out,.[[j+1]])
        }
        out
     } else{
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
  

  
  lapply(seq_len(nrow(new_metaData)),function(i){
    result_id = new_metaData$result.id[i]
    contrast_id = new_metaData$contrast.id[i]
    differential = all_differential_values[[as.character(result_id)]]
    n.DE = differential[[glue::glue()]]
    
  })
  
  
  
  
  
  saveRDS(new_metaData, file.path(RAWDIR,'metadata',paste0(taxon,'.rds')))
  return(new_metaDataTable)
}) -> species_metadata

# needs calculation for metadata
# n.DE number of differentially expressed genes per difExp




saveRDS(species_metadata,file = file.path(RAWDIR,'species_metadata.rds'))
