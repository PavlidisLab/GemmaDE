devtools::load_all()
library(gemma.R)
library(magrittr)

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


dir.create(file.path(RAWDIR,'metadata'),showWarnings = FALSE)

# objects to create for each species
# taxon,dataHolder,metaData,metaGene,goTerms
lapply(c('human', 'mouse', 'rat'), function(taxon){
  print(taxon)
  # poke call to get the number of elements
  poke_call = gemma.R::get_taxon_datasets(taxon,limit = 1)
  
  # temporary for debuggin
  attributes(poke_call)$totalElements = 99
  
  seq(0,attributes(poke_call)$totalElements,100) %>% lapply(function(x){
    gemma.R::gemma_call('taxa/{taxa}/datasets/?offset={offset}&limit=100&sort=%2Bid',taxa = 'human',offset = x)$data
  }) %>% do.call(c,.) -> all_datasets_raw
  
  seq(0,attributes(poke_call)$totalElements,100) %>% lapply(function(x){
    out = gemma.R::get_taxon_datasets(taxon,offset = x,limit = 100)
  }) %>% do.call(rbind,.) -> all_datasets
  
  
  all_datasets$experiment.ID %>% lapply(function(x){
    gemma.R::get_dataset_platforms(x)$platform.ID
  }) %>% unlist %>% unique -> all_platform_ids
  
  all_platform_annotations <- all_platform_ids %>% lapply(function(x){
    gemma.R::get_platform_annotations(x)
  })
  names(all_platform_annotations) <- all_platform_ids
  

  null_platforms <- all_platform_annotations %>% purrr::map_lgl(is.null)
  all_platform_annotations <- all_platform_annotations[!null_platforms]
  
  
  # all_genes <- all_platform_annotations %>%
  #   purrr::map('GemmaIDs') %>%
  #   unlist %>% unique %>%
  #   strsplit('|',fixed= TRUE) %>% unlist %>% unique %>% as.integer() %>% sort
  
  
  
  
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
     all_platform_annotations[[x]]
    }) %>% do.call(rbind,.)
    # data frame is manually re-created from the components to ensure
    # compatibility with downstream use by replicating nathaniel's structure
    data.frame(ee.ID = dataset$experiment.ID,
               ee.qScore = dataset$geeq.qScore,
               ee.sScore = dataset$geeq.sScore,
               rsc.ID = differential$rsc.ID,
               result.id = differential$result.ID,
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
               ee.TagLongUri = dataset_annotations$term.URL%>% paste(collapse = '; '))
    
  })  %>% do.call(rbind,.) -> new_metaData
  
  new_metaDataTable <- new_metaData %>% data.table
  
  new_metaDataTable[is.na(cf.BaseLongUri), cf.BaseLongUri := cf.Baseline]
  new_metaDataTable[is.na(cf.ValLongUri), cf.ValLongUri := cf.Val]
  
  # remove NAs
  new_metaDataTable = new_metaDataTable[!is.na(cf.BaseLongUri),]
  
  
  
  
  
  saveRDS(new_metaDataTable, file.path(RAWDIR,'metadata',paste0(taxon,'.rds')))
  return(new_metaDataTable)
}) -> species_metadata

# needs calculation for metadata
# n.DE number of differentially expressed genes per difExp




saveRDS(species_metadata,file = file.path(RAWDIR,'species_metadata.rds'))
