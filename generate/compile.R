devtools::load_all()
library(gemma.R)

# some fields require a logged in user to access
gemma_user <- readLines('generate/auth')
gemma.R::set_gemma_user(gemma_user[1],gemma_user[2])


dir.create(RAWDIR,showWarnings = FALSE,recursive = TRUE)

# set up cache and memoise so repeated calls to platforms are faster
# alternative is to process everything before that but lets stick to
# the easy way for now
options('gemma.cache' = file.path(DATADIR,'gemma_cache'))
options('gemma.memoised' = TRUE)


dir.create(file.path(RAWDIR,'metadata'),showWarnings = FALSE)
lapply(c('human', 'mouse', 'rat'), function(taxon){
  print(taxon)
  # meta.platformCoverage is a compilation species genes and platforms. this is
  # only used for 
  poke_call = gemma.R::get_taxon_datasets(taxon,limit = 1)
  
  # temporary for debuggin
  # attributes(poke_call)$totalElements = 200
  
  seq(0,attributes(poke_call)$totalElements,100) %>% lapply(function(x){
    gemma.R::gemma_call('taxa/{taxa}/datasets/?offset={offset}&limit=100&sort=%2Bid',taxa = 'human',offset = x)$data
  }) %>% do.call(c,.) -> all_datasets_raw
  
  seq(0,attributes(poke_call)$totalElements,100) %>% lapply(function(x){
    out = gemma.R::get_taxon_datasets(taxon,offset = x,limit = 100)
  }) %>% do.call(rbind,.) -> all_datasets
  
  
  seq_len(nrow(all_datasets)) %>% lapply(function(i){
    print(i)
    dataset <- all_datasets[i,]
    dataset_raw <- all_datasets_raw[[i]]
    
    dataset_annotations <- gemma.R::get_dataset_annotations(dataset$ee.ID)
    
    differential <- gemma.R::get_dataset_differential_expression_analyses(dataset$ee.ID)
    
    if(length(differential)==0){
      return(NULL)
    }
    
    samples <- tryCatch(gemma.R::get_dataset_samples(dataset$ee.ID),
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
      sf.NumSample <- dataset$ee.SampleCount
    }
    
    platform <- gemma.R::get_dataset_platforms(dataset$ee.ID)
    platform_annots <- platform$platform.ID %>% lapply(function(x){
      gemma.R::get_platform_annotations(x)
    }) %>% do.call(rbind,.)
    # data frame is manually re-created from the components to ensure
    # compatibility with downstream use by replicating nathaniel's structure
    data.frame(ee.ID = dataset$ee.ID,
               ee.qScore = dataset$geeq.qScore,
               ee.sScore = dataset$geeq.sScore,
               rsc.ID = differential$rsc.ID,
               ee.Name = dataset$ee.ShortName,
               ee.Source = dataset$ee.Database,
               ee.Scale = NA, # not accessible with the current API, not used in the dataset either
               ee.Reprocessed = NA, # not accessible with the current API, it's used in processing but we'll see...
               ef.IsBatchConfounded = dataset$geeq.batchConfound == -1,
               ad.ID = platform$platform.ID,
               ad.Type = platform$technology.Type,
               ad.NumGenes = platform_annots$GeneSymbols %>% unique %>% length(), # this counts genes a bit differently than what gemma displays as I think it should be closer to the intended purpose. gemma counts by splitting genes that are aligned to the same probeset which presumably shouldn't be used as differential expression results. need to reconsider
               ee.NumSample = dataset$ee.SampleCount,
               sf.NumSample =  sf.NumSample,
               cf.Cat = differential$cf.Cat,
               cf.CatLongUri = differential$cf.CatLongUri,
               cf.Baseline = differential$cf.Baseline,
               cf.BaselineLongUri = differential$cf.BaseLongUri,
               cf.Val = differential$cf.Val,
               cf.ValLongUri = differential$cf.ValLongUri,
               sf.Val = differential$sf.Val,
               sf.ValLongUri = differential$sf.ValLongUri,
               ee.Tag = dataset_annotations$term.Name %>% paste(collapse = '; '),
               ee.TagLongUri = dataset_annotations$term.URL%>% paste(collapse = '; '),
               n.DE = differential$stats.DE)
    
  })  %>% do.call(rbind,.) -> new_metaData
  saveRDS(new_metaData, file.path(RAWDIR,'metadata',taxon,'.rds'))
  return(new_metaData)
}) -> species_metadata

saveRDS(species_metadata,file = file.path(RAWDIR,'species_metadata.rds'))
