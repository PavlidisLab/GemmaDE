


devtools::load_all()
source(here::here("main/dependencies.R"))

genes = c('RPS4Y1','EIF1AY','DDX3Y','KDM5D','XIST')
taxa = 'human'

#* DE Search
#* 
#* @param genes
#* @param taxa
#* @param max_dist
#* @get /de_search
de_search = function(req,
                     genes = NULL,
                     taxa =NULL,
                     max_dist = 1.5,
                     confounds = FALSE,
                     multifunctionality = TRUE,
                     geeq = FALSE,
                     p_threshold = 0.05,
                     categories = c("age", "behavior", "biological process", "biological sex", 
                                    "cell type", "clinical history", "diet", "disease", "environmental history", 
                                    "environmental stress", "genotype", "medical procedure", "molecular entity", 
                                    "organism part", "phenotype", "sex", "temperature", "treatment"),
                     cores =8){
  browser()
  if(is.null(genes)){
    genes = req$HEADERS['genes'] %>% strsplit(',') %>% {.[[1]]}
  }
  if(is.null(taxa)){
    taxa = req$HEADERS['taxa'] %>% strsplit(',') %>% {.[[1]]}
  }
  if(is.null(caregories)){
    caregories = req$HEADERS['caregories'] %>% strsplit(',') %>% {.[[1]]}
  }
  
  tictoc::tic()
  genes <- processGenes(genes,taxa)
  
  experiments <- taxa %>% 
    parallel::mclapply(function(t){
      vsmSearch(genes[taxon == t, entrez.ID],
             taxa = t,
             confounds = confounds,
             filter = NULL,
             mfx = multifunctionality,
             geeq = geeq,
             p_threshold = p_threshold)
    },mc.cores = cores)
  names(experiments) = taxa
  
  conditions <- taxa %>% lapply(function(t){
    enrich(experiments[[t]], taxa = t, dist = max_dist,categories = categories,cores = cores) %>% 
      data.table::setnames(genes[taxon == t, gene.realName],
               genes[taxon == t, identifier],
               skip_absent = T)
  }) %>% data.table::rbindlist(fill = TRUE) %>% 
    reorderTags3() %>%
    .[, lapply(.SD, mean, na.rm = T), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
  
  geneInfo <- genes %>%
    data.table::copy() %>%
    data.table::setnames("identifier", "gene.Name")
  mGenes <- genes %>% data.table::copy()
  
  # components of the endSuccess function
  exps <- lapply(experiments, "[[", "rn") %>% unlist()
  
  tmp <- data.table::rbindlist(lapply(taxa, function(i) {
    DATA.HOLDER[[i]]@experiment.meta[rsc.ID %in% exps, .(rsc.ID, ee.ID, ee.Name, ee.NumSample, ef.IsBatchConfounded)]
  }))
  
  
  
  conditions[, `Condition Comparison` := paste0( cf.BaseLongUri, " vs. ", cf.ValLongUri)]
  
  tmp <- tmp %>%
    merge(getTags(taxa, exps), sort = F) %>%
    .[, N := length(unique(ee.ID)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    .[N < 0.03 * nrow(tmp)] # Get rid of contrasts that overlap in more than 3% experiments
  
  
  tmp[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, N)] %>%
    unique() %>%
    merge(conditions, by = c("cf.Cat", "cf.BaseLongUri", "cf.ValLongUri"), sort = F) %>%
    data.table::setnames(c("stat", "score", "distance"), c("Effect Size", "Test Statistic", "Ontology Steps")) ->
    conditions
  
  # out of endSuccess

  getPercentageStat <- function(x, n = 1){
    x / n
  }
  conditions[,'Test Statistic'] <- apply(conditions[,'Test Statistic'], 2, getPercentageStat, n = nrow(geneInfo))
  
  tictoc::toc()
  return(conditions)
}
