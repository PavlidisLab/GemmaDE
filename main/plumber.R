


devtools::load_all()
source(here::here("main/dependencies.R"))
crs = 16

# genes=readr::read_csv('test.csv') %>% unlist
# genes = c('RPS4Y1','EIF1AY','DDX3Y','KDM5D','XIST')
# taxa = 'human'

#* DE Search
#* 
#* @param genes
#* @param taxa
#* @param max_dist
#* @post /de_search
#* @get /de_search
#* @serializer tsv
de_search = function(req = NULL, # this is the request object
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
                     cores = crs){
  # accept variables that can be long as body parts. long urls caused issues before
  # this should keep them at reasonable lengths

  args = rlang::fn_fmls(fn = de_search)[-1]
  # read arguments from request body if body is used instead of URLs
  for (x in names(args)){
    if(!is.null(req$body[[x]])){
      assign(x,req$body[[x]])
    }
  }
  if(!is.logical(geeq)){
    geeq = as.logical(toupper(geeq))
  }
  if(!is.logical(multifunctionality)){
    multifunctionality = as.logical(toupper(multifunctionality))
  }
  if(!is.logical(confounds)){
    confounds = as.logical(toupper(confounds))
  }
  tictoc::tic()
  genes <- processGenes(genes,taxa)
  print('vsmSearch')
  tictoc::tic()
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
  tictoc::toc()
  
  print('enrich')
  tictoc::tic()
  conditions <- taxa %>% lapply(function(t){
    enrich(experiments[[t]], taxa = t, dist = max_dist,categories = categories,cores = cores) %>% 
      data.table::setnames(genes[taxon == t, gene.realName],
               genes[taxon == t, identifier],
               skip_absent = T)
  }) %>% data.table::rbindlist(fill = TRUE) %>% 
    reorderTags3() %>%
    .[, lapply(.SD, mean, na.rm = T), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
  tictoc::toc()
  
  print('ending')
  tictoc::tic()
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
  
  tictoc::toc()
  return(conditions %>% 
           data.table::setnames(c("cf.Cat", "cf.BaseLongUri", "cf.ValLongUri"), c("Category", "Baseline", "Value")))
}

# de_search(genes=genes,taxa = taxa)
