

devtools::load_all()
source(here::here("main/dependencies.R"))


#* DE Search
#* 
#* @param genes
#* @param taxa
#* @param max_dist
#* @get /de_search
de_search = function(genes,
                     taxa,
                     max_dist = 1.5,
                     confounds = FALSE,
                     multifunctionality = TRUE,
                     geeq = FALSE,
                     p_threshold = 0.05,
                     categories = c("age", "behavior", "biological process", "biological sex", 
                                    "cell type", "clinical history", "diet", "disease", "environmental history", 
                                    "environmental stress", "genotype", "medical procedure", "molecular entity", 
                                    "organism part", "phenotype", "sex", "temperature", "treatment"),
                     cores = 8){
  
  genes <- processGenes(genes,taxa)
  
  experiments <- taxa %>% 
    mclapply(function(t){
      vsmSearch(genes[taxon == t, entrez.ID],
             taxa = t,
             confounds = confounds,
             filter = NULL,
             mfx = multifunctionality,
             geeq = geeq,
             p_threshold = p_threshold)
    },cores = cores)
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
  
  
  studies <- length(unique(tmp$ee.ID))
  assays <- tmp[!duplicated(ee.ID), sum(ee.NumSample)]
  n_exp <- length(exps)
}
