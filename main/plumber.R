

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
                                    "organism part", "phenotype", "sex", "temperature", "treatment"
                     )){
  
  # delegate the checks to the user function
  genes = processGenes(genes,taxa)
  
  experiments <- taxa %>% 
    lapply(function(t){
      vsmSearch(genes[taxon == t, entrez.ID],
             taxa = t,
             confounds = confounds,
             filter = NULL,
             mfx = multifunctionality,
             geeq = geeq,
             p_threshold = p_threshold)
    })
  
  enrich(experiments, taxa = taxa, dist = max_dist,categories = categories)
    
}
