


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
de_search_plumb = function(req = NULL, # this is the request object
                     genes = NULL,
                     taxa =NULL,
                     max_dist = 1.5,
                     confounds = FALSE,
                     multifunctionality = FALSE,
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
  
  de_search(genes,
            taxa,
            max_dist,
            confounds,
            multifunctionality,
            geeq,
            p_threshold,
            categories,
            cores)
  

}

# de_search(genes=genes,taxa = taxa)
