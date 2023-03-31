# owl parsing
# test case DOID:3363 http://purl.obolibrary.org/obo/DOID_3393
library(xml2)
library(magrittr)
devtools::load_all()

data.holder <- readRDS(file.path(RAWDIR,'DATA.HOLDER.rds'))


# list and fetch ontologies form Gemma's cache -------
ontology_server <- 'moe'
ontology_path <- '/space/gemmaData/ontologyCache/ontology'

ontologies <- system2("ssh" ,paste0(ontology_server," ls ", ontology_path),stdout = TRUE)
ontologies <- ontologies[!grepl('tmp',ontologies)]

dir.create('data-raw/ontologies',showWarnings = FALSE)
ontologies %>% lapply(function(x){
  RCurl::scp(ontology_server,file.path(ontology_path,x)) %>% 
    writeBin(file.path('data-raw/ontologies',x))
})

# temporarily add obi manually until processing is complete
# file is taken from the original source for now
# ontologies %<>% c('obiOntology')


ontologies = ontologies[ontologies %in% list.files('data-raw/ontologies/')]

# get all hasAlternativeId terms from the ontologies to replace within data.holder
alternative_terms <- ontologies %>% lapply(function(x){
  print(x)
  onto <- xml2::read_xml(file.path('data-raw/ontologies/',x))
  children <- xml2::xml_children(onto)
  classes <- xml2::xml_name(children) == 'Class'
  children <- children[classes]
  
  
  term_links <- xml2::xml_attr(children,'about')
  
  alternatives <- children %>% lapply(function(y){
    term_children <- xml2::xml_children(y)
    term_contents <- xml2::xml_name(term_children)
    term_children[term_contents =='hasAlternativeId'] %>% xml2::xml_text()
  })
  
  names(alternatives) <- term_links
  return(alternatives)
})
names(alternative_terms) <- ontologies
alternatives <- alternative_terms[-13] %>% lapply(\(x){
  x %>% sapply(length) %>% rep(names(x),.)
}) %>% unlist
to_replace <- alternative_terms[-13] %>% lapply(\(x){
  x %>% unlist
}) %>% unlist


original_terms <- alternative_terms[-13] %>% lapply(names) %>% unlist %>% unique

data.holder %<>% lapply(function(x){
  val_terms <- x@experiment.meta$cf.ValLongUri %>% {.[grepl('obo',.)]} %>% unique
  base_terms <- x@experiment.meta$cf.BaseLongUri %>% {.[grepl('obo',.)]} %>% unique
  terms <- c(val_terms,base_terms) %>% unique
  missing_terms <- terms[!terms %in% original_terms]
  length(missing_terms)
  compact_names <- missing_terms %>% 
    stringr::str_extract('(?<=obo/).*') %>% 
    stringr::str_replace('_',':')
  

  replacement_terms <- alternatives[match(compact_names,to_replace)]
  names(replacement_terms) <-  missing_terms
  replacement_terms %<>% na.omit()
  
  levels(x@experiment.meta$cf.BaseLongUri) = c(levels(x@experiment.meta$cf.BaseLongUri),replacement_terms) %>% unique
  levels(x@experiment.meta$cf.ValLongUri) = c(levels(x@experiment.meta$cf.ValLongUri),replacement_terms) %>% unique
  
  x@experiment.meta$cf.BaseLongUri[x@experiment.meta$cf.BaseLongUri %in% names(replacement_terms)] %<>% 
    {replacement_terms[match(.,names(replacement_terms))]} 
  x@experiment.meta$cf.BaseLongUri %<>% droplevels()

  x@experiment.meta$cf.ValLongUri[x@experiment.meta$cf.ValLongUri%in% names(replacement_terms)] %<>% 
    {replacement_terms[match(.,names(replacement_terms))]} 
  x@experiment.meta$cf.ValLongUri %<>% droplevels()
  
  return(x)
    
})

saveRDS(data.holder, file.path(RAWDIR, 'DATA.HOLDER_fixed_terms.rds'))


