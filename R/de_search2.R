# a simplified version of DEsearch allowing more modularity
# and hopes to be easier to read
# always single taxon input, always use pre-computed scores
# always use-precomputed scores as nulls
# always return a data.table of constant size for a species
# do not rely on external configuration

# any weighting of the scores should be performed via pre-calculation

gemma_de = function(genes = NULL,
                    taxon = NULL,
                    score_source = 'score',
                    max_dist = 1.5,
                    confounds = FALSE,
                    p_threshold = 0.05,
                    categories = c("age", "behavior", "biological process", "biological sex", 
                                   "cell type", "clinical history", "diet", "disease", "environmental history", 
                                   "environmental stress", "genotype", "medical procedure", "molecular entity", 
                                   "organism part", "phenotype", "sex", "temperature", "treatment"),
                    filter_stopwords = TRUE,
                    remove_experiments = NULL,
                    remove_comparisons = NULL,
                    cache = NULL,
                    get_descriptions = TRUE, # temporary argument to allow supporting old cache files
                    cores = 8,
                    nullset = NULL){
  
  p_genes = process_genes(genes,taxon)
  
  exp_filter = DATA.HOLDER[[taxon]]@experiment.meta$ee.Name %in% remove_experiments
  comp_filter = DATA.HOLDER[[taxon]]@experiment.meta$rsc.ID %in% remove_comparisons
  confound_filter = !(DATA.HOLDER[[taxon]]@experiment.meta$ef.IsBatchConfounded | is.na(DATA.HOLDER[[taxon]]@experiment.meta$ef.IsBatchConfounded))
  
  filter = (!(exp_filter | comp_filter)) & confound_filter
  get_contrast_scores(genes = p_genes$entrez.ID,
                      taxon = taxon,
                      score_source = score_source,
                      filter = filter)
  
}
  

# mostly equivalent to old vsmSearch
# with the calculation pre-calculated, primary function of this is to simply
# subset the data. it might be prudent to add some option to weight the scores
# but if added before normalization, it'll add a substantial overhead
get_contrast_scores = function(genes,
                               taxon,
                               score_source,
                               filter){
  
  mData = DATA.HOLDER[[taxon]]
  scores = mData@data[[score_source]]
  geneMask =  which(mData@gene.meta$entrez.ID %in% genes)
  subset = scores[filter,geneMask,drop = FALSE] %>% data.table
  
  colnames(subset) = mData@gene.meta$gene.Name[geneMask]
  subset[, score := Rfast::rowsums(as.matrix(subset))]
  subset[, rn := mData@experiment.meta$rsc.ID[filter]]
  return(subset)
}



norm = function(){
  
}



process_genes = function(genes,taxon){

  clean_genes = tidy_genes(genes,taxon)
  
  tax_id = processTaxa(taxon)
  

  symbols = data.table(entrez.ID = as.character(clean_genes),taxon = taxon,taxon.ID = tax_id) %>% 
    merge(., DATA.HOLDER[[taxon]]@gene.meta[, .(.I, entrez.ID = as.character(entrez.ID), gene.Name)], by = "entrez.ID")
  
  
  return(symbols[!duplicated(entrez.ID),])
}

# return a vector of entrez IDs with names set to the original gene list
# integers are processed as entrez IDs
# go terms return matching ids
# ensemble ids, gene names also supported
# possible to return duplicates
tidy_genes = function(genes,taxon){
  oGenes <- genes
  
  # Clean numerics (interpreted as entrez IDs) and remove them from further processing.
  cleanGenes <- suppressWarnings(Filter(function(x) !is.na(as.integer(x)), genes))
  idMap <- sapply(as.character(cleanGenes), function(x) which(oGenes == x), USE.NAMES = F)
  
  genes <- genes[!(genes %in% cleanGenes)]
  
  
  # If it matches (ENSG|ENSMUS|ENSRNO)\d{11}, it's an Ensembl ID (for human, mouse or rat).
  if (length(genes) > 0) {
    ensembl <- grep("(ENSG|ENSMUS|ENSRNO)\\d{11}", genes, value = T)
    
    if (length(ensembl) != 0) {
      # Extract genes with a matching Ensembl ID and clean them too.
      ensembls <- DATA.HOLDER[[taxon]]@gene.meta[ensembl.ID %in% ensembl, .(entrez.ID, ensembl.ID)]
      cleanGenes <- c(cleanGenes, ensembls[, entrez.ID])
      idMap <- c(idMap, sapply(ensembls[, ensembl.ID], function(x) which(oGenes == x), USE.NAMES = F))
      
      genes <- genes[!(genes %in% ensembls[, ensembl.ID])]
    }
  }
  
  # Match GO identifiers
  if (length(genes > 0)) {
    go <- grep("(GO:)\\d{7}", genes, value = T)
    
    if (length(go) != 0) {
      gos <- DATA.HOLDER[[taxon]]@go[id %in% go, .(id, entrez.ID)]
      cleanGenes <- c(cleanGenes, gos[, entrez.ID])
      idMap <- c(idMap, sapply(gos[, id], function(x) which(oGenes == x), USE.NAMES = F))
      
      genes <- genes[!(genes %in% gos[, id])]
    }
  }
  
  # Try to match to gene names
  if (length(genes) > 0) {
    descriptors <- DATA.HOLDER[[taxon]]@gene.meta[gene.Name %in% genes, .(entrez.ID, gene.Name)]
    if (nrow(descriptors) != 0) {
      cleanGenes <- c(cleanGenes, descriptors[, entrez.ID])
      idMap <- c(idMap, sapply(descriptors[, gene.Name], function(x) which(oGenes == x), USE.NAMES = F))
      
      genes <- genes[!(genes %in% descriptors[, gene.Name])]
    }
  }
  
  if(length(genes)>0){
    warning('Unprocessed genes present')
  }
  if(length(cleanGenes)>0){
    cleanGenes <- cleanGenes[order(unlist(idMap))]
    names(cleanGenes) = oGenes[unlist(idMap)]
    return(cleanGenes)
  } else{
    return(character(0))
  }

}
