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
                    geeq = FALSE,
                    categories = c("age", "behavior", "biological process", "biological sex", 
                                   "cell type", "clinical history", "diet", "disease", "environmental history", 
                                   "environmental stress", "genotype", "medical procedure", "molecular entity", 
                                   "organism part", "phenotype", "sex", "temperature", "treatment"),
                    filter_stopwords = TRUE,
                    remove_experiments = NULL,
                    remove_comparisons = NULL,
                    na = 0,
                    cores = 8){
  
  p_genes = process_genes(genes,taxon)
  
  # setting up the experiment filter
  exp_filter = DATA.HOLDER[[taxon]]@experiment.meta$ee.Name %in% remove_experiments
  comp_filter = DATA.HOLDER[[taxon]]@experiment.meta$rsc.ID %in% remove_comparisons
  confound_filter = !(DATA.HOLDER[[taxon]]@experiment.meta$ef.IsBatchConfounded | is.na(DATA.HOLDER[[taxon]]@experiment.meta$ef.IsBatchConfounded))
  
  filter = (!(exp_filter | comp_filter)) & confound_filter
  
  # get the scores. these are precalculated. for a new method of score calculation
  # simply generate a new matrix
  c_scores = get_contrast_scores(genes = p_genes$entrez.ID,
                      taxon = taxon,
                      score_source = score_source,
                      filter = filter,
                      na = na)
  
  # get z scores of the scores comparing them to the null set
  cz_scores  = score_z_scores(c_scores,
                              taxon = taxon,
                              score_source=score_source)
  

  
  # fetches the annotation data for experiments and removes for unwanted categories
  # original getTags will have to be modified after the next data update due to
  # the new fixes
  # i have moved this part out of the functions to allow for better consistency
  # we previously had problems with different invocations of getTags returning different
  # results due to filtering for different things
  experiment_tags = get_tags(taxon, c_scores$rn, max_dist, filter = filter_stopwords,categories = categories) 
  
  # now we do the wilcox test based on the returned tags. p values are calculated for each group
  wc_scores <- wilcox_enrich(cz_scores,experiment_tags,cores)
  
  # normalization factor per experiment, this is the ee.q * (1 + f.IN)/(1 + 10^f.OUT) bit\
  # replace or remove this function to alter normalization
  normalization_factors = get_normalization_factors(
    genes = p_genes$entrez.ID,
    taxon = taxon,
    geeq = geeq,
    p_threshold =p_threshold,
    filter = filter
  )
  normalization_fs = normalization_factors$factor
  names(normalization_fs) = normalization_factors$rsc.ID
  
  effect_size = matrixStats::rowSums2(as.matrix(cz_scores[,p_genes$gene.Name,with = FALSE]),na.rm = TRUE)
  # normalizing effect sizes. I have misread the original code. this never actually happens. may or may not be intentional
  # effect_size = effect_size*normalization_fs[as.character(cz_scores$rn)] # reordering of normalization shouldn't be needed but just in case..
  names(effect_size) = cz_scores$rn
  # a generic function to apply a function on groups
  group_effect_size = apply_on_group(effect_size,experiment_tags,mean)
  


  test_statistic = matrixStats::rowSums2(as.matrix(wc_scores[,p_genes$gene.Name,with = FALSE]))/nrow(p_genes)
  names(test_statistic) = wc_scores$grouping
  
  # normalizing test statistics
  # median of normalization factors are used to weight the test stat
  median_normalization_fs = apply_on_group(normalization_fs,experiment_tags,median)
  test_statistic = test_statistic*median_normalization_fs$vector[match(wc_scores$grouping,median_normalization_fs$grouping)] # reordering of normalization shouldn't be needed but just in case..
  
  
  # evidence counts and other information about the groups
  info = gather_details(experiment_tags,taxon)
  
  
  # compiling the final output
  browser()
  
  out = data.table(
    `Condition Comparison` =  paste0( info$cf.Base, " vs. ", info$cf.Val),
    Category = info$cf.Cat,
    cf.Base = info$cf.Base,
    cf.Val = info$cf.Val,
    Baseline = info$cf.Base,
    Value = info$cf.Val,
    `Ontology Steps` = info$distance,
    conrast_n = info$contrast_n,
    experiment_n = info$experiment_n,
    evidence = info$evidence,
    evidence_gse = info$Evidence,
    `Test Statistic` = test_statistic[as.character(info$grouping)], # reordering shouldn't be needed but just in case
    `Effect Size` = group_effect_size$vector[match(info$grouping,group_effect_size$grouping)],
    wc_scores[,p_genes$gene.Name,with = FALSE]
  )
  
  
}


gather_details = function(experiment_tags,taxon){
  
  mData = DATA.HOLDER[[taxon]]
  
  tags = experiment_tags
  info = tags[,.(
    contrast_n = length(ee.ID),
    experiment_n = length(unique(ee.ID)),
    evidence = list(unique(ee.ID)),
    distance = round(mean(distance),2)
  ),
  .(cf.Cat, cf.BaseLongUri, cf.ValLongUri,grouping)]
  
  
  info %<>%  merge(unique(SIMPLIFIED.ONTOLOGY.DEFS[, .(Node_Long = as.character(Node_Long), cf.Base = as.character(Definition))]),
                   by.x = "cf.BaseLongUri",
                   by.y = "Node_Long", 
                   sort = F,
                   allow.cartesian = T, 
                   all.x = T) %>% 
    merge(unique(SIMPLIFIED.ONTOLOGY.DEFS[, .(Node_Long = as.character(Node_Long), cf.Val = as.character(Definition))]),
          by.x = "cf.ValLongUri",
          by.y = "Node_Long",
          sort = F, 
          allow.cartesian = T, 
          all.x = T)
  
  info$cf.Base[is.na(info$cf.Base)] <- info$cf.BaseLongUri[is.na(info$cf.Base)]
  info$cf.Val[is.na(info$cf.Val)] <- info$cf.ValLongUri[is.na(info$cf.Val)]
  
  
  info[,Evidence:= {evidence %>% unlist %>% {mData@experiment.meta$ee.Name[match(.,mData@experiment.meta$ee.ID)]} %>% relist(evidence) %>% sapply(stringi::stri_c, collapse = ',')}]
  
  info %>% setorder(grouping)
  
  
}


# a generic function to calculate group means for any given vector
apply_on_group =function(vector,experiment_tags,fun){
  tags = experiment_tags
  
  tags$vector = vector[as.character(tags$rsc.ID)]
  
  tags[,.(
    vector = fun(vector)
  ),.(cf.Cat, cf.BaseLongUri, cf.ValLongUri,grouping)] %>% data.table::setorder(grouping)
}

# a wrapper for getTags that does the category filtering.
get_tags = function(taxon,experiments,max_dist,filter,categories){
  experiment_tags = 
    getTags(taxon,experiments, max_dist, filter = filter)%>% 
    .[cf.Cat %in% categories]
  
  experiment_tags$grouping =  paste(experiment_tags$cf.Cat,experiment_tags$cf.BaseLongUri,experiment_tags$cf.ValLongUri)
  return(experiment_tags)
}

# originally, enrich function. runs the wilcox test
# note that p values test for >0 so calculation of z scores before this
# is required
wilcox_enrich = function(cz_scores,
                         experiment_tags,
                         cores){
  
  terms = cz_scores %>% merge(experiment_tags,by.x='rn',by.y = 'rsc.ID',sort = FALSE)
  
  core_split = terms$grouping %>% factor %>% as.numeric() %>% {.%%cores}
  terms[,core_split := core_split]
  
  gene_names = colnames(cz_scores) %>% {.[!. %in% 'rn']}
 
  
  grouped = terms %>% data.table:::split.data.table(by = c('core_split'))
  
  
  grouped %>% parallel::mclapply(function(t){
    t[, matrixTests::col_wilcoxon_onesample(as.matrix(.SD), alternative = "greater", exact = F) %>%
        {
          list(pv = 1 - .[, "pvalue"], gene = rownames(.)) # Small p-value if real effect (= score closer to 1)
        }, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri,grouping),
      .SDcols = gene_names
    ] %>%
      data.table::dcast(... ~ gene, value.var = "pv", fill = 0)
  },mc.cores = cores) %>% do.call(rbind,.) %>%setorder(grouping)  -> wilcox_ps
  
  return(wilcox_ps)
   
}

# return the normalization factors, returns a named array
# can be replaced with a function of equivalent output or an array of 1s
# to remove all normalization
get_normalization_factors = function(genes,
                                     taxon,
                                     geeq,
                                     p_threshold,
                                     filter){
  
  #  getOption("app.algorithm.experiment")
  
  mData = DATA.HOLDER[[taxon]]
  experimentMeta = mData@experiment.meta[filter,.(rsc.ID, ee.qScore, n.DE, ad.NumGenes,nonNa.numGenes)]
  
  geneMask <- which(mData@gene.meta$entrez.ID %in% genes)

  mDimNames <- dimnames(mData@data$adj.pv)
  pv <- mData@data$adj.pv[, geneMask, drop = F] %>%
    `dimnames<-`(list(mDimNames[[1]], mDimNames[[2]][geneMask])) %>%
    .[filter, , drop = F] %>%
    t()
  
  pv[is.na(pv)] = 1
  experimentMeta$n.DE[is.na(experimentMeta$n.DE)] = 0
  
  experimentN = Rfast::colsums(pv <= p_threshold)
  
  
  if (!geeq) {
    experimentMeta$ee.qScore <- 2
  } else {
    experimentMeta$ee.qScore <- pmax(experimentMeta$ee.qScore + 1, 0, na.rm = T)
  }
  experimentMeta$ee.qScore <- sqrt(experimentMeta$ee.qScore / 2)
  
  normalization_factors = data.table::data.table(rsc.ID =  experimentMeta$rsc.ID,
                                                 f.IN = experimentN/length(geneMask),
                                                 f.OUT =  pmax(0, experimentMeta$n.DE - experimentN)/ experimentMeta$nonNa.numGenes,
                                                 ee.q = experimentMeta$ee.qScore )
  
  normalization_factors[,factor := ee.q * (1 + f.IN)/(1 + 10^f.OUT)]

  return(normalization_factors)
}
  

# mostly equivalent to old vsmSearch
# with the calculation pre-calculated, primary function of this is to simply
# subset the data. it might be prudent to add some option to weight the scores
# but if added before normalization, it'll add a substantial overhead
get_contrast_scores = function(genes,
                               taxon,
                               score_source,
                               filter,
                               na = 0){
  
  mData = DATA.HOLDER[[taxon]]
  scores = mData@data[[score_source]]
  geneMask =  which(mData@gene.meta$entrez.ID %in% genes)
  subset = scores[filter,geneMask,drop = FALSE] %>% data.table
  subset[is.na(subset)] = na
  colnames(subset) = mData@gene.meta$gene.Name[geneMask]
  # subset[, score := Rfast::rowsums(as.matrix(subset))] # this sum isn't useful
  subset[, rn := mData@experiment.meta$rsc.ID[filter]]
  return(subset)
}



# mostly the old normalize without the normalization factor
# only calculates the z scores based on pre-calculated nulls
# for a data matrix. if the chosen matrix doesn't already have 
# a nullset, it'll be created. while this isn't very slow
# precalculation is good for you
score_z_scores = function(c_scores, taxon, score_source = 'score'){
  if(score_source %in% names(data_nulls)){
    nulls = data_nulls[[score_source]][[taxon]]
  } else{
    warning('nulls not found, recalculating. consider re-running package_data.R to generate them')
    nulls = generate_nulls(paste0(score_source,'_contrast'))
  }
  
  merged = c_scores %>% merge(nulls,by='rn',sort = F)
  merged %>%  .[, lapply(.SD, function(x) (x - score.mean) / score.sd),
                .SDcols = !c("score.mean", "score.sd", "rn")
  ]%>%
    {
      data.table::data.table(rn = merged$rn,.)
    }
  
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
