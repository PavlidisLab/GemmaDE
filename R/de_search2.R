# a simplified version of DEsearch allowing more modularity
# and hopes to be easier to read
# always single taxon input, always use pre-computed scores
# always use-precomputed scores as nulls
# always return a data.table of constant size for a species
# do not rely on external configuration

# any weighting of the scores should be performed via pre-calculation

de_search2 = function(genes = NULL,
                      taxon = NULL,
                      score_source = 'scores',
                      max_dist,
                      confounds,
                      p_threshold){
  
  genes = processGenes(genes)
  
}
  


vsmSearch2 = function(genes,
                      taxon,
                      confounds){
  
  
}

normalize2 = function(){
  
}


# use a matrix with contrasts as columns to calculate the null distribution
generate_nulls = function(scores_contrast = "score_contrast",max_cores = 8){
  tax = DATA.HOLDER %>% names
  
  parallel::mclapply(tax, function(t){
    dt = DATA.HOLDER[[t]]@data[[scores_contrast]]
    seq_len(ncol(dt)) %>% lapply(function(i){
      print(i)
      vals = dt[,i] %>% na.omit()
      data.table(rn = dimnames(dt)[[2]][i],
                 score.mean = mean(vals),
                 score.sd = sd(vals))
    }) %>% do.call(rbind,.)  -> nulls
      
  },mc.cores = min(length(tax),max_cores)
  )-> out
  
  names(out) = tax
  return(out)
}



