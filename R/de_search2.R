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





