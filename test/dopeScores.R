source('/home/jsicherman/Thesis Work/requirements.R')
source('/home/jsicherman/Thesis Work/dependencies.R')

library(lhs)
library(parallel)
options(mc.cores = 14)

mContrasts <- DATA.HOLDER$human@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique
mGraph <- simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(ChildNode_Long, ParentNode_Long)]))
graphTerms <- unique(ONTOLOGIES[, as.character(ChildNode_Long, ParentNode_Long)])

DATA.HOLDER[c('artificial', 'mouse', 'rat')] <- NULL
CACHE.BACKGROUND[c('artificial', 'mouse', 'rat')] <- NULL
NULLS[c('artificial', 'mouse', 'rat')] <- NULL

if(!exists('hypercube')) {
  print('Making hypercube')
  hypercube <- improvedLHS(6000, 3) # Groups | Genes | Experiments | Percentile
  hypercube <- cbind(1, hypercube) # Groups
  hypercube[, 2] <- 1L + as.integer(hypercube[, 2] * 99) # Genes
  hypercube[, 3] <- 1L + as.integer(hypercube[, 3] * 49) # Experiments
}

mclapply(1:nrow(hypercube), function(iter) {
  print(paste0('Iter ... ', iter, ' ... starting'))
  genes <- DATA.HOLDER$human@gene.meta[sample(1:nrow(DATA.HOLDER$human@gene.meta), hypercube[iter, 2]), .(gene.Name, entrez.ID)]
  tmp <- search(genes$entrez.ID)
  
  if(is.null(tmp))
    return(NULL)
  
  tryCatch({
  invisible(lapply(names(tmp[, !c('rn', 'score', 'f.IN', 'f.OUT', 'ee.q')]), function(.name) set(tmp, which(is.infinite(tmp[[.name]])), j = .name, value = NA)))
  
  # Insert experiments
  nExp <- letterWrap(hypercube[iter, 3])
  
  # Select contrast groups
  mGroups <- sample(1:nrow(mContrasts), hypercube[iter, 1])
  pGroups <- mContrasts[mGroups]
  mAssoc <- sample(1:hypercube[iter, 1], hypercube[iter, 3], T)
  
  # Make new experiments
  eMeta <- new('EData', taxon = 'spikein', data = list(fc = matrix(), adj.pv = matrix()),
               experiment.meta = data.table(rsc.ID = nExp,
                                            ee.ID = max(DATA.HOLDER$human@experiment.meta$ee.ID) + 1:hypercube[iter, 3],
                                            cf.Cat = pGroups[mAssoc, cf.Cat],
                                            cf.BaseLongUri = pGroups[mAssoc, cf.BaseLongUri],
                                            cf.ValLongUri = pGroups[mAssoc, cf.ValLongUri]),
               gene.meta = data.table(),
               go = data.table())
  
  # Make new cache entries
  CACHE <- CACHE.BACKGROUND %>% copy
  CACHE$human <- rbind(CACHE$human,
                       precomputeTags(mGraph = mGraph,
                                      graphTerms = graphTerms,
                                      DATA = list(human = eMeta), POST = F)) %>%
    reorderTags2 %>%
    
    # And use some heuristics to make our corpus a little more lean
    .[, N := length(unique(ee.ID)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    .[, ees := paste0(unique(ee.ID), collapse = ''), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    
    # For every chain (beginning with each rsc.ID, distance == 0), only keep the minimum distance version
    # per unique ees
    .[, .SD[1], .(rsc.ID, ees)] %>%
    
    .[, .(rsc.ID, ee.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri, reverse, distance)]
  
  # Enrich
  tmp2 <- enrich(tmp %>% normalize %>% {
    ind <- max(1, hypercube[iter, 4] * nrow(.)) %>% ceiling
    insertion <- .[sample(max(1, ind - floor(hypercube[iter, 3] / 2)):min(nrow(.), ind + floor(hypercube[iter, 3] / 2)), hypercube[iter, 3], T)] %>% copy
    rbind(., insertion[, rn := nExp])
  }, doNorm = F, CACHE = CACHE, options = getConfig(categories = getConfig('categories')$choices %>% unlist %>% unname))
  
  tmp3 <- enrich(tmp, CACHE = CACHE, options = getConfig(categories = getConfig('categories')$choices %>% unlist %>% unname))
  
  delta <- merge(tmp2[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, score, f = .I / max(.I))],
                 tmp3[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, score, f = .I / max(.I))],
                 by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri')) %>% .[, c('scoreDelta', 'fDelta') := list(score.x - score.y,
                                                                                                                f.x - f.y)]
  
  # Groups | Genes | Dropout | SD from mean | Experiments | f.OUT
  list(genes = as.integer(genes$entrez.ID),
       n_exp = hypercube[iter, 3],
       insert = hypercube[iter, 4],
       groups = mGroups,
       associations = mAssoc,
       enrich = tmp2 %>%
         .[, f := .I / max(.I)] %>%
         .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, distance, score, f)] %>%
         .[cf.Cat %in% CACHE$human[rsc.ID %in% nExp, cf.Cat] &
             cf.BaseLongUri %in% CACHE$human[rsc.ID %in% nExp, cf.BaseLongUri] &
             cf.ValLongUri %in% CACHE$human[rsc.ID %in% nExp, cf.ValLongUri]],
       diff = delta %>%
         .[cf.Cat %in% CACHE$human[rsc.ID %in% nExp, cf.Cat] &
             cf.BaseLongUri %in% CACHE$human[rsc.ID %in% nExp, cf.BaseLongUri] &
             cf.ValLongUri %in% CACHE$human[rsc.ID %in% nExp, cf.ValLongUri]] %>%
         .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, scoreDelta, fDelta)])
  }, error = function(e) {
    print(e)
    NULL
  })
}) %>% saveRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/bootstrapped_scores.rds')
