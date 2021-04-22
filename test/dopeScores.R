source('/home/jsicherman/Thesis Work/requirements.R')

source('dependencies.R')

library(lhs)
library(parallel)
options(mc.cores = 5)

mContrasts <- DATA.HOLDER$human@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique
mGraph <- simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(ChildNode_Long, ParentNode_Long)]))
graphTerms <- unique(ONTOLOGIES[, as.character(ChildNode_Long, ParentNode_Long)])

DATA.HOLDER[c('artificial', 'mouse', 'rat')] <- NULL
CACHE.BACKGROUND[c('artificial', 'mouse', 'rat')] <- NULL
NULLS.EXP[c('artificial', 'mouse', 'rat')] <- NULL

# Spike in scores by doing a search of the human corpus,
# then adding a row at a certain index with a certain SD

# Technically we need to recalculate the CACHE.BACKGROUND and prior
# every time we add a row... But there must be a shortcut...

# The CACHE can be updated by simply adding the appropriate new rows

if(!exists('hypercube')) {
  print('Making hypercube')
  hypercube <- improvedLHS(5000, 4)
  hypercube <- cbind(1, hypercube) #hypercube[, 1] <- pmin(hypercube[, 5], 1L + as.integer(hypercube[, 1] * 20)) # Groups
  hypercube[, 2] <- 1L + as.integer(hypercube[, 2] * 99) # Genes
  hypercube[, 4] <- 1 + hypercube[, 4] * rexp(nrow(hypercube), 1/15) # Mean SD from the mean
  hypercube[, 5] <- 1L + as.integer(hypercube[, 5] * 49) # Experiments
}

mclapply(1:nrow(hypercube), function(iter) {
  print(iter)
  genes <- DATA.HOLDER$human@gene.meta[sample(1:nrow(DATA.HOLDER$human@gene.meta), hypercube[iter, 2]), .(gene.Name, entrez.ID)]
  tmp <- search(genes$entrez.ID, inprod = F)
  
  if(is.null(tmp)) {
    message('No ranking on genes')
    return(NULL)
  }
  
  scores <- matrix(rnorm(hypercube[iter, 2] * hypercube[iter, 5], hypercube[iter, 4]),
                   nrow = hypercube[iter, 5]) %>%
    `*`(matrix(sample(c(rnorm(100, 0, 0.03),
                        rnorm(100, 1, 0.03)), hypercube[iter, 2] * hypercube[iter, 5], T,
                      rep(c(1 - hypercube[iter, 3], hypercube[iter, 3]), each = 100)),
               nrow = hypercube[iter, 5])) %>%
    abs
  
  # Insert experiments
  nExp <- letterWrap(hypercube[iter, 5])
  tmp <- rbind(tmp, data.table(rn = nExp, data.table(scores) %>%
                                 setnames(colnames(tmp[, !c('rn', 'is.passing', 'f.IN', 'f.OUT', 'ee.q')])),
                               is.passing = NA, f.IN = NA, f.OUT = NA, ee.q = NA))
  
  # Select contrast groups
  mGroups <- sample(1:nrow(mContrasts), hypercube[iter, 1])
  pGroups <- mContrasts[mGroups]
  mAssoc <- sample(1:hypercube[iter, 1], hypercube[iter, 5], T)
  
  # Make new experiments
  eMeta <- new('EData', taxon = 'spikein', data = list(fc = matrix(), adj.pv = matrix()),
               experiment.meta = data.table(rsc.ID = letterWrap(hypercube[iter, 5]),
                                            ee.ID = 1:hypercube[iter, 5],
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
    .[, N := .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    .[distance == 0 | (N > 1 & N < 500 & distance < 5)] %>%
    
    # This is annoying but the filter should be applied twice
    .[, N := .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    .[distance == 0 | (N > 1 & N < 500 & distance < 5)] %>%
    
    .[, !'N']
  
  # Enrich
  tmp2 <- enrich(tmp, CACHE = CACHE, inprod = F)
  
  list(mean_index = hypercube[iter, 4],
       dropout = hypercube[iter, 3],
       genes = as.integer(genes$entrez.ID),
       n_exp = hypercube[iter, 5],
       groups = mGroups,
       associations = mAssoc,
       enrich = tmp2 %>%
         .[, I := .I] %>%
         .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, distance, I)] %>%
         .[cf.Cat %in% CACHE$human[rsc.ID %in% nExp, cf.Cat] &
             cf.BaseLongUri %in% CACHE$human[rsc.ID %in% nExp, cf.BaseLongUri] &
             cf.ValLongUri %in% CACHE$human[rsc.ID %in% nExp, cf.ValLongUri]])
}) %>% saveRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/bootstrap_scores.rds')
