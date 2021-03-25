# Visualizations
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders) # From jsicherman/shinycssloaders, NOT daattali
library(htmlwidgets)
library(DT)
library(heatmaply)
library(shinyHeatmaply)
library(shinypanels) # From jsicherman/shinypanels, NOT datasketch
library(circlepackeR)
library(d3wordcloud)
library(data.tree)
# library(sparkline)
library(RColorBrewer)

library(async)
library(memoise)

# Data drivers
library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)
library(stringr)
library(bit)

# Parsing helpers
library(gemmaAPI, lib.loc = '/home/omancarci/R/x86_64-redhat-linux-gnu-library/3.6/')
library(ermineR)
library(mygene)
library(homologene)
library(jsonlite)
library(XML)

library(parallel)
library(lhs)

source('dependencies.R')

DATA.HOLDER$artificial <- NULL
DATA.HOLDER$mouse <- NULL
DATA.HOLDER$rat <- NULL
CACHE.BACKGROUND$artificial <- NULL
CACHE.BACKGROUND$mouse <- NULL
CACHE.BACKGROUND$rat <- NULL
NULLS$artificial <- NULL
NULLS$mouse <- NULL
NULLS$rat <- NULL

# Spike in scores by doing a search of the human corpus,
# then adding a row at a certain index with a certain SD

# Technically we need to recalculate the CACHE.BACKGROUND and prior
# every time we add a row... But there must be a shortcut...

# The CACHE can be updated by simply adding the appropriate new rows

letterWrap <- function(n, depth = 1) {
  x <- do.call(paste0,
               do.call(expand.grid, args = list(lapply(1:depth, function(x) return(LETTERS)), stringsAsFactors = F)) %>%
                 .[, rev(names(.[])), drop = F])
  
  if(n <= length(x)) return(x[1:n])
  
  return(c(x, letterWrap(n - length(x), depth = depth + 1)))
}

SUFFIX <- '_lhs2'
options(mc.cores = 3)

mContrasts <- DATA.HOLDER$human@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique
mGraph <- simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(ChildNode_Long, ParentNode_Long)]))
graphTerms <- unique(ONTOLOGIES[, as.character(ChildNode_Long, ParentNode_Long)])

if(!exists('hypercube')) {
  print('Making hypercube')
  hypercube <- improvedLHS(5000, 4)
  hypercube <- cbind(1, hypercube) #hypercube[, 1] <- pmin(hypercube[, 5], 1L + as.integer(hypercube[, 1] * 20)) # Groups
  hypercube[, 2] <- 1L + as.integer(hypercube[, 2] * 99) # Genes
  hypercube[, 4] <- 1L + as.integer(hypercube[, 4] * 100) # Mean SD
  hypercube[, 5] <- 1L + as.integer(hypercube[, 5] * 49) # Experiments
}

mclapply(1:nrow(hypercube), function(iter) {
  print(iter)
  genes <- DATA.HOLDER$human@gene.meta[sample(1:nrow(DATA.HOLDER$human@gene.meta), hypercube[iter, 2]), entrez.ID]
  tmp <- search(genes)
  
  if(is.null(tmp)) {
    message('No ranking on genes')
    return(NULL)
  }
  
  # Decide where to insert
  possible_indices <- 1 + as.integer(hypercube[iter, 3] * nrow(tmp)) %>% {
    max(0, . - hypercube[iter, 5] * hypercube[iter, 4] - 1):min(. + hypercube[iter, 5] * hypercube[iter, 4] + 1, nrow(tmp))
  } %>% pmin(nrow(tmp) - 1)
  
  indices <- tryCatch({
    sample(possible_indices, hypercube[iter, 5], F, dnorm(possible_indices, hypercube[iter, 3] * nrow(tmp), hypercube[iter, 4]))
  }, error = function(e) {
    message('Resolving too few probs')
    sample(possible_indices, hypercube[iter, 5], F, dnorm(possible_indices, hypercube[iter, 3] * nrow(tmp), hypercube[iter, 4]) + 1)
  })
  
  scores <- tmp[indices, score]
  
  # Insert experiments
  nExp <- letterWrap(hypercube[iter, 5])
  tmp <- rbind(tmp, lapply(1:length(indices), function(i) {
    tmp[indices[i]] %>% copy %>% .[, rn := nExp[i]]
  }) %>% rbindlist) %>% setorder(-score)
  
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
  
  # TODO compute a new prior?
  
  # Enrich
  tmp2 <- enrich(tmp, CACHE = CACHE)
  
  list(mean_index = hypercube[iter, 3],
       mean_sd = hypercube[iter, 4],
       genes = as.integer(genes),
       n_exp = hypercube[iter, 5],
       groups = mGroups,
       associations = mAssoc,
       indices = indices,
       maxIndex = nrow(tmp),
       scores = scores,
       enrich = tmp2 %>% .[, I := .I] %>% .[, f := I / max(I)] %>%
         merge(CACHE$human[, !c('reverse', 'distance')],
               by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri')) %>%
         .[rsc.ID %in% nExp,
           .(A, B, C, stat, distance, I, f,
             group = as.integer(as.character(factor(rsc.ID, levels = nExp, labels = mGroups[mAssoc]))))] %>%
         unique %>% setorder(I))
}) %>% saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/artificial',
                      SUFFIX, '/bootstrap_scores.rds'))
