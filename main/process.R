library(igraph)
library(data.table)
library(dplyr)
library(matrixStats)

#' Search
#' 
#' Uses the M-VSM to sort experiments that show at least one of the genes as DE.
#'
#' @param genes A list of Entrez Gene IDs (ie. 1, 22, 480) as characters
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
#' @param scope The ontology scope, or null. If null, no weighting by experiment prior is done
#' @param options Optional extra parameters to pass, such as:
#' * pv: A p-value cutoff (default: 0.05)
#' * fc.lower / fc.upper: Upper and lower logFC thresholds (default: 0 / 10)
#' * mfx: Whether or not to scale by gene multifunctionality
#' @param session The Shiny session
search <- function(genes, taxa = 'human', scope = NULL, options = DEFAULT_OPTIONS, session = NULL) {
  advanceProgress(session, 'Loading experiments')
  
  P_THRESHOLD <- options$pv
  FC_L_THRESHOLD <- options$fc.lower
  FC_U_THRESHOLD <- options$fc.upper
  MFX <- options$mfx
  
  # Data Extraction ---------------------------------------------------------

  # Only retain GOI
  rowFilter <- which(DATA.HOLDER[[taxa]]@gene.meta$entrez.ID %in% genes)
  n.genes <- length(rowFilter)
  if(n.genes == 0) return(NULL)
  
  # P-values for only the GOI
  pv <- DATA.HOLDER[[taxa]]@data$adj.pv[rowFilter, ]
  
  if(n.genes == 1)
    pv <- t(pv)
  
  # Only retain experiments that have at least one of the GOI as DE (pv < threshold)
  colFilter <- colSums2(pv < P_THRESHOLD, na.rm = T) > 0
  pv <- pv[, colFilter] %>% as.data.table
  
  # Number of DEs for experiments that pass thresholds
  n.DE <- DATA.HOLDER[[taxa]]@experiment.meta$n.DE[colFilter]
  
  # logFCs for only the GOI/EOI and maintain structure.
  logFC <- DATA.HOLDER[[taxa]]@data$fc[rowFilter, colFilter]
  logFC[is.na(logFC)] <- 0
  logFC <- logFC %>% abs %>% as.data.table
  
  # Only retain experiments that have at least one of the GOI with threshold_u > abs(logFC) > threshold_l
  if(FC_L_THRESHOLD != 0 || FC_U_THRESHOLD != 0) {
    exFilter <- logFC > FC_L_THRESHOLD
    
    if(FC_U_THRESHOLD != FC_L_THRESHOLD)
      exFilter <- exFilter & (logFC < FC_U_THRESHOLD)
    
    colFilter <- colSums2(exFilter) > 0
    logFC <- logFC %>% .[, colFilter, with = F]
    pv <- pv[, colFilter, with = F]
    n.DE <- n.DE[colFilter]
  }
  
  # Data Processing ---------------------------------------------------------
  advanceProgress(session, 'Ranking')
  
  tf <- logFC
  tf <- t(t(tf) / log2(1 + n.DE)) %>% as.data.table # Weight by number of DEs
  
  pv[is.na(pv)] <- 1
  tf <- tf * -log10(pv) # TODO
  
  query <- rowMaxs(tf %>% as.matrix)
  
  # Shrink by MFX /after/ extracting the query
  
  if(MFX)
    tf <- tf * (1 - DATA.HOLDER[[taxa]]@gene.meta$mfx.Rank[rowFilter] / 4) # TODO
  
  idf <- log10(ncol(tf) / (1 + rowSums2(tf != 0))) + 1
  tf[, query := query]
  tfidf <- tf * idf
  
  if(nrow(tfidf) > 1)
    tfidf <- scale(tfidf, center = F, scale = colSums2(as.matrix(tfidf))) %>% as.data.table
  
  query <- tfidf[, query]
  tfidf[, query := NULL]
  
  tfidf <- as.matrix(tfidf)
  scores <- query %*% tfidf
  
  mData <- rbind(tfidf, scores) %>% t %>% as.data.table
  colnames(mData) <- c(DATA.HOLDER[[taxa]]@gene.meta$gene.Name[rowFilter], 'score')
  rownames(mData) <- colnames(tfidf)
  mData %>% setorder(-score)
}

#' getTags
#' 
#' Expand the specified ontology for the specified experiments.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
#' @param scope The ontology scope.
#' @param rsc.IDs A list of experiment rsc IDs. Runs each list item separately and collapses into a single data table.
#' @param session The Shiny session
getTags <- function(taxa = 'human', scope = 'DO', rsc.IDs = list(NULL), session = NULL) {
  ontoScope <- ONTOLOGIES[OntologyScope %in% scope, .(ChildNode_Long, ParentNode_Long)]
  graph <- igraph::graph_from_data_frame(ontoScope)
  
  tags <- DATA.HOLDER[[taxa]]@experiment.meta[!is.na(ee.TagLongUri) |
                                                      !is.na(cf.ValLongUri) |
                                                      !is.na(cf.BaseLongUri)]
  
  allTerms <- ontoScope[, ChildNode_Long]
  lapply(1:length(rsc.IDs), function(indx) {
    IDs <- rsc.IDs[[indx]]

    advanceProgress(session, paste0('Calculating scores (', indx, '/', length(rsc.IDs), ')'))
    
    # Subset to experiment IDs if we have them
    if(!is.null(IDs)) mTags <- tags[rsc.ID %in% IDs]
    
    # Otherwise cache background populations
    else {
      # All were already cached
      if(all(scope %in% names(CACHE.BACKGROUND[[taxa]])))
        return(do.call(rbind, lapply(scope, function(id) CACHE.BACKGROUND[[taxa]][[id]])))
      
      # Not all were cached. Cache the ones that aren't
      else if(length(scope) > 1) {
        existing <- intersect(scope, names(CACHE.BACKGROUND[[taxa]]))
        if(length(existing) > 0) # Load cached
          existing <- do.call(rbind, lapply(existing, function(id) CACHE.BACKGROUND[[taxa]][[id]]))
        else
          existing <- NULL
        
        # Make sure to make `session` NULL so we don't generate extra loading bar steps.
        # TODO Loading bar won't reflect extra caches to be done...
        other <- do.call(rbind, lapply(setdiff(scope, names(CACHE.BACKGROUND[[taxa]])),
                                       function(id) getTags(taxa, id, NULL, NULL)))
        
        return(rbind(existing, other))
      }
      
      mTags <- tags
    }
    
    fineProgress(session, nrow(mTags))
    
    do.call(rbind, lapply(1:nrow(mTags), function(indx) {
      advanceProgress(session, fine = T)
      
      filteredTags <- Filter(function(tag) !is.na(tag) & tag %in% allTerms,
                             c(parseListEntry(mTags[indx, ee.TagLongUri]),
                               parseListEntry(mTags[indx, cf.BaseLongUri]),
                               parseListEntry(mTags[indx, cf.ValLongUri]))) # TODO sf?
      
      do.call(rbind, lapply(filteredTags, function(entry) {
        data.table(rsc.ID = mTags[indx, rsc.ID], tag = names(igraph::subcomponent(graph, entry, 'out')))
      }))
    })) -> mMap
    
    if(is.null(mMap)) return(NULL)
    mMap <- mMap[tag != 'NA']
    
    mData <- unique(ONTOLOGIES.DEFS[OntologyScope %in% scope & Node_Long %in% mMap[, unique(tag)],
                                    Definition, Node_Long]) %>%
      merge(mMap, by.x = 'Node_Long', by.y = 'tag', allow.cartesian = T, sort = F) %>%
      .[, .(Definition, rsc.ID)] %>% unique
    
    if(is.null(IDs)) CACHE.BACKGROUND[[taxa]][[scope]] <<- mData
    
    mData
  })
}

#' Enrich
#' 
#' Given rankings (@seealso search), generate a ranking-weighted count of all terms that can be
#' derived from tags present in the experiment (both cf.ValLongUri and ee.TagLongUri).
#'
#' @param rankings A named numeric (@seealso search).
#' @param taxa The taxon
#' @param scope The ontology scope.
#' @param session The Shiny session.
enrich <- function(rankings, taxa = 'human', scope = 'DO', session = NULL) {
  rankings <- data.table(rsc.ID = rownames(rankings), rank = rankings$score)
  
  mMaps <- getTags(taxa, scope, list(NULL, rankings[, rsc.ID]), session)
  
  mMaps[[2]] <- mMaps[[2]] %>% merge(rankings, by = 'rsc.ID', sort = F)
  
  advanceProgress(session, 'Wrapping up')
  
  prior <- function() {
    merge(mMaps[[2]][, .(N.x = .N, N.ranked = sum(rank)), Definition][, N.tags := sum(N.ranked)],
          merge(mMaps[[1]][!(rsc.ID %in% rankings[, rsc.ID])],
                do.call(rbind, lapply(scope, function(id) SUMMARY[[id]]))$experiments[, .(prior = mean(score)), rsc.ID],
                by = 'rsc.ID', sort = F) %>%
            .[, .(N.y = sum(prior)), Definition] %>% .[, N.bg := sum(N.y)],
          by = 'Definition', sort = F)
  }
  
  aprior <- function() {
    merge(mMaps[[2]][, .(N.x = .N, N.ranked = sum(rank)), Definition][, N.tags := sum(N.ranked)],
          mMaps[[1]][!(rsc.ID %in% rankings[, rsc.ID]), .(N.y = .N), Definition][, N.bg := sum(N.y)],
          by = 'Definition', sort = F)
  }
  
  # Test this
  fisher <- function(input) {
    input[, c('pv.fisher', 'OR') :=
            suppressWarnings(fisher.test(matrix(c(N.ranked, N.tags - N.ranked,
                                                  N.x + N.y, N.bg + N.x - N.y),# N.y, N.bg - N.y),
                                                nrow = 2), alternative = 'greater', conf.int = F))[c('p.value', 'estimate')], Definition]
  }
  
  chisq <- function(input) {
    input[, c('pv.chisq', 'chisq') := ifelse((N.ranked * (N.bg + N.x - N.y)) / ((N.tags - N.ranked) * (N.x + N.y)) <= 1, list(1, 1),#(N.ranked * (N.bg - N.y)) / ((N.tags - N.ranked) * N.y) <= 1, 1,
                                             suppressWarnings(chisq.test(matrix(c(N.ranked, N.tags - N.ranked,
                                                                                  N.x + N.y, N.bg + N.x - N.y),#N.y, N.bg - N.y),
                                                                                nrow = 2)))[c('p.value', 'statistic')]), Definition]
  }

  mMaps <- aprior()
  mMap <- mMaps %>% chisq %>% fisher
  
  mMap[!(Definition %in% BLACKLIST)] %>% setorder(pv.chisq)
}
