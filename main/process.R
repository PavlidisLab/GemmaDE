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
#' @param options Optional extra parameters to pass, such as:
#' * pv: A p-value cutoff (default: 0.05)
#' * fc.lower / fc.upper: Upper and lower logFC thresholds (default: 0 / 10)
#' * mfx: Whether or not to scale by gene multifunctionality
#' @param session The Shiny session
search <- function(genes, taxa = getOption('app.taxa'), options = getOption('app.all_options'), session = NULL) {
  advanceProgress(session, 'Loading experiments')
  
  P_THRESHOLD <- options$pv
  FC_L_THRESHOLD <- options$fc.lower
  FC_U_THRESHOLD <- options$fc.upper
  DIR <- options$dir
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
  # TODO For some reason these aren't the same length?
  colFilter <- colSums2(pv < P_THRESHOLD, na.rm = T) > 0 &
    DATA.HOLDER[[taxa]]@experiment.meta[, !is.na(cf.BaseLongUri) | !is.na(cf.ValLongUri)]
  pv <- pv[, colFilter]
  
  if(n.genes == 1)
    pv <- t(pv)
  
  pv <- pv %>% as.data.table
  
  # Number of DEs for experiments that pass thresholds
  geeq <- DATA.HOLDER[[taxa]]@experiment.meta$ee.qScore[colFilter]
  n.DE.exp <- DATA.HOLDER[[taxa]]@experiment.meta$n.DE[colFilter]
  n.DE.gen <- DATA.HOLDER[[taxa]]@gene.meta$n.DE[rowFilter]
  
  # logFCs for only the GOI/EOI and maintain structure.
  logFC <- DATA.HOLDER[[taxa]]@data$fc[rowFilter, colFilter]
  logFC[is.na(logFC)] <- 0
  
  if(DIR == 'Ignore')
    logFC <- abs(logFC)
  
  if(n.genes == 1)
    logFC <- t(logFC)
  
  logFC <- logFC %>% as.data.table
  
  # Only retain experiments that have at least one of the GOI with threshold_u > abs(logFC) > threshold_l
  if(FC_L_THRESHOLD != 0 || FC_U_THRESHOLD != 0) {
    exFilter <- logFC > FC_L_THRESHOLD
    
    if(FC_U_THRESHOLD != FC_L_THRESHOLD)
      exFilter <- exFilter & (logFC < FC_U_THRESHOLD)
    
    colFilter <- colSums2(exFilter) > 0
    logFC <- logFC %>% .[, colFilter, with = F]
    pv <- pv[, colFilter, with = F]
    n.DE.exp <- n.DE.exp[colFilter]
    geeq <- geeq[colFilter]
  }
  
  if(DIR != 'Ignore') {
    # Force all entries to be positive
    if(any(logFC < 0))
      logFC <- logFC + abs(min(logFC))
    
    # Flip the matrix so maximum values will be the smallest (or largest negative) ones
    if(DIR == 'Downregulated')
      logFC <- max(logFC) - logFC
  }
  
  # Data Processing ---------------------------------------------------------
  advanceProgress(session, 'Ranking')
  
  tf <- t(t(logFC) / (1 + rowMeans2(abs(logFC) %>% as.matrix, na.rm = T))) %>% as.data.table
  
  n.DE.exp[is.na(n.DE.exp)] <- 0
  tf <- t(t(tf) / log2(2 + n.DE.exp)) %>% as.data.table # Weight by number of DEGs in this experiment
  
  geeq[is.na(geeq)] <- 0
  tf <- t(t(tf) * (1 + pmax(geeq, 0))) %>% as.data.table # Weight by experiment quality score
  
  pv[is.na(pv)] <- 1
  tf <- tf * -log10(pv) # TODO

  query <- rowMaxs(tf %>% as.matrix)

  # Shrink by MFX /after/ extracting the query
  
  if(MFX)
    tf <- tf * (1 - DATA.HOLDER[[taxa]]@gene.meta$mfx.Rank[rowFilter] / 5) # TODO Fine tune? Effect similar to idf?
  
  # TODO Should we do IDF on a corpus-level or subset level? Should it recompute n.DE for the p threshold?
  idf <- log10(length(DATA.HOLDER[[taxa]]@gene.meta$n.DE) / (1 + n.DE.gen)) + 1
  # idf <- log10(ncol(tf) / (1 + rowSums2(tf != 0))) + 1
  tf[, query := query]
  tfidf <- tf * idf
  
  if(n.genes > 1)
    tfidf <- scale(tfidf, center = F, scale = colSums2(as.matrix(tfidf))) %>% as.data.table
  
  query <- tfidf[, query]
  tfidf[, query := NULL]
  
  tfidf <- as.matrix(tfidf)
  scores <- query %*% tfidf
  
  # Gene scores are meaningless for single gene queries and scores need scaling.
  if(n.genes == 1) {
    tfidf <- tfidf / tfidf
    scores <- scores / max(scores)
  }
  
  rbind(tfidf, scores) %>% t %>% as.data.table %>%
    `colnames<-`(c(DATA.HOLDER[[taxa]]@gene.meta$gene.Name[rowFilter], 'score')) %>%
    `rownames<-`(colnames(tfidf)) %>% setorder(-score)
}

#' getTags
#' 
#' Expand the specified ontology for the specified experiments.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
#' @param scope The ontology scope.
#' @param rsc.IDs A list of experiment rsc IDs or NULL for everything
#' @param session The Shiny session
getTags <- function(taxa = getOption('app.taxa'), scope = getOption('app.ontology'), rsc.IDs = NULL, session = NULL) {
  if(!exists('CACHE.BACKGROUND'))
    precomputeTags(taxa)
  
  if(is.null(rsc.IDs))
    rsc.IDs <- CACHE.BACKGROUND[[taxa]][, unique(rsc.ID)]
  
  CACHE.BACKGROUND[[taxa]] %>%
    .[rsc.ID %in% rsc.IDs & tag %in% ONTOLOGIES[OntologyScope %in% scope, as.character(ChildNode_Long, ParentNode_Long)]] %>%
    merge(ONTOLOGIES.DEFS[OntologyScope %in% scope], by.x = 'tag', by.y = 'Node_Long') %>%
    .[, .(distance = mean(distance)), .(rsc.ID, Definition)] %>% setorder(distance, rsc.ID)
}

#' Precompute Tags
#' 
#' Get all the expanded ontology tags within a given taxon.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
precomputeTags <- function(taxa = getOption('app.taxa')) {
  graph <- simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(ChildNode_Long, ParentNode_Long)]))
  graphTerms <- unique(ONTOLOGIES[, as.character(ChildNode_Long, ParentNode_Long)])
  
  availableMeta <- DATA.HOLDER[[taxa]]@experiment.meta[!is.na(cf.BaseLongUri) | !is.na(cf.ValLongUri)]
  
  availableTags <- availableMeta[, list(tags = Filter(function(tag) !is.na(tag) & tag %in% graphTerms,
                                                      unique(c(parseListEntry(as.character(cf.BaseLongUri)),
                                                               parseListEntry(as.character(cf.ValLongUri)))))), rsc.ID]
  
  expandTags <- function(entries) {
    rbindlist(lapply(entries, function(entry) {
      to <- igraph::subcomponent(graph, entry, 'out')
      data.table(tag = names(to), distance = c(igraph::distances(graph, entry, to)))
    }))
  }
  
  availableTags[, expandTags(.SD), rsc.ID] %>% .[tag != 'NA' & is.finite(distance)]
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
enrich <- function(rankings, taxa = getOption('app.taxa'), scope = getOption('app.ontology'), session = NULL) {
  rankings <- data.table(rsc.ID = rownames(rankings), rank = rankings$score)
  
  # mMaps[[1]] is the background tags, mMaps[[2]] is the foreground tags.
  mMaps <- list(getTags(taxa, scope, NULL, session), getTags(taxa, scope, rankings[, rsc.ID], session))
  
  mMaps[[2]] <- mMaps[[2]] %>% merge(rankings, by = 'rsc.ID', sort = F)
  
  tmp0 <<- rankings
  tmp1 <<- mMaps
  
  advanceProgress(session, 'Wrapping up')
  
  aprior <- function() {
    merge(mMaps[[2]][, .(N.x = .N, N.ranked = sum(rank), distance = median(distance)), Definition][, N.tags := sum(N.ranked)],
          mMaps[[1]][!(rsc.ID %in% rankings[, rsc.ID]), .(N.y = .N), Definition][, N.bg := sum(N.y)],
          by = 'Definition', sort = F)
  }
  
  fisher <- function(input) {
    input[, c('pv.fisher', 'OR') :=
            suppressWarnings(fisher.test(matrix(c(N.ranked, N.tags - N.ranked,
                                                  N.y, N.bg - N.y),
                                                nrow = 2), alternative = 'greater', conf.int = F))[c('p.value', 'estimate')], Definition]
  }
  
  chisq <- function(input) {
    input[, c('pv.chisq', 'chisq') := ifelse((N.ranked * (N.bg - N.y)) / ((N.tags - N.ranked) * N.y) <= 1, list(1, 1),
                                             suppressWarnings(chisq.test(matrix(c(N.ranked, N.tags - N.ranked,
                                                                                  N.y, N.bg - N.y),
                                                                                nrow = 2)))[c('p.value', 'statistic')]), Definition]
  }

  aprior() %>% chisq %>% fisher %>% .[!(Definition %in% BLACKLIST)] %>% setorder(pv.chisq)
}
