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
  MFX <- options$mfx
  GEEQ <- options$geeq
  
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
  
  if(n.genes == 1)
    logFC <- t(logFC)
  
  logFC <- logFC %>% as.data.table
  
  directions <- ifelse(colMeans2(logFC %>% as.matrix) < 0, -1, 1)
  logFC <- abs(logFC)
  
  # Only retain experiments that have at least one of the GOI with threshold_u > abs(logFC) > threshold_l
  if(FC_L_THRESHOLD != 0 || FC_U_THRESHOLD != 0) {
    exFilter <- logFC > FC_L_THRESHOLD
    
    if(FC_U_THRESHOLD != FC_L_THRESHOLD)
      exFilter <- exFilter & (logFC < FC_U_THRESHOLD)
    
    colFilter <- colSums2(exFilter) > 0
    logFC <- logFC %>% .[, colFilter, with = F]
    directions <- directions[colFilter]
    pv <- pv[, colFilter, with = F]
    n.DE.exp <- n.DE.exp[colFilter]
    geeq <- geeq[colFilter]
  }
  
  # Data Processing ---------------------------------------------------------
  advanceProgress(session, 'Ranking')

  tf <- t(t(logFC) / (1 + rowMeans2(abs(logFC) %>% as.matrix, na.rm = T))) %>% as.data.table
  
  n.DE.exp[is.na(n.DE.exp)] <- 0
  tf <- t(t(tf) / log2(2 + n.DE.exp)) %>% as.data.table # Weight by number of DEGs in this experiment
  
  geeq[is.na(geeq)] <- 0
  if(GEEQ)
    tf <- t(t(tf) * (1 + pmax(geeq, 0))) %>% as.data.table # Weight by experiment quality score
  
  pv[is.na(pv)] <- 1
  tf <- tf * -log10(pv + 1e-50) # TODO

  query <- matrixStats::rowQuantiles(tf %>% as.matrix, probs = 0.8)# rowMaxs(tf %>% as.matrix)

  # Shrink by MFX /after/ extracting the query
  
  # if(MFX)
  #   tf <- tf * (1 - DATA.HOLDER[[taxa]]@gene.meta$mfx.Rank[rowFilter] / 5) # TODO Fine tune? Effect similar to idf?
  
  # TODO Should we do IDF on a corpus-level or subset level? Should it recompute n.DE for the p threshold?
  idf <- log10(length(DATA.HOLDER[[taxa]]@gene.meta$n.DE) / (1 + DATA.HOLDER[[taxa]]@gene.meta$mfx.Rank[rowFilter])) + 1
  # idf <- log10(length(DATA.HOLDER[[taxa]]@gene.meta$n.DE) / (1 + n.DE.gen)) + 1
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
    .[, direction := directions] %>%
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
  if(is.null(rsc.IDs))
    rsc.IDs <- CACHE.BACKGROUND[[taxa]][, unique(rsc.ID)]
  
  graphTerms <- ONTOLOGIES[!(OntologyScope %in% scope), as.character(ChildNode_Long, ParentNode_Long)]
  
  CACHE.BACKGROUND[[taxa]] %>%
    .[rsc.ID %in% rsc.IDs & !(cf.BaseLongUri %in% graphTerms) & !(cf.ValLongUri %in% graphTerms)] %>%
    merge(ONTOLOGIES.DEFS[OntologyScope %in% scope], by.x = 'cf.BaseLongUri', by.y = 'Node_Long', sort = F, allow.cartesian = T, all.x = T) %>%
    merge(ONTOLOGIES.DEFS[OntologyScope %in% scope], by.x = 'cf.ValLongUri', by.y = 'Node_Long', sort = F, allow.cartesian = T, all.x = T) %>%
    .[is.na(Definition.x), Definition.x := cf.BaseLongUri] %>% .[is.na(Definition.y), Definition.y := cf.ValLongUri] %>%
    .[, .(rsc.ID, cf.BaseLongUri = Definition.x, cf.ValLongUri = Definition.y, distance, OntologyScope.x, OntologyScope.y)] %>%
    .[, .(distance = mean(distance)), .(rsc.ID, cf.BaseLongUri, cf.ValLongUri)] %>% setorder(distance, rsc.ID)
}

#' Precompute Tags
#' 
#' Get all the expanded ontology tags within a given taxon.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
precomputeTags <- function(taxa = getOption('app.taxa')) {
  graph <- simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(ChildNode_Long, ParentNode_Long)]))
  graphTerms <- unique(ONTOLOGIES[, as.character(ChildNode_Long, ParentNode_Long)])
  
  mTags <- DATA.HOLDER[[taxa]]@experiment.meta[, .(cf.BaseLongUri = Filter(function(tag) tag %in% graphTerms,
                                                                           unique(parseListEntry(as.character(cf.BaseLongUri))))), rsc.ID] %>%
    merge(DATA.HOLDER[[taxa]]@experiment.meta[, .(cf.ValLongUri = Filter(function(tag) tag %in% graphTerms,
                                                                         unique(parseListEntry(as.character(cf.ValLongUri))))), rsc.ID],
          by = 'rsc.ID', sort = F, allow.cartesian = T)
  
  mTags <- mTags[, rbindlist(lapply(.SD[, unique(cf.BaseLongUri)], function(uri) {
    cf.BaseLongUri <- igraph::subcomponent(graph, uri, 'out')
    distance <- igraph::distances(graph, uri, cf.BaseLongUri)
    
    data.table(tag = names(cf.BaseLongUri), type = 'cf.BaseLongUri', distance = c(distance))
  })) %>% rbind(
    rbindlist(lapply(.SD[, unique(cf.ValLongUri)], function(uri) {
      cf.ValLongUri <- igraph::subcomponent(graph, uri, 'out')
      distance <- igraph::distances(graph, uri, cf.ValLongUri)
      
      data.table(tag = names(cf.ValLongUri), type = 'cf.ValLongUri', distance = c(distance))
    }))
  ), rsc.ID] %>% .[is.finite(distance)] %>%
    rbind(
      DATA.HOLDER[[taxa]]@experiment.meta[, .(tag = Filter(function(tag) !(tag %in% graphTerms),
                                                                      unique(parseListEntry(as.character(cf.BaseLongUri))))), rsc.ID] %>%
        .[, .(tag, type = 'cf.BaseLongUri', distance = 0)] %>%
        rbind(
          DATA.HOLDER[[taxa]]@experiment.meta[, .(tag = Filter(function(tag) !(tag %in% graphTerms),
                                                               unique(parseListEntry(as.character(cf.ValLongUri))))), rsc.ID] %>%
            .[, .(tag, type = 'cf.ValLongUri', distance = 0)]
        )
    )
  
  mTags[, (expand.grid(cf.BaseLongUri = .SD[type == 'cf.BaseLongUri', tag],
                       cf.ValLongUri = .SD[type == 'cf.ValLongUri', tag])), rsc.ID] %>%
    merge(mTags[type == 'cf.BaseLongUri', .(rsc.ID, tag, distance)],
          by.x = c('rsc.ID', 'cf.BaseLongUri'), by.y = c('rsc.ID', 'tag'), sort = F, allow.cartesian = T) %>%
    merge(mTags[type == 'cf.ValLongUri', .(rsc.ID, tag, distance)],
          by.x = c('rsc.ID', 'cf.ValLongUri'), by.y = c('rsc.ID', 'tag'), sort = F, allow.cartesian = T) %>%
    .[, c('distance', 'distance.x', 'distance.y') := list(rowMeans2(cbind(distance.x, distance.y)), NULL, NULL)]
}

#' Enrich
#' 
#' Given rankings (@seealso search), generate a ranking-weighted count of all terms that can be
#' derived from tags present in the experiment (both cf.ValLongUri and ee.TagLongUri).
#'
#' @param rankings A named numeric (@seealso search).
#' @param taxa The taxon
#' @param scope The ontology scope.
#' @param options The options
#' @param session The Shiny session.
enrich <- function(rankings, taxa = getOption('app.taxa'), scope = getOption('app.ontology'),
                   options = getOption('app.all_options'), session = NULL) {
  rankings <- data.table(rsc.ID = rownames(rankings), rank = rankings$score, direction = rankings$direction)
  
  advanceProgress(session, 'Looking up ontology terms')
  
  # mMaps[[1]] is the background tags, mMaps[[2]] is the foreground tags.
  mMaps <- list(getTags(taxa, scope, NULL, session), getTags(taxa, scope, rankings[, rsc.ID], session))
  
  mMaps[[2]] <- mMaps[[2]] %>% merge(rankings, by = 'rsc.ID', sort = F, allow.cartesian = T)
  
  aprior <- function() {
    mMaps[[2]][distance <= options$distance,
               .(A = sum(rank) / (1 + sum(distance))), .(cf.BaseLongUri, cf.ValLongUri)][, B := sum(A)] %>%
      merge(mMaps[[1]][distance <= options$distance & !(rsc.ID %in% rankings[, rsc.ID]),
                       .(C = .N / (1 + sum(distance))), .(cf.BaseLongUri, cf.ValLongUri)][, D := sum(C)],
            by = c('cf.BaseLongUri', 'cf.ValLongUri'), sort = F, allow.cartesian = T)
  }
  
  chisq <- function(input) {
    advanceProgress(session, 'Running tests (1/2)')
    
    input[, c('pv.chisq', 'chisq') :=
            ifelse((A * (D - C)) / ((B - A) * C) <= 1, list(1, 1),
                   suppressWarnings(chisq.test(matrix(c(A, B - A, C, D - C), nrow = 2)))[c('p.value', 'statistic')]),
          .(cf.BaseLongUri, cf.ValLongUri)]
  }
  
  fisher <- function(input) {
    advanceProgress(session, 'Running tests (2/2)')
    
    input[, c('pv.fisher', 'OR') :=
            suppressWarnings(fisher.test(matrix(c(A, B - A, C, D - C), nrow = 2),
                                         alternative = 'greater', conf.int = F))[c('p.value', 'estimate')],
          .(cf.BaseLongUri, cf.ValLongUri)]
  }

  aprior() %>% chisq %>% fisher %>% setorder(pv.chisq)
}
