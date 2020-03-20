library(igraph)
library(data.table)
library(dplyr)
library(matrixStats)

#' Search
#' 
#' Uses the M-VSM to sort experiments that show at least one of the genes
#' as DE.
#'
#' @param genes A list of Entrez Gene IDs (ie. 1, 22, 480) as characters
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
#' @param options Optional extra parameters to pass, such as:
#' * pv: A p-value cutoff (default: 0.05)
#' * fc: A logFC threshold (default: 0)
#' * fc.lower / fc.upper: Upper and lower logFC thresholds (default: 0 / 0)
#'
#' @return A named numeric with experiment IDs as names and scores as entries, or NULL if no GOI
#' or EOI are found in the underlying data structure for @param taxa and @param options$pv / @param options$fc.
search <- function(genes, taxa = 'human', options = list()) {
  P_THRESHOLD <- ifelse(is.null(options$pv), 0.05, options$pv)
  if(is.null(options$fc)) {
    FC_L_THRESHOLD <- ifelse(is.null(options$fc.lower), 0, options$fc.lower)
    FC_U_THRESHOLD <- options$fc.upper
  } else {
    FC_L_THRESHOLD <- ifelse(is.null(options$fc), 0, options$fc)
    FC_U_THRESHOLD <- NULL
  }
  
  # Data Extraction ---------------------------------------------------------

  # Only retain GOI
  rowFilter <- which(rownames(DATA.HOLDER[[taxa]]@data$fc) %in% genes)
  n.genes <- length(rowFilter)
  if(n.genes == 0) return(NULL)
  
  # P-values for only the GOI
  pv <- DATA.HOLDER[[taxa]]@data$adj.pv[rowFilter, ]
  
  # Only retain experiments that have at least one of the GOI as DE (pv < threshold)
  colFilter <- colSums2(pv < P_THRESHOLD, na.rm = T) > 0
  pv <- pv[, colFilter] %>% as.data.table
  
  # logFCs for only the GOI/EOI and maintain structure.
  logFC <- DATA.HOLDER[[taxa]]@data$fc[rowFilter, colFilter] %>% as.data.table %>% abs
  
  # Only retain experiments that have at least one of the GOI with abs(logFC) > threshold
  if(FC_L_THRESHOLD != 0 || !is.null(FC_U_THRESHOLD)) {
    exFilter <- logFC > FC_L_THRESHOLD
    if(!is.null(FC_U_THRESHOLD))
      exFilter <- exFilter & logFC < FC_U_THRESHOLD
    
    colFilter <- colSums2(exFilter, na.rm = T) > 0
    logFC <- logFC[, colFilter, with = F]
    pv <- pv[, colFilter, with = F]
  }
  
  # Data Processing ---------------------------------------------------------
  
  tf <- logFC# * -log10(pv)
  query <- rowMaxs(tf %>% as.matrix, na.rm = T)
  
  idf <- log10(ncol(tf) / (1 + rowSums2(tf != 0, na.rm = T))) + 1
  tf[, query := query]
  tfidf <- tf * idf
  
  if(nrow(tfidf) > 1)
    tfidf <- scale(tfidf, center = F, scale = sqrt(colSums2(as.matrix(tfidf^2), na.rm = T))) %>% as.data.table
  
  query <- tfidf[, query]
  tfidf <- tfidf[, query := NULL]
  
  tfidf <- as.matrix(tfidf)
  tfidf[is.nan(tfidf)] <- 0
  scores <- query %*% tfidf
  
  names(scores) <- colnames(tfidf)
  scores[order(-scores)]
}

#' Condition Count
#' 
#' Given rankings (@seealso search), generate a ranking-weighted count of all terms that can be
#' derived from tags present in the experiment (both cf.ValLongUri and ee.TagLongUri).
#'
#' @param rankings A named numeric (@seealso search).
#' @param taxa 
#' @param scope 
#'
#' @return
conditionCount <- function(rankings, taxa = 'human', scope = ONTOLOGIES[, unique(OntologyScope)]) {
  rankings <- data.table(rsc.ID = names(rankings), rank = rankings)
  
  ontoScope <- ONTOLOGIES[OntologyScope %in% scope, .(ChildNode_Long, ParentNode_Long)]
  graph <- igraph::graph_from_data_frame(ontoScope)
  
  tags <- DATA.HOLDER[[taxa]]@experiment.meta[rsc.ID %in% rankings[, rsc.ID] &
                                                !is.na(ee.TagUri) | !is.na(cf.ValUri)] %>%
    .[is.na(ee.TagLongUri) | ee.TagLongUri %in% ontoScope[, ChildNode_Long]] %>%
    .[is.na(cf.ValLongUri) | cf.ValLongUri %in% ontoScope[, ChildNode_Long]]
  
  do.call(rbind, lapply(1:nrow(tags), function(indx) {
    do.call(rbind, lapply(Filter(Negate(is.na), c(parseListEntry(tags[indx, ee.TagLongUri]),
                                                  parseListEntry(tags[indx, cf.ValLongUri]))), function(entry) {
      data.table(rsc.ID = tags[indx, rsc.ID], tag = names(igraph::subcomponent(graph, entry, 'out')))
    }))
  }))[tag != 'NA'] -> mMap
  
  # Merge in definitions and ranks...
  mMap <- unique(ONTOLOGIES.DEFS[OntologyScope %in% scope & Node_Long %in% mMap[, unique(tag)],
                                 Definition, Node_Long]) %>%
    merge(mMap, by.x = 'Node_Long', by.y = 'tag') %>% unique %>% merge(rankings, by = 'rsc.ID')
  
  mMap[, .(.N, `ranked N` = .N * sum(rank)), Definition] %>% .[!(Definition %in% BLACKLIST)] %>% setorder(-N)
}

#' Parse List Entry
#' 
#' Given an entry (that may or may not be NA) of semicolon delimited (; ) list (in which, any may be NA),
#' return a vector of those values, coercing NA strings into true NAs.
#'
#' @param entry NA or a character of semicolon delimited (; ) values.
#'
#' @return A character vector of values in the input character
parseListEntry <- function(entry) {
  if(is.na(entry) || !grepl('; ', entry, fixed = T)) return(entry)
  
  cleaned <- strsplit(entry, '; ', fixed = T) %>% unlist
  ifelse('NA' == cleaned, NA, cleaned)
}

# DATA.HOLDER$human@gene.meta[gene.Name %in% c('XIST', 'RPS4Y1', 'KDM5D'), entrez.ID] -> genes
