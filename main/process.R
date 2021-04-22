#' Weighted correlation
#'
#' @param x A vector or matrix (whose columns will be correlated)
#' @param y A vector with compatible dimensions of x to correlate against
#' @param w A vector of weights, the same length as y
#'
#' @return The weighted Pearson correlation of x (or the columns of x) and y
cor.wt <- function(x, y, w = rep(1, length(y))) {
  if(is.matrix(x))
    apply(x, 2, cor.wt, y = y, w = w)
  else {
    sw <- sum(w)
    
    mx <- sum(w * x) / sw
    my <- sum(w * y) / sw
    
    (sum(w * (x - mx) * (y - my)) / sw) / sqrt((sum(w * (x - mx)^2) / sw) * (sum(w * (y - my)^2) / sw))
  }
}

#' Search
#' 
#' Uses the M-VSM to sort experiments that show at least one of the genes as DE.
#'
#' @param genes A list of Entrez Gene IDs (ie. 1, 22, 480) as characters
#' @param options Optional extra parameters to pass, such as:
#' * pv: An FDR cutoff (default: 0.05)
#' * fc.lower / fc.upper: Upper and lower logFC thresholds (default: 0 / 10)
#' * mfx: Whether or not to scale by gene multifunctionality
#' * geeq: Whether or not to scale by GEEQ score
#' @param inprod Whether to use the generated null distribution or not
search <- function(genes, options = getConfig(), inprod = T, DATA = NULL) {
  if(is.null(DATA)) # TODO Remove after doing scores as this might trigger a full copy
    DATA <- DATA.HOLDER
  
  mData <- DATA[[options$taxa$value]]
  
  # Only retain GOI
  geneMask <- which(mData@gene.meta$entrez.ID %in% genes)
  
  n.genes <- length(geneMask)
  if(n.genes == 0)
    return(NULL)
  
  # No searching if we have this shortcut
  # A more elegant solution would be nice
  if(inprod && options$method$value == 'diff')
    return(geneMask)
  
  query <- suppressWarnings(options$sig$value %>% as.numeric)
  # TODO Should signal why it stopped (mismatch of signature length)
  if(length(query) > 1 && length(query) != n.genes)
    return(NULL)
  
  # P-values for only the GOI
  pv <- mData@data$adj.pv[geneMask, ]
  
  pv[is.na(pv)] <- 1
  
  if(n.genes == 1)
    pv <- t(pv)
  
  # Only retain experiments that have at least one of the GOI as DE (pv < threshold)
  experimentMask <- pv <= options$pv$value
  logFC <- mData@data$fc[geneMask, ]
  zScore <- mData@data$zscore[geneMask, ]
  
  logFC[is.na(logFC)] <- 0
  zScore[is.na(zScore)] <- 0
  
  if(n.genes == 1) {
    logFC <- t(logFC)
    zScore <- t(zScore)
  }
  
  logFC <- logFC %>% as.data.table
  zScore <- zScore %>% as.data.table
  pv <- pv %>% as.data.table
  
  # Mask experiments that have at least n of the GOI with threshold_u > abs(logFC) > threshold_l
  if(options$fc$value[1] != 0 || options$fc$value[2] != 0) {
    exMask <- abs(logFC) >= options$fc$value[1]
    
    if(options$fc$value[1] != options$fc$value[2])
      exMask <- exMask & (abs(logFC) <= options$fc$value[2])
    
    experimentMask <- experimentMask & exMask
  }
  
  # Fail if 0 experiments show these genes as DE
  # TODO Not technically necessary anymore?
  if(sum(experimentMask) == 0)
    return(NULL)
  
  # Number of DEGs for experiments that pass thresholds
  experimentMeta <- mData@experiment.meta %>% copy %>% as.data.frame %>%
    `rownames<-`(.[, 'rsc.ID']) %>% .[, c('ee.qScore', 'n.DE', 'n.detect')]
  
  experimentMeta$n.DE[is.na(experimentMeta$n.DE)] <- 0
  
  if(!options$geeq$value)
    experimentMeta$ee.qScore <- 2
  else
    experimentMeta$ee.qScore <- pmax(experimentMeta$ee.qScore + 1, 0, na.rm = T)
  experimentMeta$ee.qScore <- sqrt(experimentMeta$ee.qScore / 2)
  
  if(length(query) == 0 || any(is.na(query))) {
    # TODO Consider warning if length(query) > 1 (some number was invalid)
    zScore <- abs(zScore)
    
    # TODO maybe this should come after pv weight
    if(n.genes == 1)
      query <- max(zScore %>% as.matrix %>% `*`(experimentMask))
    else
      query <- Rfast::rowMaxs(zScore %>% as.matrix %>% `*`(experimentMask), value = T)
  } else
    zScore <- logFC # Doesn't make sense to query a z-score
  # TODO it might if we only care about directionality
  
  if(options$mfx$value)
    MFX_WEIGHT <- 1 - mData@gene.meta$mfx.Rank[geneMask]
  else
    MFX_WEIGHT <- rep(1, n.genes)
  
  zScore <- zScore * -log10(Rfast::Pmax(matrix(1e-10, ncol = ncol(pv), nrow = nrow(pv)), as.matrix(pv)))
  
  mNames <- c(mData@gene.meta$gene.Name[geneMask], 'score')
  
  if(options$method$value == 'diff') {
    ret <- (zScore - query) %>% `*`(MFX_WEIGHT) %>% t %>% as.data.table
    mNames <- mData@gene.meta$gene.Name[geneMask]
  } else if(options$method$value == 'cor') { # TODO Other methods should also mfx weight the gene matrix
    ret <- zScore %>% t %>% as.data.table %>%
      .[, score := abs(cor.wt(zScore %>% as.matrix, query, MFX_WEIGHT))] %>% # TODO abs?
      .[is.nan(score), score := 0]
  } else if(options$method$value == 'mvsm') {
    idf <- Rfast::Log(1 / MFX_WEIGHT) + 1
    zScore[, query := query]
    zScore <- zScore * idf
    
    # Extract the query
    query <- zScore[, query]
    zScore[, query := NULL]
    
    # Cosine similarity
    zScore <- as.matrix(zScore)
    
    cross_x <- crossprod(query)
    scores <- sapply(1:ncol(zScore), function(i) {
      cross_q <- crossprod(query, zScore[, i])
      if(cross_q == 0) 0
      else cross_q / sqrt(cross_x * crossprod(zScore[, i]))
    })
    
    ret <- zScore %>% t %>% as.data.table %>%
      .[, score := abs(scores)] %>%
      .[is.nan(score), score := 0]
  }
  
  # TODO it may be more appropriate to add and divide by constants
  # (ie. 1e5, sufficiently large enough to guarantee it's always positive)
  # to prevent imprudent score deviations.
  # This seems to actually make it worse
  experimentN <- Rfast::colsums(experimentMask)
  
  ret %>% setnames(mNames) %>%
    .[, is.passing := experimentN >= options$req$value] %>%
    .[, f.IN := experimentN / n.genes] %>%
    .[, f.OUT := pmax(0, experimentMeta$n.DE - experimentN) / experimentMeta$n.detect] %>%
    .[, ee.q := experimentMeta$ee.qScore] %>%
    .[, rn := colnames(zScore)]
  
  if('score' %in% colnames(ret)) {
    ret[, score := (score + abs(min(score))) * ee.q * (1 + f.IN) / (1 + 10^(f.OUT))] %>%
      .[, score := score / max(score)] %>%
      setorder(-score)
  } else {
    ret <- data.table(
      rn = ret[, rn],
      ret[, lapply(1:n.genes, function(i) {
        col <- .SD[, i, with = F]
        ((col + abs(min(col))) * ee.q * (1 + experimentMask[i, ] / n.genes) /
            (1 + 10^(pmax(0, experimentMeta$n.DE - experimentMask[i, ]) / experimentMeta$n.detect))) %>%
          `/`(max(.))
      }), .SDcols = !c('rn', 'is.passing', 'f.IN', 'f.OUT', 'ee.q')],
      ret[, .(is.passing, f.IN, f.OUT, ee.q)]
    )
  }
  
  # TODO What happens when these stats are negative? Taking abs for now
  if(exists('NULLS.EXP', envir = globalenv()) && options$method$value == 'diff') {
    ret %>%
      merge(NULLS.EXP[[options$taxa$value]], by = 'rn', all = T, sort = F) %>% {
        data.table(
          rn = .[, rn],
          .[, lapply(.SD, function(col) {
            abs((col - M1) / S1)
          }), .SDcols = !c('rn', 'is.passing', 'f.IN', 'f.OUT', 'ee.q', 'M1', 'S1')],
          .[, .(is.passing, f.IN, f.OUT, ee.q)]
        )
      }
  } else
    ret
}

#' getTags
#' 
#' Get tags for the specified experiments within a specified distance.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat, any].
#' @param rsc.IDs A list of experiment rsc IDs or NULL for everything
#' @param max.distance The maximum tree traversal distance to include
#' @param inv Whether or not to inverse the selected rscs
#' @param CACHE A cache to use. If null, uses the global CACHE.BACKGROUND
getTags <- function(taxa = getConfig(key = 'taxa')$value,
                    rsc.IDs = NULL, max.distance = Inf, inv = F, CACHE = NULL) {
  if(length(taxa) > 1)
    return(rbindlist(lapply(taxa, getTags, rsc.IDs, max.distance, inv, CACHE)))
  
  if(is.null(CACHE)) # TODO Remove after doing scores as this might trigger a full copy
    CACHE <- CACHE.BACKGROUND
  
  if(is.null(rsc.IDs))
    rsc.IDs <- CACHE[[taxa]][, unique(as.character(rsc.ID))]
  
  if(inv)
    rsc.IDs <- CACHE[[taxa]][!(rsc.ID %in% rsc.IDs), unique(as.character(rsc.ID))]
  
  CACHE[[taxa]][distance <= max.distance & rsc.ID %in% rsc.IDs]
}

#' Precompute Tags
#' 
#' Get all the expanded ontology tags within a given taxon.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat, any].
#' @param mGraph The igraph form of the ontologies. Computes internally if not supplied
#' @param graphTerms The unique ontology terms. Computes internally if not supplied
#' @param DATA Data to use. If null, uses the global DATA.HOLDER
#' @param POST Whether to do "post-pre-processing" like reordering and trimming
precomputeTags <- function(taxa = getConfig(key = 'taxa')$value, mGraph = NULL, graphTerms = NULL,
                           DATA = NULL, POST = T) {
  if(is.null(mGraph))
    mGraph <- simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(as.character(ChildNode_Long), as.character(ParentNode_Long))]))
  if(is.null(graphTerms))
    graphTerms <- unique(ONTOLOGIES[, as.character(ChildNode_Long, ParentNode_Long)])
  
  # Tags that are simple ontology terms
  mSimpleTags <- DATA[[taxa]]@experiment.meta[, .(tag = unique(as.character(cf.BaseLongUri)), type = 'cf.BaseLongUri'), .(rsc.ID, ee.ID)] %>%
    .[tag %in% graphTerms] %>%
    rbind(DATA[[taxa]]@experiment.meta[, .(tag = unique(as.character(cf.ValLongUri)), type = 'cf.ValLongUri'), .(rsc.ID, ee.ID)] %>%
            .[tag %in% graphTerms])
  
  # Tags that are "bagged" (ie. have multiple tags)
  bagged <- DATA[[taxa]]@experiment.meta[, grepl('; ', as.character(cf.BaseLongUri), fixed = T) |
                                           grepl('; ', as.character(cf.ValLongUri), fixed = T)]
  
  # Structured text entries (ie. non-ontology terms)
  mStructuredTags <- DATA[[taxa]]@experiment.meta[!bagged, .(tag = unique(as.character(cf.BaseLongUri)),
                                                             type = 'cf.BaseLongUri'), .(rsc.ID, ee.ID)] %>%
    .[!(tag %in% graphTerms)] %>%
    rbind(DATA[[taxa]]@experiment.meta[!bagged, .(tag = unique(as.character(cf.ValLongUri)),
                                                  type = 'cf.ValLongUri'), .(rsc.ID, ee.ID)] %>%
            .[!(tag %in% graphTerms)]) %>% .[, distance := 0]
  
  # Expand the bagged tags so there's one row per entry
  if(sum(bagged) == 0)
    mBagOfWords <- data.table()
  else
    mBagOfWords <- DATA[[taxa]]@experiment.meta[bagged, .(tag = unique(as.character(cf.BaseLongUri)),
                                                          type = 'cf.BaseLongUri'), .(rsc.ID, ee.ID)] %>%
    .[!(tag %in% graphTerms)] %>%
    rbind(DATA[[taxa]]@experiment.meta[bagged, .(tag = unique(as.character(cf.ValLongUri)),
                                                 type = 'cf.ValLongUri'), .(rsc.ID, ee.ID)] %>%
            .[!(tag %in% graphTerms)]) %>%
    .[, lapply(.SD, function(x) parseListEntry(as.character(x))), .(rsc.ID, ee.ID, type)] %>%
    .[, ID := 1:length(tag), .(rsc.ID, ee.ID, type)]
  
  if(nrow(mBagOfWords) > 0)
    mComputable <- union(mSimpleTags[, unique(tag)], mBagOfWords[tag %in% graphTerms, tag])
  else
    mComputable <- mSimpleTags[, unique(tag)]
  
  # Compute ontology expansions on ontology terms
  # Applies to both simple tags and bag of word tags that expanded into ontology terms
  mComputedTags <- rbindlist(lapply(mComputable, function(uri) {
    tag <- igraph::subcomponent(mGraph, uri, 'out')
    distance <- igraph::distances(mGraph, uri, tag)
    data.table(startTag = uri, tag = names(tag), distance = c(distance))
  }))
  
  if(nrow(mComputedTags) > 0)
    # Simple (ontology) tags get associated with their parents from mComputedTags
    mTags <- mSimpleTags %>% merge(mComputedTags[startTag %in% mSimpleTags[, unique(tag)]],
                                   by.x = 'tag', by.y = 'startTag', sort = F, allow.cartesian = T) %>%
    .[, c('tag', 'tag.y') := list(tag.y, NULL)] %>% .[, ID := NA]
  else
    mTags <- data.table()
  
  # Structured tags just get inserted
  mTags <- mTags %>% rbind(mStructuredTags[, ID := NA] %>% .[, .(tag, rsc.ID, ee.ID, type, distance, ID)])
  
  if(nrow(mBagOfWords) > 0 && nrow(mComputedTags) > 0)
    # Ontology-expanded bag of word tags get associated with their parents from mComputedTags
    mTags <- mTags %>% rbind(
      mBagOfWords %>%
        merge(mComputedTags[startTag %in% mBagOfWords[, unique(tag)]],
              by.x = 'tag', by.y = 'startTag', all = T, sort = F, allow.cartesian = T) %>%
        .[is.na(tag.y), c('tag.y', 'distance') := list(tag, 0)] %>%
        .[, c('tag', 'tag.y') := list(tag.y, NULL)]
    )
  
  mTags <- mTags %>%
    merge(unique(ONTOLOGIES.DEFS[, .(Node_Long = as.character(Node_Long), Definition = as.character(Definition))]),
          by.x = 'tag', by.y = 'Node_Long', sort = F, allow.cartesian = T, all.x = T) %>%
    .[is.na(Definition), Definition := tag] %>%
    .[, .(rsc.ID, ee.ID, type, tag = Definition, distance, ID)]
  
  # Do a grid expansion of bagged terms, maintaining indexed ordering
  if(nrow(mBagOfWords) > 0 && nrow(mTags[!is.na(ID) & distance < 1]) > 0)
    mTagsExpanded <- mTags # TODO this loses tags?
  #  mTagsExpanded <- mTags %>% .[!is.na(ID) & distance < 1, expand.grid(aggregate(.SD[, tag], by = list(.SD[, ID]), FUN = list)[[-1]]) %>%
  #                                 apply(1, paste0, collapse = '; ') %>% unique, .(rsc.ID, ee.ID, type)] %>%
  #  .[, c('tag', 'V1') := list(V1, NULL)] %>% rbind(mTags[is.na(ID) | distance >= 1, .(rsc.ID, ee.ID, type, tag)])
  else
    mTagsExpanded <- mTags
  
  mTagsExpanded %>% .[, .(rsc.ID, ee.ID, type, tag)] %>%
    
    # Now expand these expanded terms
    .[, expand.grid(cf.BaseLongUri = .SD[type == 'cf.BaseLongUri', tag],
                    cf.ValLongUri = .SD[type == 'cf.ValLongUri', tag]), .(rsc.ID, ee.ID)] %>% unique %>%
    
    # Add back distances
    merge(mTags[type == 'cf.BaseLongUri', .(rsc.ID, ee.ID, tag, distance)],
          by.x = c('rsc.ID', 'ee.ID', 'cf.BaseLongUri'), by.y = c('rsc.ID', 'ee.ID', 'tag'), all.x = T, sort = F, allow.cartesian = T) %>%
    merge(mTags[type == 'cf.ValLongUri', .(rsc.ID, ee.ID, tag, distance)],
          by.x = c('rsc.ID', 'ee.ID', 'cf.ValLongUri'), by.y = c('rsc.ID', 'ee.ID', 'tag'), all.x = T, sort = F, allow.cartesian = T) %>%
    
    # Recombinants have unknown distances. Let's just make them 0
    .[is.na(distance.x), distance.x := 0] %>% .[is.na(distance.y), distance.y := 0] %>%
    
    # Compute net distance as the mean of distances to the base and contrasting factor
    .[, c('distance', 'distance.x', 'distance.y') := list(Rfast::rowmeans(cbind(distance.x, distance.y)), NULL, NULL)] %>%
    
    # And add back categories
    merge(DATA[[taxa]]@experiment.meta[, .(rsc.ID, ee.ID, cf.Cat)], by = c('rsc.ID', 'ee.ID'), all.x = T, sort = F, allow.cartesian = T) %>%
    .[cf.BaseLongUri != cf.ValLongUri] %>% {
      if(!POST) {
        .[, .(rsc.ID, ee.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri, reverse = F, distance)]
      } else {
        # Finally reorder
        reorderTags(.) %>%
          
          # And use some heuristics to make our corpus a little more lean
          .[, N := .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
          .[distance == 0 | (N > 1 & N < 500 & distance < 5)] %>%
          
          # This is annoying but the filter should be applied twice
          .[, N := .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
          .[distance == 0 | (N > 1 & N < 500 & distance < 5)] %>%
          
          .[, !'N']
      }
    }
}

#' Reorder Tags
#' 
#' Reorders the baseline and contrast so that the most common one is in the baseline.
#'
#' @param cache The cache entry
#'
#' @return A reordered cache entry
reorderTags <- function(cache) {
  vals <- cache[, .N, cf.BaseLongUri] %>%
    merge(cache[, .N, cf.ValLongUri],
          by.x = 'cf.BaseLongUri', by.y = 'cf.ValLongUri', all = T, sort = F, allow.cartesian = T) %>%
    .[is.na(N.x), N.x := 0] %>% .[is.na(N.y), N.y := 0] %>%
    .[, N := N.x + N.y, cf.BaseLongUri] %>%
    .[, .(tag = cf.BaseLongUri, N)]
  
  cache[, .(rsc.ID, ee.ID, cf.Cat, b = cf.BaseLongUri, v = cf.ValLongUri, distance)] %>%
    merge(vals, by.x = 'b', by.y = 'tag', sort = F, allow.cartesian = T) %>%
    merge(vals, by.x = 'v', by.y = 'tag', sort = F, allow.cartesian = T) %>%
    .[, reverse := (N.y > N.x) | (N.y == N.x & as.character(b) < as.character(v))] %>%
    .[reverse == T, c('b', 'v') := list(v, b)] %>%
    .[, .(rsc.ID = as.factor(rsc.ID), ee.ID,
          cf.Cat = as.factor(cf.Cat), cf.BaseLongUri = as.factor(b), cf.ValLongUri = as.factor(v),
          reverse, distance)]
}

#' Reorder Tags
#' 
#' Reorders the baseline and contrast so that the most common one is in the baseline for the expanded
#' set of tags.
#'
#' @param cache The cache entry
#'
#' @return A reordered cache entry
reorderTags2 <- function(cache) {
  vals <- cache[, .N, cf.BaseLongUri] %>%
    merge(cache[, .N, cf.ValLongUri],
          by.x = 'cf.BaseLongUri', by.y = 'cf.ValLongUri', all = T, sort = F, allow.cartesian = T) %>%
    .[is.na(N.x), N.x := 0] %>% .[is.na(N.y), N.y := 0] %>%
    .[, N := N.x + N.y, cf.BaseLongUri] %>%
    .[, .(tag = cf.BaseLongUri, N)]
  
  cache[, .(rsc.ID, ee.ID, cf.Cat, b = cf.BaseLongUri, v = cf.ValLongUri, distance, oreverse = reverse)] %>%
    merge(vals, by.x = 'b', by.y = 'tag', sort = F, allow.cartesian = T) %>%
    merge(vals, by.x = 'v', by.y = 'tag', sort = F, allow.cartesian = T) %>%
    .[, reverse := (N.y > N.x) | (N.y == N.x & as.character(b) < as.character(v))] %>%
    .[reverse == T, c('b', 'v') := list(v, b)] %>%
    .[, reverse := ifelse(oreverse, !reverse, reverse)] %>%
    .[, .(rsc.ID = as.factor(rsc.ID), ee.ID,
          cf.Cat = as.factor(cf.Cat), cf.BaseLongUri = as.factor(b), cf.ValLongUri = as.factor(v),
          reverse, distance)]
}

# TODO this could be merged with reorderTags2
reorderTags3 <- function(data) {
  vals <- data[, .N, cf.BaseLongUri] %>%
    merge(data[, .N, cf.ValLongUri],
          by.x = 'cf.BaseLongUri', by.y = 'cf.ValLongUri', all = T, sort = F, allow.cartesian = T) %>%
    .[is.na(N.x), N.x := 0] %>% .[is.na(N.y), N.y := 0] %>%
    .[, N := N.x + N.y, cf.BaseLongUri] %>%
    .[, .(tag = cf.BaseLongUri, N)]
  
  data %>% copy %>%
    setnames(c('cf.BaseLongUri', 'cf.ValLongUri'), c('b', 'v')) %>%
    merge(vals, by.x = 'b', by.y = 'tag', sort = F, allow.cartesian = T) %>%
    merge(vals, by.x = 'v', by.y = 'tag', sort = F, allow.cartesian = T) %>%
    .[, reverse := (N.y > N.x) | (N.y == N.x & as.character(b) < as.character(v))] %>%
    .[reverse == T, c('b', 'v') := list(v, b)] %>%
    .[, c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri') :=
        list(as.factor(cf.Cat), as.factor(b), as.factor(v))] %>%
    .[, !c('reverse', 'b', 'v', 'N.x', 'N.y')]
}

# TODO this could be memoized when inprod && is.integer(rankings)

#' Enrich
#' 
#' Given rankings (@seealso search), generate a ranking-weighted count of all terms that can be
#' derived from tags present in the experiment.
#'
#' @param rankings A named numeric (@seealso search).
#' @param options The options
#' @param inprod Whether to use the generated null distribution or not
#' @param keepopen How many of the entire files (used when @param inprod is true) to keep in memory (per species)
#' @param CACHE A cache to use. If null, uses the global CACHE.BACKGROUND
enrich <- function(rankings, options = getConfig(), inprod = T, keepopen = 0, CACHE = NULL) {
  if(inprod && is.numeric(rankings))
    return(enrichMem(rankings, options, keepopen))
  
  mMaps <- getTags(options$taxa$value, rankings$rn, CACHE = CACHE)
  
  # Stat of 0 means it's as expected (0 SD from mean)
  mMaps <- mMaps[as.character(cf.BaseLongUri) != as.character(cf.ValLongUri)] %>%
    merge(rankings[, !c('is.passing', 'f.IN', 'f.OUT', 'ee.q')],
          by.x = 'rsc.ID', by.y = 'rn', sort = F, allow.cartesian = T)
  
  for(i in 8:ncol(mMaps))
    set(mMaps, which(is.na(mMaps[[i]])), i, 0)
  
  mMaps[, lapply(.SD, function(stat) { sum(stat + 1, na.rm = T) / .N }),
        .(cf.Cat, cf.BaseLongUri, cf.ValLongUri), .SDcols = !c('reverse', 'distance', 'rsc.ID', 'ee.ID')] %>%
    merge(mMaps[, .(distance = round(mean(distance, na.rm = T), 2)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
          by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri')) %>%
    .[, stat := Rfast::rowmeans(as.matrix(.SD[, !c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri', 'distance')]))] %>%
    setorder(-stat, distance) %>%
    .[, .SD[1], stat]
}

maybeMemoise <- function(FUN, fname) {
  if(!exists(fname, envir = globalenv())) {
    if(Sys.getenv('RSTUDIO') == '1' || getOption('force.memoise', F))
      FUN <- memoise::memoise(FUN)
    
    assign(fname, FUN, envir = globalenv())
  }
}

enrichPreMem <- function(rankings, options, keepopen = 0, forceCache = NULL) {
  enrichMem(rankings, options, keepopen)
}

maybeMemoise(function(rankings, options, keepopen = 0) {
  if(exists('forceCache', envir = parent.frame(2)) && !is.null(get('forceCache', envir = parent.frame(2))))
    return(get('forceCache', envir = parent.frame(2)))
  
  if(keepopen > 0 && !exists('SINGLES'))
    assign('SINGLES', new.env(parent = globalenv()), envir = globalenv())
  
  # TODO it would be nice to have smaller clusters or cluster commonly searched genes together
  lapply(unique(ceiling(rankings / getOption('chunk.size', 200))), function(filenum) {
    mFile <- paste0(options$taxa$value, '_', filenum)
    if(exists('SINGLES', envir = globalenv()) &&
       exists(mFile, envir = SINGLES))
      fileConnection <- get(mFile, envir = SINGLES)
    else
      fileConnection <- readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/singlegene/', mFile, '.rds'))
    
    fileConnection %>% {
      if(keepopen > 0 && !exists(mFile, envir = SINGLES)) {
        rm(list = ls(pattern = paste0(options$taxa$value, '_.*'), envir = SINGLES) %>%
             head(ifelse(length(.) == keepopen, 1, 0)), envir = SINGLES)
        assign(mFile, ., envir = SINGLES)
      }
      .
    } %>% .[I %in% rankings] %>% dcast(... ~ I, value.var = 'stat')
  }) %>% Reduce(function(...) merge(..., by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri', 'distance')), .) %>%
    .[, stat := Rfast::rowmeans(as.matrix(.SD[, !c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri', 'distance')]))] %>%
    setorder(-stat, distance) %>%
    .[, .SD[1], stat]
}, 'enrichMem')
