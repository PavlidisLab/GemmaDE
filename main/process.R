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
#' @param taxa A taxa scope. Can be one of [human, mouse, rat, any].
#' @param signature A DE signature to search or null to compute an "ideal" one
#' @param options Optional extra parameters to pass, such as:
#' * pv: An FDR cutoff (default: 0.05)
#' * fc.lower / fc.upper: Upper and lower logFC thresholds (default: 0 / 10)
#' * mfx: Whether or not to scale by gene multifunctionality
#' * geeq: Whether or not to scale by GEEQ score
search <- (function(genes, options = getConfig()) {
  if(options$taxa$value == 'any') {
    mData <- new('EData')
    
    for(x in options$taxa$core) {
      ID <- paste0(unname(options$taxa$mapping[x]), '_ID')
      IDs <- genes[, ..ID] %>% unlist %>% unname
      
      gMeta <- DATA.HOLDER[[x]]@gene.meta %>% copy %>% .[, I := .I]
      gMap <- merge(genes[, c('key', ID), with = F],
                    gMeta %>% .[entrez.ID %in% IDs, .(entrez.ID, I)],
                    by.x = ID, by.y = 'entrez.ID', all = T) %>%
        .[, c(ID) := list(NULL)]
      
      pvData <- DATA.HOLDER[[x]]@data$adj.pv[gMap$I, ] %>% `rownames<-`(gMap$key)
      pvData[is.na(pvData)] <- 1
      
      fcData <- DATA.HOLDER[[x]]@data$fc[gMap$I, ] %>% `rownames<-`(gMap$key)
      fcData[is.na(fcData)] <- 0
      
      zScoreData <- DATA.HOLDER[[x]]@data$zscore[gMap$I, ] %>% `rownames<-`(gMap$key)
      zScoreData[is.na(zScoreData)] <- 0
      
      mData@experiment.meta <- rbind(mData@experiment.meta, DATA.HOLDER[[x]]@experiment.meta)
      
      if(is.null(mData@data$adj.pv)) {
        mData@data$adj.pv <- pvData
        mData@data$fc <- fcData
        mData@data$zscore <- zScoreData

        mData@gene.meta <- gMeta[gMap$I, ] %>% .[, entrez.ID := gMap$key]
      } else {
        mData@data$adj.pv <- merge(mData@data$adj.pv, pvData, by = 0, all = T) %>% `rownames<-`(.[, 1]) %>% .[, -1]
        mData@data$fc <- merge(mData@data$fc, fcData, by = 0, all = T) %>% `rownames<-`(.[, 1]) %>% .[, -1]
        mData@data$zscore <- merge(mData@data$zscore, zScoreData, by = 0, all = T) %>% `rownames<-`(.[, 1]) %>% .[, -1]

        mData@gene.meta$gene.Name <- coalesce(mData@gene.meta$gene.Name, gMeta[gMap$I, gene.Name])
        mData@gene.meta$mfx.Rank <- rowMeans2(cbind(mData@gene.meta$mfx.Rank, gMeta[gMap$I, mfx.Rank]), na.rm = T)
      }
    }
    
    genes <- genes$key
    geneMask <- which(!is.na(mData@gene.meta$mfx.Rank))
  } else {
    mData <- DATA.HOLDER[[options$taxa$value]]

    # Only retain GOI
    geneMask <- which(mData@gene.meta$entrez.ID %in% genes)
  }
  
  n.genes <- length(geneMask)
  if(n.genes == 0)
    return(NULL)
  
  query <- suppressWarnings(options$sig$value %>% as.numeric)
  # TODO Should signal why it stopped (mismatch of signature length)
  if(length(query) > 1 && length(query) != n.genes)
    return(NULL)
  
  # P-values for only the GOI
  pv <- mData@data$adj.pv[geneMask, ]
  
  if(options$taxa$value != 'any') # For taxa == 'any', we do this on the fly
    pv[is.na(pv)] <- 1
  
  if(n.genes == 1)
    pv <- t(pv)
  
  # Only retain experiments that have at least one of the GOI as DE (pv < threshold)
  experimentN <- Rfast::colsums(pv <= options$pv$value) %>% as.integer
  experimentMask <- Rfast::colsums(pv <= options$pv$value) >= min(options$req$value, n.genes) %>% as.bit
  logFC <- mData@data$fc[geneMask, ]
  zScore <- mData@data$zscore[geneMask, ]
  
  if(options$taxa$value != 'any') {
    logFC[is.na(logFC)] <- 0
    zScore[is.na(zScore)] <- 0
  }
  
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
    
    experimentMask <- experimentMask & as.bit(Rfast::colsums(exMask) >= min(options$req$value, n.genes))
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
      query <- max(zScore %>% as.matrix %>% .[, as.which(experimentMask)])
    else
      query <- Rfast::rowMaxs(zScore %>% as.matrix %>% .[, as.which(experimentMask)], value = T)
  } else
    zScore <- logFC # Doesn't make sense to query a z-score
  # TODO it might if we only care about directionality
  
  if(options$mfx$value)
    MFX_WEIGHT <- 1 - mData@gene.meta$mfx.Rank[geneMask]
  else
    MFX_WEIGHT <- rep(1, n.genes)
  
  zScore <- zScore * -log10(Rfast::Pmax(matrix(1e-10, ncol = ncol(pv), nrow = nrow(pv)), as.matrix(pv)))
  
  if(options$method$value == 'diff') {
    ret <- (zScore - query) %>% t %>% as.data.table %>% `*`(MFX_WEIGHT) %>%
      .[, score := Rfast::rowsums(as.matrix(.))]
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
  ret %>% setnames(c(mData@gene.meta$gene.Name[geneMask], 'score')) %>%
    .[, is.passing := experimentMask %>% as.booltype] %>%
    .[, f.IN := experimentN / n.genes] %>%
    .[, f.OUT := pmax(0, experimentMeta$n.DE - experimentN) / experimentMeta$n.detect] %>%
    .[, ee.q := experimentMeta$ee.qScore] %>%
    .[, score := (score + 1e5) * ee.q * (1 + f.IN) / (1 + 10^(f.OUT))] %>% # abs(min(score))
    #.[, score := score / max(score)] %>%
    .[, rn := colnames(zScore)] %>%
    setorder(-score)
})

#' getTags
#' 
#' Expand the specified ontology for the specified experiments.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat, any].
#' @param rsc.IDs A list of experiment rsc IDs or NULL for everything
#' @param max.distance The maximum tree traversal distance to include
#' @param inv Whether or not to inverse the selected rscs
getTags <- function(taxa = getConfig(key = 'taxa')$value,
                    rsc.IDs = NULL, max.distance = Inf, inv = F) {
  if(is.null(rsc.IDs))
    rsc.IDs <- CACHE.BACKGROUND[[taxa]][, unique(as.character(rsc.ID))]
  
  if(inv)
    rsc.IDs <- CACHE.BACKGROUND[[taxa]][!(rsc.ID %in% rsc.IDs), unique(as.character(rsc.ID))]
  
  if(exists('TAGS', envir = globalenv()) && !is.null(TAGS[[taxa]]))
    return(TAGS[[taxa]][rsc.ID %in% rsc.IDs])
    
  CACHE.BACKGROUND[[taxa]] %>%
    .[distance <= max.distance] %>%
    .[, .(rsc.ID = as.character(rsc.ID), cf.Cat = as.character(cf.Cat), cf.BaseLongUri = as.character(cf.BaseLongUri),
          cf.ValLongUri = as.character(cf.ValLongUri), distance, reverse)] %>%
    .[rsc.ID %in% rsc.IDs] %>%
    .[as.character(cf.BaseLongUri) != as.character(cf.ValLongUri)] %>%
    .[, .(distance = mean(distance), reverse = data.table::first(reverse)),
      .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% setorder(distance, rsc.ID)
}

#' Precompute Tags
#' 
#' Get all the expanded ontology tags within a given taxon.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat, any].
#' @param mGraph The igraph form of the ontologies. Computes internally if not supplied
#' @param graphTerms The unique ontology terms. Computes internally if not supplied
precomputeTags <- function(taxa = getConfig(key = 'taxa')$value, mGraph = NULL, graphTerms = NULL) {
  if(is.null(mGraph))
    mGraph <- simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(ChildNode_Long, ParentNode_Long)]))
  if(is.null(graphTerms))
    graphTerms <- unique(ONTOLOGIES[, as.character(ChildNode_Long, ParentNode_Long)])
  
  # Tags that are simple ontology terms
  mSimpleTags <- DATA.HOLDER[[taxa]]@experiment.meta[, .(tag = unique(cf.BaseLongUri), type = 'cf.BaseLongUri'), .(rsc.ID, ee.ID)] %>%
    .[tag %in% graphTerms] %>%
    rbind(DATA.HOLDER[[taxa]]@experiment.meta[, .(tag = unique(cf.ValLongUri), type = 'cf.ValLongUri'), .(rsc.ID, ee.ID)] %>%
            .[tag %in% graphTerms])
  
  # Tags that are "bagged" (ie. have multiple tags)
  bagged <- DATA.HOLDER[[taxa]]@experiment.meta[, grepl('; ', cf.BaseLongUri, fixed = T) |
                                                  grepl('; ', cf.ValLongUri, fixed = T)]
  
  # Structured text entries (ie. non-ontology terms)
  mStructuredTags <- DATA.HOLDER[[taxa]]@experiment.meta[!bagged, .(tag = unique(cf.BaseLongUri),
                                                             type = 'cf.BaseLongUri'), .(rsc.ID, ee.ID)] %>%
    .[!(tag %in% graphTerms)] %>%
    rbind(DATA.HOLDER[[taxa]]@experiment.meta[!bagged, .(tag = unique(cf.ValLongUri),
                                                  type = 'cf.ValLongUri'), .(rsc.ID, ee.ID)] %>%
            .[!(tag %in% graphTerms)]) %>% .[, distance := 0]
  
  # Expand the bagged tags so there's one row per entry
  if(sum(bagged) == 0)
    mBagOfWords <- data.table()
  else
    mBagOfWords <- DATA.HOLDER[[taxa]]@experiment.meta[bagged, .(tag = unique(cf.BaseLongUri),
                                                                 type = 'cf.BaseLongUri'), .(rsc.ID, ee.ID)] %>%
      .[!(tag %in% graphTerms)] %>%
      rbind(DATA.HOLDER[[taxa]]@experiment.meta[bagged, .(tag = unique(cf.ValLongUri),
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
    merge(unique(ONTOLOGIES.DEFS[, .(Node_Long, Definition)]),
          by.x = 'tag', by.y = 'Node_Long', sort = F, allow.cartesian = T, all.x = T) %>%
    .[is.na(Definition), Definition := tag] %>%
    .[, .(rsc.ID, ee.ID, type, tag = as.character(Definition), distance, ID)]
  
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
    merge(DATA.HOLDER[[taxa]]@experiment.meta[, .(rsc.ID, ee.ID, cf.Cat)], by = c('rsc.ID', 'ee.ID'), all.x = T, sort = F, allow.cartesian = T) %>%
    
    # Finally reorder
    reorderTags %>%
    
    # And use some heuristics to make our corpus a little more lean
    .[, N := .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    .[distance == 0 | (N > 1 & N < 500 & distance < 5)] %>%
    .[, !'N']
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
#' set of tags
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
  
  cache[, .(rsc.ID, cf.Cat, b = cf.BaseLongUri, v = cf.ValLongUri, distance, oreverse = reverse)] %>%
    merge(vals, by.x = 'b', by.y = 'tag', sort = F, allow.cartesian = T) %>%
    merge(vals, by.x = 'v', by.y = 'tag', sort = F, allow.cartesian = T) %>%
    .[, reverse := (N.y > N.x) | (N.y == N.x & as.character(b) < as.character(v))] %>%
    .[reverse == T, c('b', 'v') := list(v, b)] %>%
    .[, reverse := ifelse(oreverse, !reverse, reverse)] %>%
    .[, .(rsc.ID = as.factor(rsc.ID),
          cf.Cat = as.factor(cf.Cat), cf.BaseLongUri = as.factor(b), cf.ValLongUri = as.factor(v),
          reverse, distance)]
}

#' Enrich
#' 
#' Given rankings (@seealso search), generate a ranking-weighted count of all terms that can be
#' derived from tags present in the experiment.
#'
#' @param rankings A named numeric (@seealso search).
#' @param options The options
#' @param inprod Whether to use the generated null distribution or not
enrich <- (function(rankings, options = getConfig(), inprod = T) {
  if(options$taxa$value == 'any') { # TODO This will include even species for which there were no homologs
    mMaps <- list(rbindlist(lapply(options$taxa$core, function(x) getTags(x, NULL))),
                  rbindlist(lapply(options$taxa$core, function(x) getTags(x, rankings$rn))))
  } else
    mMaps <- list(getTags(options$taxa$value, NULL),
                  getTags(options$taxa$value, rankings$rn))
  
  mMaps[[2]] <- mMaps[[2]] %>% merge(rankings[, .(rsc.ID = rn, score)],
                                     by = 'rsc.ID', sort = F, allow.cartesian = T)
  
  aprior <- function() {
    tmp <- mMaps[[2]][, .(distance = round(mean(distance, na.rm = T), 2),
                          A = sum(score, na.rm = T) / (1 + sum(distance, na.rm = T))),
                      .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)][, B := sum(A, na.rm = T)] %>%
      merge(mMaps[[1]][, .(C = .N / (1 + sum(distance, na.rm = T))),
                       .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],#[, D := sum(C, na.rm = T)],
            by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), sort = F, all = T, allow.cartesian = T)
    
    tmp[is.na(tmp)] <- 0
    tmp[, B := max(B)]
  }
  
  enrichTest <- function(input) {
    if(inprod && exists('NULLS', envir = globalenv())) {
      whichCols <- c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri',
                     paste0(c('M', 'S'), which(colnames(rankings) == 'score') - 1))
      
      input[, stat := A / B]
      
      #if(options$taxa$value == 'any') {
      #  for(tax in options$taxa$core) {
      #    input <- input %>%
      #      merge(NULLS[[tax]][, whichCols, with = F],
      #            by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all = T, sort = F)
      #  }
        
      #  sdCols <- grepl('S\\d+', colnames(input))
      #  mnCols <- grepl('M\\d+', colnames(input))
      #  sdPooled <- (input[, sdCols, with = F] %>% as.matrix %>% rowSums2) / sum(sdCols)
      #  mnPooled <- (input[, mnCols, with = F] %>% as.matrix %>% rowSums2) / sum(mnCols)
        
      #  input[, 1:8] %>% .[, c('M', 'S') := list(mnPooled, sdPooled)] %>% .[, stat := (stat - M) / S]
      #} else {
      input %>%
        merge(NULLS[[options$taxa$value]][, whichCols, with = F],
              by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all = T, sort = F) %>%
        .[, stat := (stat - .SD[, 9]) / .SD[, 10]]
      #}
    } else
      input[, stat := A / B] # suppressWarnings(phyper(A - 1, B, D, C + A, F)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)
  }
  
  aprior() %>% enrichTest %>% .[is.finite(stat)] %>% setorder(-stat)
})
