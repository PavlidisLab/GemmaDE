#' Search
#' 
#' Uses the M-VSM to sort experiments that show at least one of the genes as DE.
#'
#' @param genes A list of Entrez Gene IDs (ie. 1, 22, 480) as characters
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
#' @param options Optional extra parameters to pass, such as:
#' * pv: An FDR cutoff (default: 0.05)
#' * fc.lower / fc.upper: Upper and lower logFC thresholds (default: 0 / 10)
#' * mfx: Whether or not to scale by gene multifunctionality
#' * geeq: Whether or not to scale by GEEQ score
#' @param verbose Whether or not to print messages to the progress bar (if it can be found)
search <- function(genes, taxa = getOption('app.taxa'), options = getOption('app.all_options'), verbose = T) {
  if(verbose)
    advanceProgress('Loading experiments')
  
  P_THRESHOLD <- options$pv
  FC_L_THRESHOLD <- options$fc.lower
  FC_U_THRESHOLD <- options$fc.upper
  MFX <- options$mfx
  GEEQ <- options$geeq
  METHOD <- options$method
  REQ_ALL <- options$reqall
  MEANVAL <- options$meanval
  
  # Data Extraction ---------------------------------------------------------
  
  if(taxa == 'any') {
    mData <- new('EData')
    taxOptions <- Filter(function(x) !(x %in% c('artificial', 'dope', 'any')), getOption('app.all_taxa'))
    taxIDs <- c(human = 9606, mouse = 10090, rat = 10116)
    
    for(x in taxOptions) {
      message(x)
      ID <- paste0(unname(taxIDs[x]), '_ID')
      IDs <- genes[, ..ID] %>% unlist %>% unname
      
      gMeta <- DATA.HOLDER[[x]]@gene.meta %>% copy %>% .[, I := .I]
      gMap <- merge(genes[, c('key', ID), with = F],
                    gMeta %>% .[entrez.ID %in% IDs, .(entrez.ID, I)],
                    by.x = ID, by.y = 'entrez.ID', all = T) %>%
        .[, c(ID) := list(NULL)]
      
      pvData <- DATA.HOLDER[[x]]@data$adj.pv[gMap$I, ] %>% `rownames<-`(gMap$key)
      pvData[is.na(pvData)] <- 1
      fcData <- DATA.HOLDER[[x]]@data$fc[gMap$I, ] %>% `rownames<-`(gMap$key)
      pvzData <- DATA.HOLDER[[x]]@data$pvz[gMap$I, ] %>% `rownames<-`(gMap$key)
      meanvalData <- DATA.HOLDER[[x]]@data$meanval[gMap$I, ] %>% `rownames<-`(gMap$key)
      
      mData@experiment.meta <- rbind(mData@experiment.meta, DATA.HOLDER[[x]]@experiment.meta)
      
      if(is.null(mData@data$adj.pv)) {
        mData@data$adj.pv <- pvData
        mData@data$fc <- fcData
        mData@data$pvz <- pvzData
        mData@data$meanval <- meanvalData
        
        mData@gene.meta <- gMeta[gMap$I, ] %>% .[, entrez.ID := gMap$key]
      } else {
        mData@data$adj.pv <- merge(mData@data$adj.pv, pvData, by = 0, all = T) %>%
          `rownames<-`(.[, 1]) %>% .[, -1]
        mData@data$fc <- merge(mData@data$fc, fcData, by = 0, all = T) %>%
          `rownames<-`(.[, 1]) %>% .[, -1]
        mData@data$pvz <- merge(mData@data$pvz, pvzData, by = 0, all = T) %>%
          `rownames<-`(.[, 1]) %>% .[, -1]
        mData@data$meanval <- merge(mData@data$meanval, meanvalData, by = 0, all = T) %>%
          `rownames<-`(.[, 1]) %>% .[, -1]
        
        mData@gene.meta$gene.Name <- coalesce(mData@gene.meta$gene.Name, gMeta[gMap$I, gene.Name])
        mData@gene.meta$mfx.Rank <- rowMeans2(cbind(mData@gene.meta$mfx.Rank, gMeta[gMap$I, mfx.Rank]),
                                              na.rm = T)
      }
    }
    
    genes <- genes$key
    rowFilter <- which(!is.na(mData@gene.meta$mfx.Rank))
  } else {
    mData <- DATA.HOLDER[[taxa]]

    # Only retain GOI
    rowFilter <- which(mData@gene.meta$entrez.ID %in% genes)
  }
  
  n.genes <- length(rowFilter)
  if(n.genes == 0) {
    setProgress(T)
    return(NULL)
  }
  
  # P-values for only the GOI
  pv <- mData@data$adj.pv[rowFilter, ]
  
  if(n.genes == 1)
    pv <- t(pv)
  
  # Only retain experiments that have at least one of the GOI as DE (pv < threshold)
  if(taxa != 'any')
    pv[is.na(pv)] <- 1 # Turns out this is significantly faster on smaller matrices
  
  colFilter <- Rfast::colsums(pv <= P_THRESHOLD) >= min(REQ_ALL, n.genes)
  
  # logFCs for only the GOI/EOI and maintain structure.
  logFC <- mData@data$fc[rowFilter, colFilter]
  zscore <- mData@data$pvz[rowFilter, colFilter]
  mx <- mData@data$meanval[rowFilter, colFilter]
  
  if(n.genes == 1) {
    logFC <- t(logFC)
    zscore <- t(zscore)
    mx <- t(mx)
  }
  
  logFC <- logFC %>% as.data.table
  zscore <- zscore %>% as.data.table
  mx <- mx %>% as.data.table
  
  # Only retain experiments that have at least one of the GOI with threshold_u > abs(logFC) > threshold_l
  if(FC_L_THRESHOLD != 0 || FC_U_THRESHOLD != 0) {
    exFilter <- abs(logFC) >= FC_L_THRESHOLD
    exFilter[is.na(exFilter)] <- F
    
    if(FC_U_THRESHOLD != FC_L_THRESHOLD)
      exFilter <- exFilter & (abs(logFC) <= FC_U_THRESHOLD)
    exFilter[is.na(exFilter)] <- F
    
    colFilter <- Rfast::colsums(exFilter) >= min(REQ_ALL, n.genes)
    
    if(sum(colFilter) == 0) {
      setProgress(T)
      return(NULL)
    }
    
    logFC <- logFC[, colFilter, with = F]
    zscore <- zscore[, colFilter, with = F]
    mx <- mx[, colFilter, with = F]
  }
  
  # Data Processing ---------------------------------------------------------
  if(verbose)
    advanceProgress('Ranking')
  
  # Number of DEGs for experiments that pass thresholds
  mFilter <- mData@experiment.meta %>% as.data.frame %>%
    `rownames<-`(.[, 'rsc.ID']) %>% .[colnames(zscore), c('ee.qScore', 'n.DE', 'ee.Scale')]
  
  geeq <- mFilter$ee.qScore
  n.DE.exp <- mFilter$n.DE
  scales <- mFilter$ee.Scale
  
  n.DE.exp[is.na(n.DE.exp)] <- 0
  geeq[is.na(geeq)] <- -1

  # Put everything on a linear scale
  mx <- as.matrix(mx)
  mx[, scales == 'LOG2'] <- 2^mx[, scales == 'LOG2']
  mx[, scales == 'LOG10'] <- 10^mx[, scales == 'LOG10']
  mx[, scales == 'LN'] <- exp(1)^mx[, scales == 'LN']
  # What to do with LINEAR, UNSCALED, OTHER, COUNT?
  
  # TODO Capping randomly at 10k
  mx <- Rfast::Pmin(matrix(1e4, nrow(mx), ncol(mx)), mx)
  
  mx[is.na(mx) | is.infinite(mx)] <- 0
  
  if(!GEEQ)
    geeq <- 2
  else
    geeq <- pmax(geeq + 1, 0)
  
  geeq <- c(geeq / Rfast::Log(2 + n.DE.exp)^(5/4)) / 2
  
  logFC[is.na(logFC)] <- 0
  directions <- ifelse(Rfast::colmeans(logFC %>% as.matrix) < 0, F, T)
  
  # TODO This should have an equivalent in artificial.
  # TODO Justify weighting
  if(MEANVAL && taxa != 'artificial')
    zscore <- zscore * (1 + mx^(1/4))
  
  zscore[is.na(zscore)] <- 0
  
  query <- Rfast::rowMaxs(zscore %>% as.matrix, value = T)
  
  if(MFX)
    MFX_WEIGHT <- 1 - mData@gene.meta$mfx.Rank[rowFilter]
  else
    MFX_WEIGHT <- rep(1, n.genes)
  
  if(METHOD == 'zscore') {
    (zscore / query) %>% t %>% as.data.table %>%
      .[, score := Rfast::rowsums(as.matrix(. * MFX_WEIGHT)) / sum(MFX_WEIGHT)] %>%
      setnames(c(mData@gene.meta$gene.Name[rowFilter], 'score')) %>%
      .[, direction := directions] %>%
      .[, rn := colnames(zscore)] %>%
      .[is.nan(score), score := 0] %>%
      .[, score := 1 + scale(c(geeq) * score, F)] %>%
      setorder(-score)
  } else if(METHOD == 'mvsm') {
    #idf <- Rfast::Log(1 / MFX_WEIGHT) + 1 # Gene-wise IDF
    #idf2 <- c(geeq) # Experiment-wise IDF
    
    #tfidf <- t(t(zscore * idf) * idf2) %>% as.matrix
    
    #query <- Rfast::rowMaxs(tfidf, value = T)
    
    idf <- Rfast::Log(1 / MFX_WEIGHT) + 1
    zscore[, query := query]
    tfidf <- zscore * idf
    
    # Extract the query
    query <- tfidf[, query]
    tfidf[, query := NULL]
    
    # Cosine similarity
    tfidf <- as.matrix(tfidf)
    
    cross_x <- crossprod(query)
    scores <- sapply(1:ncol(tfidf), function(i) {
      cross_q <- crossprod(query, tfidf[, i])
      if(cross_q == 0) 0
      else cross_q / sqrt(cross_x * crossprod(tfidf[, i]))
    })
    scores[is.nan(scores)] <- 0
    
    ret <- rbind(tfidf, scores) %>% t %>% as.data.table %>%
      setnames(c(mData@gene.meta$gene.Name[rowFilter], 'score')) %>%
      .[, direction := directions] %>%
      .[, rsc.ID := colnames(tfidf)] %>%
      .[, score := 1 + scale(c(geeq) * score, F)]
  
    ret %>% setnames('rsc.ID', 'rn') %>% setorder(-score)
  }
}

#' getTags
#' 
#' Expand the specified ontology for the specified experiments. This could theoretically be done
#' before and cached like this, but it seems to take up much more space in RAM.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
#' @param scope The ontology scope.
#' @param rsc.IDs A list of experiment rsc IDs or NULL for everything
#' @param max.distance The maximum tree traversal distance to include
getTags <- function(taxa = getOption('app.taxa'), scope = getOption('app.ontology'),
                    rsc.IDs = NULL, max.distance = Inf) {
  if(is.null(rsc.IDs))
    rsc.IDs <- CACHE.BACKGROUND[[taxa]][, unique(as.character(rsc.ID))]
  
  graphTerms <- ONTOLOGIES.DEFS[!(OntologyScope %in% scope), as.character(Definition)]
  
  if(exists('TAGS') && !is.null(TAGS[[taxa]]))
    return(TAGS[[taxa]][rsc.ID %in% rsc.IDs & !(cf.BaseLongUri %in% graphTerms) & !(cf.ValLongUri %in% graphTerms)])
    
  CACHE.BACKGROUND[[taxa]] %>%
    .[distance <= max.distance] %>%
    .[, .(rsc.ID = as.character(rsc.ID), cf.Cat = as.character(cf.Cat), cf.BaseLongUri = as.character(cf.BaseLongUri),
          cf.ValLongUri = as.character(cf.ValLongUri), distance, reverse)] %>%
    .[rsc.ID %in% rsc.IDs & !(cf.BaseLongUri %in% graphTerms) & !(cf.ValLongUri %in% graphTerms)] %>%
    .[as.character(cf.BaseLongUri) != as.character(cf.ValLongUri)] %>%
    .[, .(distance = mean(distance), reverse = data.table::first(reverse)),
      .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% setorder(distance, rsc.ID)
}

#' Precompute Tags
#' 
#' Get all the expanded ontology tags within a given taxon.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
precomputeTags <- function(taxa = getOption('app.taxa')) {
  mGraph <- simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(ChildNode_Long, ParentNode_Long)]))
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
            .[!(tag %in% graphTerms)]) %>% na.omit %>% .[, distance := 0]
  
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
      na.omit %>%
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
  
  if(nrow(mBagOfWords) > 0 & nrow(mComputedTags) > 0)
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
    mTagsExpanded <- mTags %>% .[!is.na(ID) & distance < 1, expand.grid(aggregate(.SD[, tag], by = list(.SD[, ID]), FUN = list)[[-1]]) %>%
                                   apply(1, paste0, collapse = '; ') %>% unique, .(rsc.ID, ee.ID, type)] %>%
    .[, c('tag', 'V1') := list(V1, NULL)] %>% rbind(mTags[is.na(ID), .(rsc.ID, ee.ID, type, tag)])
  else
    mTagsExpanded <- mTags
  
  mTagsExpanded %>%
    
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
    reorderTags
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

#' Enrich
#' 
#' Given rankings (@seealso search), generate a ranking-weighted count of all terms that can be
#' derived from tags present in the experiment (both cf.ValLongUri and ee.TagLongUri).
#'
#' @param rankings A named numeric (@seealso search).
#' @param scope The ontology scope.
#' @param options The options
#' @param verbose Whether or not to print messages to the progress bar (if it can be found)
enrich <- function(rankings, taxa = getOption('app.taxa'), scope = getOption('app.ontology'),
                   options = getOption('app.all_options'), verbose = T) {
  if(verbose)
    advanceProgress('Looking up ontology terms')
  
  if(taxa == 'any') {
    taxOptions <- Filter(function(x) !(x %in% c('artificial', 'dope', 'any')), getOption('app.all_taxa'))
    
    mMaps <- list(rbindlist(lapply(taxOptions, function(x) getTags(x, scope, NULL, options$distance))),
                  rbindlist(lapply(taxOptions, function(x) getTags(x, scope, rankings$rn, options$distance))))
  } else
    mMaps <- list(getTags(taxa, scope, NULL, options$distance),
                  getTags(taxa, scope, rankings$rn, options$distance))
  
  mMaps[[2]] <- mMaps[[2]] %>%
    merge(rankings[, .(rsc.ID = rn, score, direction)],
          by = 'rsc.ID', sort = F, allow.cartesian = T)
  
  aprior <- function() {
    tmp <- mMaps[[2]][, .(A = sum(score, na.rm = T) / (1 + sum(distance, na.rm = T))),
                      .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)][, B := sum(A, na.rm = T)] %>%
      merge(mMaps[[1]][, .(C = .N / (1 + sum(distance, na.rm = T))),
                       .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)][, D := sum(C, na.rm = T)],
            by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), sort = F, all = T, allow.cartesian = T)
    
    tmp[is.na(tmp)] <- 0
    tmp[, B := max(B)]
  }
  
  fisher <- function(input) {
    if(verbose)
      advanceProgress('Running tests')
    
    input[, pv.fisher := suppressWarnings(phyper(A - 1, B, D, C + A, F)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
  }
  
  ret <- aprior() %>% fisher
  
  ret %>% setorder(pv.fisher) %>% .[C >= max(0, options$min.tags)]
}
