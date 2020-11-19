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
  # GEEQ <- options$geeq
  METHOD <- options$method
  
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
  pv[is.na(pv)] <- 1
  colFilter <- Rfast::colsums(pv < P_THRESHOLD) > 0
  
  # Number of DEGs for experiments that pass thresholds
  # geeq <- DATA.HOLDER[[taxa]]@experiment.meta$ee.qScore[colFilter]
  # n.DE.exp <- DATA.HOLDER[[taxa]]@experiment.meta$n.DE[colFilter]
  
  # logFCs for only the GOI/EOI and maintain structure.
  logFC <- DATA.HOLDER[[taxa]]@data$fc[rowFilter, colFilter]
  logFC[is.na(logFC)] <- 0
  
  zscore <- DATA.HOLDER[[taxa]]@data$pvz[rowFilter, colFilter]
  zscore[is.na(zscore)] <- 0
  
  if(n.genes == 1) {
    logFC <- t(logFC)
    zscore <- t(zscore)
  }
  
  logFC <- logFC %>% as.data.table
  zscore <- zscore %>% as.data.table
  
  # Only retain experiments that have at least one of the GOI with threshold_u > abs(logFC) > threshold_l
  if(FC_L_THRESHOLD != 0 || FC_U_THRESHOLD != 0) {
    exFilter <- abs(logFC) > FC_L_THRESHOLD
    
    if(FC_U_THRESHOLD != FC_L_THRESHOLD)
      exFilter <- exFilter & (abs(logFC) < FC_U_THRESHOLD)
    
    colFilter <- Rfast::colsums(exFilter) > 0
    logFC <- logFC %>% .[, colFilter, with = F]
    # n.DE.exp <- n.DE.exp[colFilter]
    # geeq <- geeq[colFilter]
  }
  
  # Data Processing ---------------------------------------------------------
  if(verbose)
    advanceProgress('Ranking')
  
  #n.DE.exp[is.na(n.DE.exp)] <- 0
  #tf <- t(t(tf) / c(Rfast::Log(2 + n.DE.exp))) %>% as.data.table # Weight by number of DEGs in this experiment
  
  #geeq[is.na(geeq)] <- 0
  #if(GEEQ)
  #  tf <- t(t(tf) * (1 + Rfast::Pmax(geeq, 0))) %>% as.data.table # Weight by experiment quality score
  
  query <- Rfast::rowMaxs(zscore %>% as.matrix, value = T)
  directions <- ifelse(Rfast::colmeans(zscore %>% as.matrix) < 0, F, T)
  
  if(MFX)
    MFX_WEIGHT <- 1 - DATA.HOLDER[[taxa]]@gene.meta$mfx.Rank[rowFilter]
  else
    MFX_WEIGHT <- rep(1, n.genes)

  if(METHOD == 'zscore') {
    (zscore / query) %>% t %>% as.data.table %>%
      .[, score := Rfast::rowsums(as.matrix(. * MFX_WEIGHT)) / sum(MFX_WEIGHT)] %>%
      setnames(c(DATA.HOLDER[[taxa]]@gene.meta$gene.Name[rowFilter], 'score')) %>%
      .[, direction := directions] %>%
      .[, rn := colnames(zscore)] %>%
      .[, score := scale(score, F)] %>%
      setorder(-score)
  } else if(METHOD == 'mvsm') {
    idf <- log(1 / MFX_WEIGHT) + 1
    zscore[, query := query]
    tfidf <- zscore * idf
    
    # Extract the query
    query <- tfidf[, query]
    tfidf[, query := NULL]
    
    # Cosine similarity
    tfidf <- as.matrix(tfidf)
    
    cross_x <- crossprod(query)
    scores <- sapply(1:ncol(tfidf), function(i) {
      crossprod(query, tfidf[, i]) / sqrt(cross_x * crossprod(tfidf[, i]))
    })
    
    ret <- rbind(tfidf, scores) %>% t %>% as.data.table %>%
      setnames(c(DATA.HOLDER[[taxa]]@gene.meta$gene.Name[rowFilter], 'score')) %>%
      .[, direction := directions] %>%
      .[, rsc.ID := colnames(tfidf)] %>%
      .[, score := scale(score, F)]
  
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
    
  CACHE.BACKGROUND[[taxa]] %>%
    .[distance <= max.distance] %>%
    .[, .(rsc.ID = as.character(rsc.ID), cf.Cat = as.character(cf.Cat), cf.BaseLongUri = as.character(cf.BaseLongUri),
          cf.ValLongUri = as.character(cf.ValLongUri), distance, reverse)] %>%
    .[rsc.ID %in% rsc.IDs & !(cf.BaseLongUri %in% graphTerms) & !(cf.ValLongUri %in% graphTerms)] %>%
    .[as.character(cf.BaseLongUri) != as.character(cf.ValLongUri)] %>%
    .[, .(distance = mean(distance), reverse = first(reverse)),
      .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% setorder(distance, rsc.ID)
}

#' Precompute Tags
#' 
#' Get all the expanded ontology tags within a given taxon.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
precomputeTags <- function(taxa = getOption('app.taxa')) {
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
      .[, reverse := (N.y > N.x) | (N.y == N.x & b < v)] %>%
      .[reverse == T, c('b', 'v') := list(v, b)] %>%
      .[, .(rsc.ID = as.factor(rsc.ID), ee.ID,
            cf.Cat = as.factor(cf.Cat), cf.BaseLongUri = as.factor(b), cf.ValLongUri = as.factor(v),
            reverse, distance)]
  }
  
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
  mBagOfWords <- DATA.HOLDER[[taxa]]@experiment.meta[bagged, .(tag = unique(cf.BaseLongUri),
                                                                type = 'cf.BaseLongUri'), .(rsc.ID, ee.ID)] %>%
    .[!(tag %in% graphTerms)] %>%
    rbind(DATA.HOLDER[[taxa]]@experiment.meta[bagged, .(tag = unique(cf.ValLongUri),
                                                         type = 'cf.ValLongUri'), .(rsc.ID, ee.ID)] %>%
            .[!(tag %in% graphTerms)]) %>%
    na.omit %>%
    .[, lapply(.SD, function(x) parseListEntry(as.character(x))), .(rsc.ID, ee.ID, type)] %>%
    .[, ID := 1:length(tag), .(rsc.ID, ee.ID, type)]
  
  # Compute ontology expansions on ontology terms
  # Applies to both simple tags and bag of word tags that expanded into ontology terms
  mComputedTags <- rbindlist(lapply(union(mSimpleTags[, unique(tag)], mBagOfWords[tag %in% graphTerms, tag]), function(uri) {
    tag <- igraph::subcomponent(mGraph, uri, 'out')
    distance <- igraph::distances(mGraph, uri, tag)
    data.table(startTag = uri, tag = names(tag), distance = c(distance))
  }))
  
  # Simple (ontology) tags get associated with their parents from mComputedTags
  mTags <- mSimpleTags %>% merge(mComputedTags[startTag %in% mSimpleTags[, unique(tag)]],
                                 by.x = 'tag', by.y = 'startTag', sort = F, allow.cartesian = T) %>%
    .[, c('tag', 'tag.y') := list(tag.y, NULL)] %>% .[, ID := NA] %>%
    
    # Structured tags just get inserted
    rbind(mStructuredTags[, ID := NA] %>% .[, .(tag, rsc.ID, ee.ID, type, distance, ID)]) %>%
    
    # Ontology-expanded bag of word tags get associated with their parents from mComputedTags
    rbind(
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
  mTags %>% .[!is.na(ID) & distance < 1, expand.grid(aggregate(.SD[, tag], by = list(.SD[, ID]), FUN = list)[[-1]]) %>%
                apply(1, paste0, collapse = '; ') %>% unique, .(rsc.ID, ee.ID, type)] %>%
    .[, c('tag', 'V1') := list(V1, NULL)] %>% rbind(mTags[is.na(ID), .(rsc.ID, ee.ID, type, tag)]) %>%
    
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
  
  mMaps <- list(getTags(taxa, scope, NULL, options$distance),
                getTags(taxa, scope, rankings$rn, options$distance))
  
  # TODO We're losing tags?
  
  mMaps[[2]] <- mMaps[[2]] %>%
    merge(rankings[, .(rsc.ID = rn, rank = ifelse(options$method == 'mvsm', score^2, score), direction)],
          by = 'rsc.ID', sort = F, allow.cartesian = T)
  
  aprior <- function() {
    tmp <- mMaps[[2]][, .(A = sum(rank) / (1 + sum(distance))),
                      .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)][, B := sum(A)] %>%
      merge(mMaps[[1]][, .(C = .N / (1 + sum(distance))),
                       .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)][, D := sum(C)],
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
  
  PATTERN <- 'http://purl.org/commons/record/ncbi_gene/(.*)'
  ret %>% .[, pv.fisher.adj := p.adjust(pv.fisher, 'BH')] %>% setorder(pv.fisher)
}
