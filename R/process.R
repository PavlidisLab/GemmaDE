#' Weighted correlation
#'
#' @param x A vector or matrix (whose columns will be correlated)
#' @param y A vector with compatible dimensions of x to correlate against
#' @param w A vector of weights, the same length as y
#'
#' @return The weighted Pearson correlation of x (or the columns of x) and y
cor.wt <- function(x, y, w = rep(1, length(y))) {
  if (is.matrix(x)) {
    apply(x, 2, cor.wt, y = y, w = w)
  } else {
    sw <- sum(w)
    
    mx <- sum(w * x) / sw
    my <- sum(w * y) / sw
    
    (sum(w * (x - mx) * (y - my)) / sw) / sqrt((sum(w * (x - mx)^2) / sw) * (sum(w * (y - my)^2) / sw))
  }
}

#' getTags
#'
#' Get tags for the specified experiments within a specified distance.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
#' @param rsc.IDs A list of experiment rsc IDs or NULL for everything
#' @param max.distance The maximum tree traversal distance to include
#' @param inv Whether or not to inverse the selected rscs
#' @param CACHE A cache to use. If null, uses the global CACHE.BACKGROUND
getTags <- function(taxa,
                    rsc.IDs = NULL, max.distance = Inf, inv = F, CACHE = NULL) {
  if (length(taxa) > 1) {
    return(data.table::rbindlist(lapply(taxa, getTags, rsc.IDs, max.distance, inv, CACHE)))
  }
  
  if (is.null(CACHE)) {
    CACHE <- CACHE.BACKGROUND
  }
  
  if (is.null(rsc.IDs)) {
    rsc.IDs <- CACHE[[taxa]][, unique(as.character(rsc.ID))]
  }
  
  if (inv) {
    rsc.IDs <- CACHE[[taxa]][!(rsc.ID %in% rsc.IDs), unique(as.character(rsc.ID))]
  }
  
  CACHE[[taxa]][distance <= max.distance & rsc.ID %in% rsc.IDs]
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
          by.x = "cf.BaseLongUri", by.y = "cf.ValLongUri", all = T, sort = F, allow.cartesian = T
    ) %>%
    .[is.na(N.x), N.x := 0] %>%
    .[is.na(N.y), N.y := 0] %>%
    .[, N := N.x + N.y, cf.BaseLongUri] %>%
    .[, .(tag = cf.BaseLongUri, N)]
  
  cache[, .(rsc.ID, ee.ID, cf.Cat, b = cf.BaseLongUri, v = cf.ValLongUri, distance)] %>%
    merge(vals, by.x = "b", by.y = "tag", sort = F, allow.cartesian = T) %>%
    merge(vals, by.x = "v", by.y = "tag", sort = F, allow.cartesian = T) %>%
    .[, reverse := (N.y > N.x) | (N.y == N.x & as.character(b) < as.character(v))] %>%
    .[reverse == T, c("b", "v") := list(v, b)] %>%
    .[, .(
      rsc.ID = as.factor(rsc.ID), ee.ID,
      cf.Cat = as.factor(cf.Cat), cf.BaseLongUri = as.factor(b), cf.ValLongUri = as.factor(v),
      reverse, distance
    )]
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
          by.x = "cf.BaseLongUri", by.y = "cf.ValLongUri", all = T, sort = F, allow.cartesian = T
    ) %>%
    .[is.na(N.x), N.x := 0] %>%
    .[is.na(N.y), N.y := 0] %>%
    .[, N := N.x + N.y, cf.BaseLongUri] %>%
    .[, .(tag = cf.BaseLongUri, N)]
  
  cache[, .(rsc.ID, ee.ID, cf.Cat, b = cf.BaseLongUri, v = cf.ValLongUri, distance, oreverse = reverse)] %>%
    merge(vals, by.x = "b", by.y = "tag", sort = F, allow.cartesian = T) %>%
    merge(vals, by.x = "v", by.y = "tag", sort = F, allow.cartesian = T) %>%
    .[, reverse := (N.y > N.x) | (N.y == N.x & as.character(b) < as.character(v))] %>%
    .[reverse == T, c("b", "v") := list(v, b)] %>%
    .[, reverse := ifelse(oreverse, !reverse, reverse)] %>%
    .[, .(
      rsc.ID = as.factor(rsc.ID), ee.ID,
      cf.Cat = as.factor(cf.Cat), cf.BaseLongUri = as.factor(b), cf.ValLongUri = as.factor(v),
      reverse, distance
    )]
}

# TODO this could be merged with reorderTags2
reorderTags3 <- function(data) {
  vals <- data[, .N, cf.BaseLongUri] %>%
    merge(data[, .N, cf.ValLongUri],
          by.x = "cf.BaseLongUri", by.y = "cf.ValLongUri", all = T, sort = F, allow.cartesian = T
    ) %>%
    .[is.na(N.x), N.x := 0] %>%
    .[is.na(N.y), N.y := 0] %>%
    .[, N := N.x + N.y, cf.BaseLongUri] %>%
    .[, .(tag = cf.BaseLongUri, N)]
  
  data %>%
    data.table::copy() %>%
    data.table::setnames(c("cf.BaseLongUri", "cf.ValLongUri"), c("b", "v")) %>%
    merge(vals, by.x = "b", by.y = "tag", sort = F, allow.cartesian = T) %>%
    merge(vals, by.x = "v", by.y = "tag", sort = F, allow.cartesian = T) %>%
    .[, reverse := (N.y > N.x) | (N.y == N.x & as.character(b) < as.character(v))] %>%
    .[reverse == T, c("b", "v") := list(v, b)] %>%
    .[, c("cf.Cat", "cf.BaseLongUri", "cf.ValLongUri") :=
        list(as.factor(cf.Cat), as.factor(b), as.factor(v))] %>%
    .[, !c("reverse", "b", "v", "N.x", "N.y")]
}

#' Precompute Tags
#'
#' Get all the expanded ontology tags within a given taxon.
#'
#' @param taxa A taxa scope. Can be one of [human, mouse, rat].
#' @param mGraph The igraph form of the ontologies. Computes internally if not supplied
#' @param graphTerms The unique ontology terms. Computes internally if not supplied
#' @param DATA Data to use. If null, uses the global DATA.HOLDER
#' @param POST Whether to do "post-pre-processing" like reordering and trimming
precomputeTags <- function(taxa, mGraph = NULL, graphTerms = NULL,
                           DATA = NULL, ONTOLOGIES = NULL, ONTOLOGIES.DEFS = NULL , POST = T) {
  if (is.null(DATA)) {
    DATA <- DATA.HOLDER
  }
  
  if(is.null(ONTOLOGIES)){
    ONTOLOGIES = parent.frame()$ONTOLOGIES
  }
  
  if(is.null(ONTOLOGIES.DEFS)){
    ONTOLOGIES.DEFS = parent.frame()$ONTOLOGIES.DEFS
  }
  
  if (is.null(mGraph)) {
    mGraph <- igraph::simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(as.character(ChildNode_Long), as.character(ParentNode_Long))]))
  }
  if (is.null(graphTerms)) {
    graphTerms <- unique(ONTOLOGIES[, as.character(ChildNode_Long, ParentNode_Long)])
  }
  
  # Tags that are simple ontology terms
  mSimpleTags <- DATA[[taxa]]@experiment.meta[, .(tag = unique(as.character(cf.BaseLongUri)), type = "cf.BaseLongUri"), .(rsc.ID, ee.ID)] %>%
    .[tag %in% graphTerms] %>%
    rbind(DATA[[taxa]]@experiment.meta[, .(tag = unique(as.character(cf.ValLongUri)), type = "cf.ValLongUri"), .(rsc.ID, ee.ID)] %>%
            .[tag %in% graphTerms])
  
  # Tags that are "bagged" (ie. have multiple tags)
  bagged <- DATA[[taxa]]@experiment.meta[, grepl("; ", as.character(cf.BaseLongUri), fixed = T) |
                                           grepl("; ", as.character(cf.ValLongUri), fixed = T)]
  
  # Structured text entries (ie. non-ontology terms)
  mStructuredTags <- DATA[[taxa]]@experiment.meta[!bagged, .(
    tag = unique(as.character(cf.BaseLongUri)),
    type = "cf.BaseLongUri"
  ), .(rsc.ID, ee.ID)] %>%
    .[!(tag %in% graphTerms)] %>%
    rbind(DATA[[taxa]]@experiment.meta[!bagged, .(
      tag = unique(as.character(cf.ValLongUri)),
      type = "cf.ValLongUri"
    ), .(rsc.ID, ee.ID)] %>%
      .[!(tag %in% graphTerms)]) %>%
    .[, distance := 0]
  
  # Expand the bagged tags so there's one row per entry
  if (sum(bagged) == 0) {
    mBagOfWords <- data.table()
  } else {
    mBagOfWords <- DATA[[taxa]]@experiment.meta[bagged, .(
      tag = unique(as.character(cf.BaseLongUri)),
      type = "cf.BaseLongUri"
    ), .(rsc.ID, ee.ID)] %>%
      # .[!(tag %in% graphTerms)] %>%
      rbind(DATA[[taxa]]@experiment.meta[bagged, .(
        tag = unique(as.character(cf.ValLongUri)),
        type = "cf.ValLongUri"
      ), .(rsc.ID, ee.ID)]) %>%
      # .[!(tag %in% graphTerms)]) %>%
      .[, lapply(.SD, function(x) parseListEntry(as.character(x))), .(rsc.ID, ee.ID, type)] %>%
      .[, ID := 1:length(tag), .(rsc.ID, ee.ID, type)]
  }
  
  if (nrow(mBagOfWords) > 0) {
    mComputable <- union(mSimpleTags[, unique(tag)], mBagOfWords[tag %in% graphTerms, tag])
  } else {
    mComputable <- mSimpleTags[, unique(tag)]
  }
  
  # Compute ontology expansions on ontology terms
  # Applies to both simple tags and bag of word tags that expanded into ontology terms
  mComputedTags <- data.table::rbindlist(lapply(mComputable, function(uri) {
    tag <- igraph::subcomponent(mGraph, uri, "out")
    distance <- igraph::distances(mGraph, uri, tag)
    data.table(startTag = uri, tag = names(tag), distance = c(distance))
  }))
  
  if (nrow(mSimpleTags) > 0 && nrow(mComputedTags) > 0) {
    # Simple (ontology) tags get associated with their parents from mComputedTags
    mTags <- mSimpleTags %>%
      merge(mComputedTags, by.x = "tag", by.y = "startTag", sort = F, allow.cartesian = T) %>%
      .[, c("tag", "tag.y") := list(tag.y, NULL)] %>%
      .[, ID := NA]
  } else if (nrow(mSimpleTags) > 0) {
    mTags <- mSimpleTags[, .(tag, rsc.ID, ee.ID, type, distance = 0, ID = NA)]
  } else {
    mTags <- data.table()
  }
  
  # Structured tags just get inserted
  if (nrow(mStructuredTags) > 0) {
    mTags <- mTags %>% rbind(mStructuredTags[, .(tag, rsc.ID, ee.ID, type, distance, ID = NA)])
  }
  
  if (nrow(mBagOfWords) > 0 && nrow(mComputedTags) > 0) {
    # Ontology-expanded bag of word tags get associated with their parents from mComputedTags
    mTags <- mTags %>% rbind(
      mBagOfWords %>%
        merge(mComputedTags[startTag %in% mBagOfWords[, unique(tag)]],
              by.x = "tag", by.y = "startTag", all = T, sort = F, allow.cartesian = T
        ) %>%
        .[is.na(tag.y), c("tag.y", "distance") := list(tag, 0)] %>%
        .[, c("tag", "tag.y") := list(tag.y, NULL)]
    )
  }
  
  mTags <- mTags %>%
    merge(unique(ONTOLOGIES.DEFS[, .(Node_Long = as.character(Node_Long), Definition = as.character(Definition))]),
          by.x = "tag", by.y = "Node_Long", sort = F, allow.cartesian = T, all.x = T
    ) %>%
    .[is.na(Definition), Definition := tag] %>%
    .[, .(rsc.ID, ee.ID, type, tag = Definition, distance, ID)] %>%
    unique()
  
  # Do a grid expansion of bagged terms, maintaining indexed ordering
  # TODO This only works with 0 dist tags. Increasing the dist threshold is too memory hungry
  if (nrow(mBagOfWords) > 0 && nrow(mTags[!is.na(ID) & distance < 1]) > 0) {
    mTagsExpanded <- mTags %>% rbind(
      .[
        !is.na(ID) & distance < 1, expand.grid(aggregate(tag, by = list(ID), FUN = list)[[-1]]) %>% apply(1, paste0, collapse = "; "),
        .(rsc.ID, ee.ID, type)
      ] %>%
        cbind(
          distance =
            mTags %>%
            .[
              !is.na(ID) & distance < 1, expand.grid(aggregate(distance, by = list(ID), FUN = list)[[-1]]) %>% apply(1, sum),
              .(rsc.ID, ee.ID, type)
            ] %>% .[, V1]
        ) %>% data.table::setnames("V1", "tag") %>%
        .[, .(tag, rsc.ID, ee.ID, type, distance, ID = NA)]
    )
  } else {
    mTagsExpanded <- mTags
  }
  
  mTagsExpanded %>%
    .[, .(rsc.ID, ee.ID, type, tag)] %>%
    # Now expand these expanded terms
    .[, expand.grid(
      cf.BaseLongUri = .SD[type == "cf.BaseLongUri", tag],
      cf.ValLongUri = .SD[type == "cf.ValLongUri", tag]
    ), .(rsc.ID, ee.ID)] %>%
    unique() %>%
    # Add back distances
    merge(mTags[type == "cf.BaseLongUri", .(rsc.ID, ee.ID, tag, distance)],
          by.x = c("rsc.ID", "ee.ID", "cf.BaseLongUri"), by.y = c("rsc.ID", "ee.ID", "tag"), all.x = T, sort = F, allow.cartesian = T
    ) %>%
    merge(mTags[type == "cf.ValLongUri", .(rsc.ID, ee.ID, tag, distance)],
          by.x = c("rsc.ID", "ee.ID", "cf.ValLongUri"), by.y = c("rsc.ID", "ee.ID", "tag"), all.x = T, sort = F, allow.cartesian = T
    ) %>%
    # Recombinants have unknown distances. Let's just make them 0 since that's
    # the memory prohibitive limit
    .[is.na(distance.x), distance.x := 0] %>%
    .[is.na(distance.y), distance.y := 0] %>%
    # Compute net distance as the mean of distances to the base and contrasting factor
    .[, c("distance", "distance.x", "distance.y") := list(Rfast::rowmeans(cbind(distance.x, distance.y)), NULL, NULL)] %>%
    # And add back categories
    merge(DATA[[taxa]]@experiment.meta[, .(rsc.ID, ee.ID, cf.Cat)], by = c("rsc.ID", "ee.ID"), all.x = T, sort = F, allow.cartesian = T) %>%
    .[cf.BaseLongUri != cf.ValLongUri] %>%
    unique() %>%
    {
      if (!POST) {
        .[, .(rsc.ID, ee.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri, reverse = F, distance)]
      } else {
        # Finally reorder
        reorderTags(.) %>%
          # And use some heuristics to make our corpus a little more lean
          .[, N := length(unique(ee.ID)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
          .[, ees := paste0(unique(ee.ID), collapse = ""), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
          # For every chain (beginning with each rsc.ID, distance == 0), only keep the minimum distance version
          # per unique ees
          .[, .SD[1], .(rsc.ID, ees)] %>%
          # .[distance == 0 | (N > 1 & N < 200 & distance < 3)] %>%
          
          .[, .(rsc.ID, ee.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri, reverse, distance)]
      }
    }
}


#' Normalize scores
#'
#' Convert raw rankings (@seealso search) to z-scores, and also calculate a column needed for experiment-wise
#' normalization.
#'
#' @param scores The output of (@seealso search).
#' @param taxa The taxon
normalize <- function(scores, #taxa = getConfig(key = "taxa")$value
                      taxa) {
  scores %>%
    merge(NULLS[[taxa]], by = "rn", sort = F) %>%
    .[, lapply(.SD, function(x) (x - score.mean) / score.sd),
      .SDcols = !c("score.mean", "score.sd", "rn", "score", "f.IN", "f.OUT", "ee.q")
    ] %>%
    {
      data.table::data.table(rn = scores$rn, ., scores[, .(score, f.IN, f.OUT, ee.q)])
    } %>%
    .[, score := rowSums2(as.matrix(.SD), na.rm = T), .SDcols = !c("rn", "score", "f.IN", "f.OUT", "ee.q")] %>%
    .[, normalization := eval(getOption("app.algorithm.experiment"))] %>%
    .[, sn := score * normalization] %>%
    data.table::setorder(-sn) %>%
    .[, !"sn"]
}


#' Enrich
#'
#' Given rankings (@seealso search), generate a ranking-weighted count of all terms that can be
#' derived from tags present in the experiment.
#'
#' @param rankings The output of (@seealso search).
#' @param options The options
#' @param doNorm Whether or not to normalize scores
#' @param CACHE A cache to use. If null, uses the global CACHE.BACKGROUND
enrich <- function(rankings, # options = getConfig(),
                   taxa,
                   dist,
                   categories,
                   doNorm = T, CACHE = NULL, cores = 8) {
  tictoc::tic()
  # terms <- getTags(options$taxa$value, rankings$rn, options$dist$value, CACHE = CACHE)
  terms <- getTags(taxa, rankings$rn, dist, CACHE = CACHE)
  terms <- rankings %>%
    {
      if (doNorm) {
        # normalize(., options$taxa$value)
        normalize(., taxa)
      } else {
        .
      }
    } %>%
    merge(terms, by.x = "rn", by.y = "rsc.ID", sort = F) %>%
    .[score > 0] %>%
    # .[cf.Cat %in% options$categories$value]
    .[cf.Cat %in% categories]
  # .[, N := length(unique(ee.ID)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
  # .[N > 1, !'N']
  
  # TODO Excluding singles because their p-value will always be 0.5
  # and so score highly driven by fIN/fOUT... This is maybe okay
  
  terms_bckp<<-terms
  
  # filtering singles
  grouping = paste(terms$cf.Cat,terms$cf.BaseLongUri,terms$cf.ValLongUri)
  valid_groups = grouping %>% table %>% {names(.[.>1])}
  terms = terms[grouping %in% valid_groups,]
  
  gene_names = 
    colnames(terms)[!colnames(terms) %in% c("reverse", "distance", "ee.ID", "rn", "score", "f.IN", "f.OUT", "ee.q", "normalization",'cf.Cat','cf.BaseLongUri','cf.ValLongUri')]
  grouping_vars = c('cf.Cat','cf.BaseLongUri','cf.ValLongUri')
  
  
  keys<- terms %>% 
    dplyr::group_by(cf.Cat,cf.BaseLongUri,cf.ValLongUri) %>% group_keys()
  term_ps <- terms %>% 
    dplyr::group_by(cf.Cat,cf.BaseLongUri,cf.ValLongUri) %>% 
    {
      .[,c(gene_names,grouping_vars)]
    } %>% group_split()  %>%  mclapply(function(x){
      out = 1-matrixTests::col_wilcoxon_onesample(as.matrix(x[gene_names]),alternative= 'greater', exact = FALSE)$pvalue
    },mc.cores = cores) %>% do.call(rbind,.)
  colnames(term_ps) = gene_names
  
  wilcox_ps = cbind(keys,term_ps)
  
  
  terms[, .(
    distance = round(mean(distance, na.rm = T), 2),
    stat = mean(score),
    normalization = median(normalization)
  ), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    merge(
      wilcox_ps,
      by = c("cf.Cat", "cf.BaseLongUri", "cf.ValLongUri"), sort = F
    ) %>%
    .[, score := Rfast::rowsums(as.matrix(.SD) * normalization), .SDcols = !c("stat", "cf.Cat", "cf.BaseLongUri", "cf.ValLongUri", "distance", "normalization")] %>%
    .[, !"normalization"] %>%
    data.table::setorder(-stat, distance) %>%
    .[, .SD[1], stat] %>%
    data.table::setorder(-score) -> out
  tictoc::toc()
  return(out)
}





#' Search
#'
#' Uses the M-VSM to sort experiments that show at least one of the genes as DE.
#'
#' @param genes A list of Entrez Gene IDs (ie. 1, 22, 480) as characters
#' @param options Optional extra parameters to pass from @seealso(getConfig)
search <- function(genes, 
                   taxa,
                   confounds,
                   filter = NULL,
                   mfx,
                   geeq,
                   p_threshold,
                   DATA = NULL) {
  if (is.null(DATA)) {
    DATA <- DATA.HOLDER
  }
  
  # mData <- DATA[[options$taxa$value]]
  mData <- DATA[[taxa]]
  
  # if (!options$confounds$value) {
  if(!confounds){
    experimentMask <- !mData@experiment.meta$ef.IsBatchConfounded
  } else {
    experimentMask <- rep(T, nrow(mData@experiment.meta))
  }
  
  # if (!is.null(options$filter$value)) {
  #   experimentMask <- experimentMask & options$filter$value
  # }
  if(!is.null(filter)){
    experimentMask <- experimentMask & filter
  }
  
  # Only retain GOI
  geneMask <- which(mData@gene.meta$entrez.ID %in% genes)
  
  n.genes <- length(geneMask)
  if (n.genes == 0) {
    return(NULL)
  }
  # the signature input appears to be disabled in the upstream menus which means
  # this piece of code isn't doing anything. there are other conflicting references in
  # the shiny code too either to "input$sig" or "input$signature" if we bring it
  # back this check should be moved into the server code to keep the function independent
  # from options - ogan
  # query <- suppressWarnings(options$sig$value %>% as.numeric())
  # print(query)
  # # TODO Should signal why it stopped (mismatch of signature length)
  # if (length(query) > 1 && length(query) != n.genes) {
  #   return(NULL)
  # }
  
  # if (length(query) == 0 && options$method$value != "diff") {
  #   return(NULL)
  # }
  
  # P-values for only the GOI
  mDimNames <- dimnames(mData@data$adj.pv)
  pv <- mData@data$adj.pv[, geneMask, drop = F] %>%
    `dimnames<-`(list(mDimNames[[1]], mDimNames[[2]][geneMask])) %>%
    .[experimentMask, , drop = F] %>%
    t()
  zScore <- mData@data$zscore[, geneMask, drop = F] %>%
    `dimnames<-`(list(mDimNames[[1]], mDimNames[[2]][geneMask])) %>%
    .[experimentMask, , drop = F] %>%
    t()
  
  pv[is.na(pv)] <- 1
  zScore[is.na(zScore)] <- 0
  
  zScore <- data.table::as.data.table(zScore)
  
  # Number of DEGs for experiments that pass thresholds
  experimentMeta <- mData@experiment.meta[experimentMask, .(rsc.ID, ee.qScore, n.DE, ad.NumGenes)] %>%
    as.data.frame() %>%
    `rownames<-`(.[, "rsc.ID"])
  
  experimentMeta$n.DE[is.na(experimentMeta$n.DE)] <- 0
  
  if (!geeq) {
    experimentMeta$ee.qScore <- 2
  } else {
    experimentMeta$ee.qScore <- pmax(experimentMeta$ee.qScore + 1, 0, na.rm = T)
  }
  experimentMeta$ee.qScore <- sqrt(experimentMeta$ee.qScore / 2)
  if (mfx) {
    MFX_WEIGHT <- 1.5 - mData@gene.meta$mfx.Rank[geneMask]
  } else {
    MFX_WEIGHT <- rep(1, n.genes)
  }
  
  # Adding 1 to not exclude CCs where p.adj is 1
  zScore <- getOption("app.algorithm.gene.pre") %>% eval()
  
  ret <- getOption("app.algorithm.gene.post") %>%
    eval() %>%
    .[, score := Rfast::rowsums(as.matrix(.))]
  
  experimentN <- Rfast::colsums(pv <= p_threshold)
  
  ret %>%
    data.table::setnames(c(mData@gene.meta$gene.Name[geneMask], "score")) %>%
    .[, f.IN := experimentN / n.genes] %>%
    .[, f.OUT := pmax(0, experimentMeta$n.DE - experimentN) / experimentMeta$ad.NumGenes] %>%
    .[, ee.q := experimentMeta$ee.qScore] %>%
    .[, rn := colnames(zScore)] %>%
    data.table::setorder(-score)
}


