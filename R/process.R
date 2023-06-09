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
#' @param filter Should the tags be filtered for recorded stop words
getTags <- function(taxa,
                    rsc.IDs = NULL, max.distance = Inf, inv = F, CACHE = NULL,filter = FALSE) {
  if (length(taxa) > 1) {
    return(data.table::rbindlist(lapply(taxa, getTags, rsc.IDs, max.distance, inv, CACHE,filter)))
  }
  
  if (is.null(CACHE)) {
    CACHE <- CACHE.BACKGROUND
  }
  
  if(filter){
    CACHE <- CACHE %>% lapply(function(x){
      x %>% dplyr::filter(!(cf.ValLongUri %in% filters$universal_filter | cf.BaseLongUri %in% filters$universal_filter)) %>%
        dplyr::filter(!cf.ValLongUri %in% filters$val_filter) %>%
        dplyr::filter(!cf.BaseLongUri %in% filters$base_filter)
    })
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
  
  # if(is.null(ONTOLOGIES.DEFS)){
  #   ONTOLOGIES.DEFS = parent.frame()$ONTOLOGIES.DEFS
  # }
  
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
  
  # filter mComputable to remove specific terms from propagating
  
  
  mComputedTags <- data.table::rbindlist(parallel::mclapply(mComputable, function(uri) {
    tag <- igraph::subcomponent(mGraph, uri, "out")
    distance <- igraph::distances(mGraph, uri, tag)
    data.table(startTag = uri, tag = names(tag), distance = c(distance))
  },mc.cores = 6))
  # saveRDS(mComputedTags,paste0(taxa,"computedTags.rds"))
  
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
  
  # this part used to replace URIs with plain text names. it is moved to the end
  # to allow retaining both URIs and plain text names in the output
  # mTags <- mTags %>%
  #   merge(unique(ONTOLOGIES.DEFS[, .(Node_Long = as.character(Node_Long), Definition = as.character(Definition))]),
  #         by.x = "tag", by.y = "Node_Long", sort = F, allow.cartesian = T, all.x = T
  #   ) %>%
  #   .[is.na(Definition), Definition := tag] %>%
  #   .[, .(rsc.ID, ee.ID, type, tag = Definition, distance, ID, uri = tag)] %>%
  #   unique()
  
  
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
          by.x = c("rsc.ID", "ee.ID", "cf.BaseLongUri"), by.y = c("rsc.ID", "ee.ID","tag"), all.x = T, sort = F, allow.cartesian = T
    ) %>%
    merge(mTags[type == "cf.ValLongUri", .(rsc.ID, ee.ID, tag, distance)],
          by.x = c("rsc.ID", "ee.ID", "cf.ValLongUri"), by.y = c("rsc.ID", "ee.ID","tag"), all.x = T, sort = F, allow.cartesian = T
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
    } -> out
  
  
  # this is a bit of wrangling to preserve the URIs. The terminology used in the column names
  # is a bit misleading due to historical reasons. The cf.ValLongUri and cf.BaseLongUri always
  # included the plain text definitions rather than the ontology uris. This piece
  # adds the URIs to the output - Ogan
  # out %>%  merge(unique(ONTOLOGIES.DEFS[, .(Node_Long = as.character(Node_Long), cf.BaseOriginal = as.character(Definition))]),
  #                by.x = "cf.BaseLongUri", by.y = "Node_Long", sort = F, allow.cartesian = T, all.x = T
  # ) %>% 
  #   merge(unique(ONTOLOGIES.DEFS[, .(Node_Long = as.character(Node_Long), cf.ValOriginal = as.character(Definition))]),
  #        by.x = "cf.ValLongUri", by.y = "Node_Long", sort = F, allow.cartesian = T, all.x = T
  #   ) %>% 
  #   .[,.(cf.ValLongUri = cf.ValOriginal, cf.BaseLongUri = cf.BaseOriginal, rsc.ID, ee.ID,cf.Cat, reverse, distance, cf.BaseOriginal = cf.BaseLongUri,cf.ValOriginal = cf.ValLongUri)]
  return(out)
}


#' Normalize scores
#'
#' Convert raw rankings (@seealso search) to z-scores, and also calculate a column needed for experiment-wise
#' normalization.
#'
#' @param scores The output of (@seealso search).
#' @param taxa The taxon
normalize <- function(scores, #taxa = getConfig(key = "taxa")$value
                      taxa,nullset) {
  if(is.null(nullset)){
    nullset = NULLS
  }
  # using the merged object allows some experiments to be missing in the nullset
  # this only happens if the nullset is created with a blocksize which causes
  # null scores for some experiments (14 contrasts in total) with low gene coverage to fail
  merged <- scores %>%
    merge(nullset[[taxa]], by = "rn", sort = F)
  merged %>%
    .[, lapply(.SD, function(x) (x - score.mean) / score.sd),
      .SDcols = !c("score.mean", "score.sd", "rn", "score", "f.IN", "f.OUT", "ee.q")
    ] %>%
    {
      data.table::data.table(rn = merged$rn, ., merged[, .(score, f.IN, f.OUT, ee.q)])
    } %>%
    .[, score := matrixStats::rowSums2(as.matrix(.SD), na.rm = T), .SDcols = !c("rn", "score", "f.IN", "f.OUT", "ee.q")] %>%
    .[, normalization := eval(getOption("app.algorithm.experiment"))] %>%
    .[, sn := score * normalization] %>%
    data.table::setorder(-sn) %>% # order by normalized scores then remove them?
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
                   filter = TRUE,
                   doNorm = T, CACHE = NULL, cores = 8,nullset = NULL) {
  if(is.null(nullset)){
    nullset = NULLS
  }
  # terms <- getTags(options$taxa$value, rankings$rn, options$dist$value, CACHE = CACHE)
  terms <- getTags(taxa, rankings$rn, dist, CACHE = CACHE, filter = filter)
  terms <- rankings %>%
    {
      if (doNorm) {
        # normalize(., options$taxa$value)
        normalize(., taxa,nullset)
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
  grouping = grouping[grouping %in% valid_groups]
  
  
  # assign rows belonging to same groups to same cores
  core_split = grouping %>% factor %>% as.numeric() %>% {.%%cores}
  gene_names = 
    colnames(terms)[!colnames(terms) %in% c("reverse", "distance", "ee.ID", "rn", "score", "f.IN", "f.OUT", "ee.q", "normalization",'cf.Cat','cf.BaseLongUri','cf.ValLongUri')]
  grouping_vars = c('cf.Cat','cf.BaseLongUri','cf.ValLongUri')
  
  terms[,core_split := core_split]
  grouped = terms %>% data.table:::split.data.table(by = c('core_split'))
  
  grouped %>% parallel::mclapply(function(t){
    t[, matrixTests::col_wilcoxon_onesample(as.matrix(.SD), alternative = "greater", exact = F) %>%
        {
          list(pv = 1 - .[, "pvalue"], gene = rownames(.)) # Small p-value if real effect (= score closer to 1)
        }, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri),
      .SDcols = !c("reverse", "distance", "ee.ID", "rn", "score", "f.IN", "f.OUT", "ee.q", "normalization","core_split")
    ] %>%
      data.table::dcast(... ~ gene, value.var = "pv", fill = 0)
  },mc.cores = cores) %>% do.call(rbind,.) -> wilcox_ps
  
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
  enrich_out <<- out
  return(out)
}





#' Search
#'
#' Uses the M-VSM to sort experiments that show at least one of the genes as DE.
#'
#' @param genes A list of Entrez Gene IDs (ie. 1, 22, 480) as characters
#' @param options Optional extra parameters to pass from @seealso(getConfig)
vsmSearch <- function(genes, 
                   taxa,
                   confounds,
                   diff_exp_source = 'zscore',
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
    experimentMask <- !mData@experiment.meta$ef.IsBatchConfounded | is.na(mData@experiment.meta$ef.IsBatchConfounded)
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
  if('FBM' %in% class(mData@data[[diff_exp_source]])){
    zScore <- mData@data[[diff_exp_source]][, geneMask, drop = FALSE] %>%
      `dimnames<-`(list(mDimNames[[1]], mDimNames[[2]][geneMask])) %>%
      .[experimentMask, , drop = F] %>%
      t()
  } else if('matrix' %in% class(mData@data[[diff_exp_source]])){
    zScore <- mData@data[[diff_exp_source]][geneMask,,drop = FALSE] %>%
      `dimnames<-`(list(mDimNames[[2]][geneMask],mDimNames[[1]])) %>%
      .[,experimentMask , drop = F]
  }
  
  
  # remove experiments where none of the query genes are represented
  all_nas = apply(pv,2,function(x){all(is.na(x))})
  pv = pv[,!all_nas,drop = FALSE]
  zScore = zScore[,!all_nas, drop = FALSE]
  pv[is.na(pv)] <- 1
  zScore[is.na(zScore)] <- 0
  
  zScore <- data.table::as.data.table(zScore)
  # Number of DEGs for experiments that pass thresholds
  experimentMeta <- mData@experiment.meta[experimentMask, .(rsc.ID, ee.qScore, n.DE, ad.NumGenes,nonNa.numGenes)] %>%
    as.data.frame() %>%
    `rownames<-`(.[, "rsc.ID"])
  

  experimentMeta  = experimentMeta[!all_nas,]
  
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
  ret<<- ret
  experimentN <- Rfast::colsums(pv <= p_threshold)
  ret %>%
    data.table::setnames(c(mData@gene.meta$gene.Name[geneMask], "score")) %>%
    .[, f.IN := experimentN / n.genes] %>%
    .[, f.OUT := pmax(0, experimentMeta$n.DE - experimentN) / experimentMeta$nonNa.numGenes] %>%
    .[, ee.q := experimentMeta$ee.qScore] %>%
    .[, rn := colnames(zScore)] %>%
    data.table::setorder(-score)
}


#' Process input taxa
#' @param taxa taxonomy ids, common or scientific names'
#' @return A vector with the taxanomic IDs of the input taxa
processTaxa <- function(taxa){
  apply(TAX.DATA, 1, function(x){
    (taxa %in% x)
  }) %>% matrix(ncol = ncol(TAX.DATA)) %>% apply(2,any) %>%
    {out=TAX.DATA$id[.];names(out)=TAX.DATA$common_names[.];out}
}

#' Process input genes
#' 
#' Replacement for the tidyGenes function for the API. 
#' Only difference is that it always returns the data.table
#' 
#' @param genes a vector of genes. Symbols, ensembl ids or 
#' @param taxa Taxanomy ids
processGenes <- function(genes, taxa){
  taxIDs <- processTaxa(taxa)
  
  symbols = taxIDs %>%
    lapply(function(x){
      tax_name <- TAX.DATA$common_names[TAX.DATA$id == x]
      
      symbols <- data.table(entrez.ID = as.character(tidyGenes(genes, tax_name)$genes), taxon = tax_name) %>%
        {
          if (is.null(.) || ncol(.) == 1) {
            NULL
          } else {
            merge(., DATA.HOLDER[[tax_name]]@gene.meta[, .(.I, entrez.ID = as.character(entrez.ID), gene.Name)], by = "entrez.ID")
          }
        }
    }) %>% data.table::rbindlist()
  
  if (!is.null(symbols) && nrow(symbols) > 0) {
    symbols <- symbols %>%
      merge(data.table(taxon.ID = taxIDs, taxon = names(taxIDs)), by = "taxon") %>%
      .[, .SD[1], gene.Name]
  } else {
    return(symbols)
  }
  
  
  # Identify homologs/orthologs/whatever
  orthologs <- symbols[, list(lapply(taxIDs[taxIDs != unique(taxon.ID)], function(t) {
    homologs <- homologene::homologene(entrez.ID, unique(taxon.ID), t)
    if (nrow(homologs) == 0) {
      NULL
    } else {
      homologs %>%
        data.table::as.data.table() %>%
        data.table::melt(measure.vars = which(!grepl("_ID", colnames(.))), value.name = "gene.realName", variable.name = "taxon.ID") %>%
        .[, entrez.ID := unlist(lapply(1:nrow(.), function(I) {
          .SD[I, grepl(paste0(taxon.ID[I], "_ID"), colnames(.SD)), with = F]
        }))] %>%
        .[, !grepl("_ID", colnames(.)), with = F] %>%
        .[, taxon.ID := as.integer(levels(taxon.ID))[taxon.ID]] %>%
        cbind(identifier = homologs[, as.character(unique(taxon.ID))]) %>%
        merge(DATA.HOLDER[[names(taxIDs)[taxIDs == t]]]@gene.meta[, .(.I, entrez.ID = as.integer(entrez.ID))], by = "entrez.ID")
    }
  })), taxon] %>%
    .[, V1] %>%
    data.table::rbindlist() %>%
    rbind(symbols[, .(entrez.ID, taxon.ID, gene.realName = gene.Name, identifier = gene.Name, I)]) %>%
    unique() %>%
    merge(data.table(taxon.ID = taxIDs, taxon = names(taxIDs)), by = "taxon.ID") %>%
    .[, entrez.ID := as.character(entrez.ID)] %>%
    .[, identifier := make.names(identifier, T), taxon.ID]
  
  # subsetting here fixes inputting same gene orthologue for different genes in one query
  return(orthologs[!duplicated(I),.SD])  
}


#' Convert gene identifiers/symbols/names/etc to Entrez IDs and generate suggestions
#' for any mismatches or homologs if taxa is plural
#'
#' @param genes A set of genes in any supported format
#' @param taxa A (set of) taxa these genes may belong to
tidyGenes <- function(genes, taxa) {
  if (length(taxa) > 1) {
    taxIDs <- processTaxa(taxa)
    
    symbols <- lapply(taxa, function(x) {
      data.table(entrez.ID = tidyGenes(genes, x)$genes, taxon = x) %>%
        {
          if (is.null(.) || ncol(.) == 1) {
            NULL
          } else {
            merge(., DATA.HOLDER[[x]]@gene.meta[, .(.I, entrez.ID = as.character(entrez.ID), gene.Name)], by = "entrez.ID")
          }
        }
    }) %>% data.table::rbindlist()
    
    if (!is.null(symbols) && nrow(symbols) > 0) {
      symbols <- symbols %>%
        merge(data.table(taxon.ID = taxIDs, taxon = names(taxIDs)), by = "taxon") %>%
        .[, .SD[1], gene.Name]
    } else {
      return(tidyGenes(genes, taxa[1]))
    }
    
    # Identify homologs/orthologs/whatever
    orthologs <- symbols[, list(lapply(taxIDs[taxIDs != unique(taxon.ID)], function(t) {
      homologs <- homologene::homologene(entrez.ID, unique(taxon.ID), t)
      if (nrow(homologs) == 0) {
        NULL
      } else {
        homologs %>%
          data.table::as.data.table() %>%
          data.table::melt(measure.vars = which(!grepl("_ID", colnames(.))), value.name = "gene.realName", variable.name = "taxon.ID") %>%
          .[, entrez.ID := unlist(lapply(1:nrow(.), function(I) {
            .SD[I, grepl(paste0(taxon.ID[I], "_ID"), colnames(.SD)), with = F]
          }))] %>%
          .[, !grepl("_ID", colnames(.)), with = F] %>%
          .[, taxon.ID := as.integer(levels(taxon.ID))[taxon.ID]] %>%
          cbind(identifier = homologs[, as.character(unique(taxon.ID))]) %>%
          merge(DATA.HOLDER[[names(taxIDs)[taxIDs == t]]]@gene.meta[, .(.I, entrez.ID = as.integer(entrez.ID))], by = "entrez.ID")
      }
    })), taxon] %>%
      .[, V1] %>%
      data.table::rbindlist() %>%
      rbind(symbols[, .(entrez.ID, taxon.ID, gene.realName = gene.Name, identifier = gene.Name, I)]) %>%
      unique() %>%
      merge(data.table(taxon.ID = taxIDs, taxon = names(taxIDs)), by = "taxon.ID") %>%
      .[, entrez.ID := as.character(entrez.ID)] %>%
      .[, identifier := make.names(identifier, T), taxon.ID]
    # subsetting here fixes inputting same gene orthologue for different genes in one query
    return(orthologs[!duplicated(I),.SD])
  }
  
  oGenes <- genes
  
  # Clean numerics (interpreted as entrez IDs) and remove them from further processing.
  cleanGenes <- suppressWarnings(Filter(function(x) !is.na(as.integer(x)), genes))
  idMap <- sapply(as.character(cleanGenes), function(x) which(oGenes == x), USE.NAMES = F)
  
  genes <- genes[!(genes %in% cleanGenes)]
  
  # If it matches (ENSG|ENSMUS|ENSRNO)\d{11}, it's an Ensembl ID (for human, mouse or rat).
  if (length(genes) > 0) {
    ensembl <- grep("(ENSG|ENSMUS|ENSRNO)\\d{11}", genes, value = T)
    
    if (length(ensembl) != 0) {
      # Extract genes with a matching Ensembl ID and clean them too.
      ensembls <- DATA.HOLDER[[taxa]]@gene.meta[ensembl.ID %in% ensembl, .(entrez.ID, ensembl.ID)]
      cleanGenes <- c(cleanGenes, ensembls[, entrez.ID])
      idMap <- c(idMap, sapply(ensembls[, ensembl.ID], function(x) which(oGenes == x), USE.NAMES = F))
      
      genes <- genes[!(genes %in% ensembls[, ensembl.ID])]
    }
  }
  
  # Match GO identifiers
  if (length(genes > 0)) {
    go <- grep("(GO:)\\d{7}", genes, value = T)
    
    if (length(go) != 0) {
      gos <- DATA.HOLDER[[taxa]]@go[id %in% go, .(id, entrez.ID)]
      cleanGenes <- c(cleanGenes, gos[, entrez.ID])
      idMap <- c(idMap, sapply(gos[, id], function(x) which(oGenes == x), USE.NAMES = F))
      
      genes <- genes[!(genes %in% gos[, id])]
    }
  }
  
  # Try to match to gene names
  if (length(genes) > 0) {
    descriptors <- DATA.HOLDER[[taxa]]@gene.meta[gene.Name %in% genes, .(entrez.ID, gene.Name)]
    if (nrow(descriptors) != 0) {
      cleanGenes <- c(cleanGenes, descriptors[, entrez.ID])
      idMap <- c(idMap, sapply(descriptors[, gene.Name], function(x) which(oGenes == x), USE.NAMES = F))
      
      genes <- genes[!(genes %in% descriptors[, gene.Name])]
    }
  }
  
  # Make suggestions for approximate matches
  mSuggestions <- NULL
  if (length(genes) > 0) {
    mSuggestions <- sapply(genes, function(gene) {
      stringdist::stringdist(tolower(gene), tolower(DATA.HOLDER[[taxa]]@gene.meta[, gene.Name])) %>%
        {
          DATA.HOLDER[[taxa]]@gene.meta[which.min(.), gene.Name]
        }
    }) %>% stats::setNames(genes)
  }
  
  # Return the new genes in the same order they were provided (hopefully)
  mGenes <- NULL
  if (length(cleanGenes) > 0) {
    mGenes <- cleanGenes[order(unlist(idMap))]
  }
  
  if (length(mGenes) > 0 || length(mSuggestions) > 0) {
    list(genes = mGenes, suggestions = mSuggestions)
  } else {
    NULL
  }
}


# given a cache, return the list of possible results
get_possible_results = function(cache=NULL, filter = TRUE){
  if(is.null(cache)){
    cache = CACHE.BACKGROUND
  }
  
  possible_results = cache %>% 
    lapply(function(x){
      x <- x %>% 
        dplyr::group_by(cf.Cat, cf.BaseLongUri, cf.ValLongUri) %>% 
        dplyr::summarise(n = dplyr::n()) %>% 
        merge(unique(SIMPLIFIED.ONTOLOGY.DEFS[, .(Node_Long = as.character(Node_Long), cf.Val = as.character(Definition))]),
              by.x = "cf.ValLongUri",
              by.y = "Node_Long",
              sort = F, 
              allow.cartesian = T, 
              all.x = T
        ) %>%
        merge(unique(SIMPLIFIED.ONTOLOGY.DEFS[, .(Node_Long = as.character(Node_Long), cf.Base = as.character(Definition))]),
              by.x = "cf.BaseLongUri",
              by.y = "Node_Long", 
              sort = F,
              allow.cartesian = T, 
              all.x = T) %>%  dplyr::arrange(dplyr::desc(n)) %>%
        dplyr::select(cf.Cat,n,cf.ValLongUri,cf.BaseLongUri,cf.Val,cf.Base)
      
      # x$cfBase[is.na(x$cf.Base)] <- x$cf.BaseLongUri[is.na(x$cf.Base)]
      # x$cf.Val[is.na(x$cf.Val)] <- x$cf.ValLongUri[is.na(x$cf.Val)]
      
      return(x)
    })
  
  if(filter){
    possible_results <- possible_results %>% lapply(function(x){
      x %>% dplyr::filter(!(cf.ValLongUri %in% filters$universal_filter | cf.BaseLongUri %in% filters$universal_filter)) %>%
        dplyr::filter(!cf.ValLongUri %in% filters$val_filter) %>%
        dplyr::filter(!cf.BaseLongUri %in% filters$base_filter)
    })
  }
  
  
   return(possible_results)
}


# unified function to run the whole test
de_search = function(genes = NULL,
                     taxa =NULL,
                     diff_exp_source  = 'zscore',
                     max_dist = 1.5,
                     confounds = FALSE,
                     multifunctionality = FALSE,
                     geeq = FALSE,
                     p_threshold = 0.05,
                     categories = c("age", "behavior", "biological process", "biological sex", 
                                    "cell type", "clinical history", "diet", "disease", "environmental history", 
                                    "environmental stress", "genotype", "medical procedure", "molecular entity", 
                                    "organism part", "phenotype", "sex", "temperature", "treatment"),
                     filter_stopwords = TRUE,
                     remove_experiments = NULL,
                     remove_comparisons = NULL,
                     cache = NULL,
                     get_descriptions = TRUE, # temporary argument to allow supporting old cache files
                     cores = 8,
                     nullset = NULL){
  if(is.null(nullset)){
    nullset = NULLS
  }

  if(!is.logical(geeq)){
    geeq = as.logical(toupper(geeq))
  }
  if(!is.logical(multifunctionality)){
    multifunctionality = as.logical(toupper(multifunctionality))
  }
  if(!is.logical(confounds)){
    confounds = as.logical(toupper(confounds))
  }
  
  
  # tictoc::tic()
  genes <- processGenes(genes,taxa)
  # print('vsmSearch')
  # tictoc::tic()
  experiments <- taxa %>% 
    parallel::mclapply(function(t){
      
      exp_filter = DATA.HOLDER[[t]]@experiment.meta$ee.Name %in% remove_experiments
      comp_filter = DATA.HOLDER[[t]]@experiment.meta$rsc.ID %in% remove_comparisons
      
      
      vsmSearch(genes[taxon == t, entrez.ID],
                taxa = t,
                confounds = confounds,
                diff_exp_source = diff_exp_source,
                filter = !(exp_filter | comp_filter),
                mfx = multifunctionality,
                geeq = geeq,
                p_threshold = p_threshold)
    },mc.cores = cores)
  names(experiments) = taxa
  # tictoc::toc()
  
  # print('enrich')
  # tictoc::tic()
  conditions <- taxa %>% lapply(function(t){
    enrich(experiments[[t]], taxa = t, dist = max_dist,categories = categories,cores = cores,filter = filter_stopwords,CACHE = cache,nullset = nullset) %>% 
      data.table::setnames(genes[taxon == t, gene.realName],
                           genes[taxon == t, identifier],
                           skip_absent = T)
  }) %>% data.table::rbindlist(fill = TRUE) %>% 
    # reorderTags3() %>% # appears to be redundant. cache tags are already re-ordered
    .[, lapply(.SD, mean, na.rm = T), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
  # tictoc::toc()
  
  # print('ending')
  # tictoc::tic()
  geneInfo <- genes %>%
    data.table::copy() %>%
    data.table::setnames("identifier", "gene.Name")
  mGenes <- genes %>% data.table::copy()
  
  # components of the endSuccess function
  exps <- lapply(experiments, "[[", "rn") %>% unlist()
  
  tmp <- data.table::rbindlist(lapply(taxa, function(i) {
    DATA.HOLDER[[i]]@experiment.meta[rsc.ID %in% exps, .(rsc.ID, ee.ID, ee.Name, ee.NumSample, ef.IsBatchConfounded)]
  }))
  

  if (get_descriptions){
    conditions %<>%  
      merge(unique(SIMPLIFIED.ONTOLOGY.DEFS[, .(Node_Long = as.character(Node_Long), cf.Base = as.character(Definition))]),
            by.x = "cf.BaseLongUri",
            by.y = "Node_Long", 
            sort = F,
            allow.cartesian = T, 
            all.x = T) %>% 
      merge(unique(SIMPLIFIED.ONTOLOGY.DEFS[, .(Node_Long = as.character(Node_Long), cf.Val = as.character(Definition))]),
            by.x = "cf.ValLongUri",
            by.y = "Node_Long",
            sort = F, 
            allow.cartesian = T, 
            all.x = T
      )
    
    # if a field is filled with free text, the processing code moves them to 
    # URIs as is. here we fill those blank spaces with the free texts
    # again so we don't have blank results
    # pending a discussion with Paul, and Neera these might be removed entirely
    # in the future
    conditions$cf.Base[is.na(conditions$cf.Base)] <- conditions$cf.BaseLongUri[is.na(conditions$cf.Base)]
    conditions$cf.Val[is.na(conditions$cf.Val)] <- conditions$cf.ValLongUri[is.na(conditions$cf.Val)]
    
    
    conditions[, `Condition Comparison` := paste0( cf.Base, " vs. ", cf.Val)]
    
    conditions %>% data.table::setcolorder(c('cf.Base','cf.Val'))
    
  } else{
    conditions[, `Condition Comparison` := paste0( cf.BaseLongUri, " vs. ", cf.ValLongUri)]
  }
  
  
  tmp <- tmp %>%
    merge(getTags(taxa, exps), sort = F) %>%
    .[, N := length(unique(ee.ID)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    .[N < 0.03 * nrow(tmp)] # Get rid of contrasts that overlap in more than 3% experiments
  
  tmp[, Evidence := {
    dedup <- !duplicated(ee.ID)
    stringi::stri_c(ee.Name[dedup], collapse = ",")
  },.(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
  
  
  tmp[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, N, Evidence)] %>%
    unique() %>%
    merge(conditions, by = c("cf.Cat", "cf.BaseLongUri", "cf.ValLongUri"), sort = F) %>%
    data.table::setnames(c("stat", "score", "distance"), c("Effect Size", "Test Statistic", "Ontology Steps")) ->
    conditions
  
  # out of endSuccess
  
  getPercentageStat <- function(x, n = 1){
    x / n
  }
  conditions[,'Test Statistic'] <- apply(conditions[,'Test Statistic'], 2, getPercentageStat, n = nrow(geneInfo))
  # tictoc::toc()
  
  # tictoc::toc()
  return(conditions %>% 
           data.table::setcolorder(c('Condition Comparison',"cf.Cat", 'Test Statistic')) %>% 
           data.table::setnames(c("cf.Cat", "cf.BaseLongUri", "cf.ValLongUri"), c("Category", "Baseline", "Value")))
}

contribs = function(data){
  mData <- data %>%
    setorder(-`Test Statistic`) %>%
    .[, !c(
      "Test Statistic", "Effect Size", "Ontology Steps", "N", "Evidence","EvidencePlain",
      "cf.Cat", "cf.BaseLongUri", "cf.ValLongUri"
    )]
}


get_parents = function(terms){
  mGraph <- igraph::simplify(igraph::graph_from_data_frame(ONTOLOGIES[, .(as.character(ChildNode_Long), as.character(ParentNode_Long))]))
  terms %>% lapply(function(x){
    tryCatch( igraph::subcomponent(mGraph, x, "out") %>% names, error = function(e){
      return(x)
    })
  }) %>% unlist %>% unique
}
