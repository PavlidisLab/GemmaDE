dopeScores <- function(taxa,
                       genes,
                       geneIDs = NULL,
                       examples,
                       N_conditions = 1,
                       mConditions = NULL,
                       regularizer = 5, # Higher numbers make it more linear. Negative numbers make it reverse. 0 is linear
                       seed = NULL) {
  if(!is.null(seed))
    set.seed(seed)
  
  # Select genes that we'll work with
  nGenes <- nrow(DATA.HOLDER[[taxa]]@gene.meta)
  
  if(!is.null(geneIDs)) {
    nSelected <- which(DATA.HOLDER[[taxa]]@gene.meta$entrez.ID %in% geneIDs)
    genes <- length(nSelected)
  } else
    nSelected <- sample(1:nGenes, max(genes))
  
  mGenes <- DATA.HOLDER[[taxa]]@gene.meta[nSelected]
  mGenes.dist <- DATA.HOLDER[[taxa]]@data$meanval[nSelected, ] %>% as.matrix %>% {
    tmp <- .
    # TODO Let's ignore this detail for now
    #tmp <- switch(mData@experiment.meta$ee.Scale,
    #               LOG2 = .,
    #               LOG10 = . / log10(2),
    #               LN = . / ln(2),
    #               log2(. + 1)
    #)
    tmp[is.nan(tmp) | is.na(tmp)] <- 0
    list(mean = rowMeans2(tmp, na.rm = T),
         sd = Rfast::rowVars(tmp, na.rm = T, std = T))
  }
  
  # Select our conditions
  if(is.null(mConditions)) {
    mConditions <- DATA.HOLDER[[taxa]]@experiment.meta[, .(cf.Cat, cf.CatLongUri, cf.Val,
                                                           cf.ValLongUri, cf.Baseline, cf.BaseLongUri)] %>%
      unique %>% .[sample(1:.N, N_conditions)]
  }
  
  # Search ahead of time for N genes
  message('Searching...')
  searched <- lapply(genes, function(i) {
    search(mGenes$entrez.ID[sample(1:nrow(mGenes), i)], taxa, verbose = F)
  })
  
  list(data = lapply(examples, function(N_examples) {
    message(paste0('Simulating ', N_examples, ' examples...'))
    lapply(genes, function(N_genes) {
      message(paste0('Simulating ', N_genes, ' genes...'))
      
      mSelected <- nSelected[1:N_genes]
      
      CACHE.BACKGROUND$dope <<- NULL
      TAGS$dope <<- NULL
      
      mName <- taxa
      searched2 <- searched
      if(N_examples > 0) {
        mData <- DATA.HOLDER[[taxa]]
        mData@taxon <- 'dope'
        mData@go <- data.table()
        
        experimentNames <- paste0('RSCID.DOPE.', 1:N_examples)
        
        mMeta <- data.table(ee.ID = mData@experiment.meta[, max(ee.ID) + 1] %>% seq(., length.out = N_examples),
                            ee.qScore = 1,
                            ee.sScore = 1,
                            rsc.ID = experimentNames,
                            ee.Troubled = F,
                            ee.Public = T,
                            ee.Name = paste0('DOPE', 1:N_examples),
                            ee.Source = 'DOPE',
                            ee.Scale = 'LOG2',
                            ee.NumSamples = 1, # TODO These aren't used in search/enrich so we don't need them (yet)
                            ee.TagLongUri = NA,
                            ad.Name = 'DOPE',
                            ad.Company = 'DOPE',
                            ad.Sequencing = T,
                            sf.Subset = F,
                            sf.Cat = NA,
                            sf.CatLongUri = NA,
                            sf.ValLongUri = NA,
                            mConditions[sample(1:N_conditions, N_examples, T)],
                            n.DE = 1,
                            mean.fc = 1)
        
        message('Precomputing cache...')
        DATA.HOLDER$dope <<- mData
        DATA.HOLDER$dope@experiment.meta <<- mMeta
        CACHE.BACKGROUND$dope <<- precomputeTags('dope') %>% reorderTags
        TAGS$dope <<- getTags('dope')
        
        DATA.HOLDER$dope <<- mData
        CACHE.BACKGROUND$dope <<- rbind(CACHE.BACKGROUND[[taxa]], CACHE.BACKGROUND$dope) %>% reorderTags
        TAGS$dope <<- rbind(TAGS[[taxa]], TAGS$dope)
        mName <- 'dope'
        
        message('Spoofing scores...')
        searched2 <- lapply(searched, function(search) {
          mSelect <- sample(1:nrow(search), N_examples, N_examples > nrow(search),
                            prob = exp(-1/regularizer * 1:nrow(search)) %>% `/`(sum(.)))
          message(paste0('Regularizer: ', regularizer, ' inserting at: ', paste0(mSelect, collapse = ', '), ' of ', nrow(search)))
          rbind(search, copy(search)[mSelect] %>% .[, rn := experimentNames]) %>% setorder(-score)
        })
      }
      
      message(paste0('Enriching ', length(searched2), ' times...'))
      enriched <- lapply(searched2, function(i) {
        if(is.null(i) || nrow(i) == 0)
          data.table()
        else
          enrich(i, mName, verbose = F)
      })
      
      list(searched = searched2,
           enriched = enriched,
           N_examples = N_examples,
           N_searched = genes,
           N_genes = N_genes)
    })
  }),
  genes = mGenes,
  conditions = mConditions,
  regularizer = regularizer,
  N_conditions = N_conditions)
}

dope <- function(taxa,
                 genes,
                 geneIDs = NULL, # Specify how many genes or which gene IDs
                 examples,
                 N_conditions = 1,
                 mConditions = NULL, # Specify how many conditions or which conditions
                 fc = list(mean = 2, sd = 0.2),
                 pv = list(mean = 0.001, sd = 0.05),
                 dropout = list(severity = 0.2, frequency = 0.1),
                 geeq = list(mean = 1, sd = 0),
                 gene_perms = 1, min_search = genes, seed = NULL) {
  if(!is.null(seed))
    set.seed(seed)
  
  mData <- DATA.HOLDER[[taxa]]
  mData@taxon <- 'dope'
  mData@go <- data.table()
  
  # Select genes that we'll work with
  nGenes <- nrow(mData@gene.meta)
  
  if(!is.null(geneIDs)) {
    nSelected <- which(mData@gene.meta$entrez.ID %in% geneIDs)
    genes <- length(nSelected)
  } else
    nSelected <- sample(1:nGenes, max(genes))
  
  mGenes <- mData@gene.meta[nSelected]
  mGenes.dist <- mData@data$meanval[nSelected, ] %>% as.matrix %>% {
    tmp <- .
    # TODO Let's ignore this detail for now
    #tmp <- switch(mData@experiment.meta$ee.Scale,
    #               LOG2 = .,
    #               LOG10 = . / log10(2),
    #               LN = . / ln(2),
    #               log2(. + 1)
    #)
    tmp[is.nan(tmp) | is.na(tmp)] <- 0
    list(mean = rowMeans2(tmp, na.rm = T),
         sd = Rfast::rowVars(tmp, na.rm = T, std = T))
  }
  
  # Select our conditions
  if(is.null(mConditions))
    mConditions <- mData@experiment.meta[, .(cf.Cat, cf.CatLongUri, cf.Val, cf.ValLongUri, cf.Baseline, cf.BaseLongUri)] %>%
    unique %>% .[sample(1:.N, N_conditions)]
  
  list(data = lapply(examples, function(N_examples) {
    message(paste0('Simulating ', N_examples, ' examples...'))
    lapply(genes, function(N_genes) {
      message(paste0('Simulating ', N_genes, ' genes...'))
      
      mSelected <- nSelected[1:N_genes]
      
      DATA.HOLDER$dope <<- NULL
      CACHE.BACKGROUND$dope <<- NULL
      TAGS$dope <<- NULL
      
      mName <- taxa
      if(N_examples > 0) {
        mData <- DATA.HOLDER[[taxa]]
        mData@taxon <- 'dope'
        mData@go <- data.table()
        
        experimentNames <- paste0('RSCID.DOPE.', 1:N_examples)
        
        # Initialize as random columns
        mFC <- mData@data$fc[, sample(1:ncol(mData@data$fc), N_examples, T)] %>% { if(N_examples == 1) t(t(.)) else . } %>% `colnames<-`(experimentNames)
        mMean <- mData@data$meanval[, sample(1:ncol(mData@data$meanval), N_examples, T)] %>% { if(N_examples == 1) t(t(.)) else . } %>% `colnames<-`(experimentNames)
        mPV <- mData@data$adj.pv[, sample(1:ncol(mData@data$adj.pv), N_examples, T)] %>% { if(N_examples == 1) t(t(.)) else . } %>% `colnames<-`(experimentNames)
        mZ <- mData@data$zscore[, sample(1:ncol(mData@data$zscore), N_examples, T)] %>% { if(N_examples == 1) t(t(.)) else . } %>% `colnames<-`(experimentNames)
        mPVZ <- mData@data$pvz[, sample(1:ncol(mData@data$pvz), N_examples, T)] %>% { if(N_examples == 1) t(t(.)) else . } %>% `colnames<-`(experimentNames)
        
        # Update with random values for the genes we care about
        mFC[mSelected, ] <- matrix(rowQuantiles(mData@data$fc[mSelected, ], na.rm = T, probs = dropout$quantile) %>%
                                     rep(N_examples) %>%
                                     `*`(sample(c(1, dropout$fc_severity), N_genes * N_examples, T, c(1 - dropout$frequency, dropout$frequency))),
                                   ncol = N_examples)
        mPV[mSelected, ] <- matrix(rnorm(N_genes * N_examples, pv$mean, pv$sd) +
                                     sample(c(0, dropout$p_severity), N_genes * N_examples, T, c(1 - dropout$frequency, dropout$frequency)),
                                   N_genes, N_examples) %>% pmax(0) %>% pmin(1)
        mMean[mSelected, ] <- matrix(rnorm(N_genes * N_examples, mGenes.dist$mean[1:N_genes], mGenes.dist$sd[1:N_genes]), N_genes, N_examples)
        mZ[mSelected, ] <- (mFC[mSelected, ] - mData@gene.meta$dist.Mean[mSelected]) / mData@gene.meta$dist.SD[mSelected]
        mPVZ[mSelected, ] <- mZ[mSelected, ] %>% {
          tmp <- mPV[mSelected, ]
          tmp[tmp < 1e-20] <- 1e-20
          abs(.) * -log(tmp, 100)
        }
        
        mMeta <- data.table(ee.ID = mData@experiment.meta[, max(ee.ID) + 1] %>% seq(., length.out = N_examples),
                            ee.qScore = rnorm(N_examples, geeq$mean, geeq$sd) %>% pmax(0) %>% pmin(1),
                            ee.sScore = 1,
                            rsc.ID = experimentNames,
                            ee.Troubled = F,
                            ee.Public = T,
                            ee.Name = paste0('DOPE', 1:N_examples),
                            ee.Source = 'DOPE',
                            ee.Scale = 'LOG2',
                            ee.NumSamples = 1, # TODO These aren't used in search/enrich so we don't need them (yet)
                            ee.TagLongUri = NA,
                            ad.Name = 'DOPE',
                            ad.Company = 'DOPE',
                            ad.Sequencing = T,
                            sf.Subset = F,
                            sf.Cat = NA,
                            sf.CatLongUri = NA,
                            sf.ValLongUri = NA,
                            mConditions[sample(1:N_conditions, N_examples, T)],
                            n.DE = sum(mPV < 0.05, na.rm = T),
                            mean.fc = mean(mMean, na.rm = T))
        
        message('Precomputing cache...')
        DATA.HOLDER$dope <<- mData
        DATA.HOLDER$dope@experiment.meta <<- mMeta
        CACHE.BACKGROUND$dope <<- precomputeTags('dope') %>% reorderTags
        TAGS$dope <<- getTags('dope')
        
        mData@experiment.meta <- rbind(mData@experiment.meta, mMeta)
        mData@data$fc <- cbind(mData@data$fc, mFC)
        mData@data$meanval <- cbind(mData@data$meanval, mMean)
        mData@data$adj.pv <- cbind(mData@data$adj.pv, mPV)
        mData@data$zscore <- cbind(mData@data$zscore, mZ)
        mData@data$pvz <- cbind(mData@data$pvz, mPVZ)
        
        DATA.HOLDER$dope <<- mData
        rm(mData, mMeta, mFC, mMean, mPV, mZ, mPVZ)
        CACHE.BACKGROUND$dope <<- rbind(CACHE.BACKGROUND[[taxa]], CACHE.BACKGROUND$dope) %>% reorderTags
        TAGS$dope <<- rbind(TAGS[[taxa]], TAGS$dope)
        mName <- 'dope'
      }
      
      message('Searching...')
      search_N <- seq(min(N_genes, min_search), N_genes, by = gene_perms)
      searched <- lapply(search_N, function(i) { search(mGenes$entrez.ID[sample(1:N_genes, i)], mName, verbose = F) })
      message('Enriching...')
      enriched <- lapply(searched, function(i) {
        if(is.null(i) || nrow(i) == 0)
          data.table()
        else
          enrich(i, mName, verbose = F)
      })
      
      list(searched = searched,
           enriched = enriched,
           N_examples = N_examples,
           N_searched = search_N,
           N_genes = N_genes)
    })
  }),
  genes = mGenes,
  conditions = mConditions,
  fc = fc,
  pv = pv,
  dropout = dropout,
  geeq = geeq,
  N_conditions = N_conditions)
}

evalDope <- function(dopeResults, mTruthy = NULL) {
  if(is.null(mTruthy)) {
    cond <- dopeResults$conditions
    
    mBaseOptions <- parseListEntry(as.character(cond$cf.BaseLongUri))
    mValOptions <- parseListEntry(as.character(cond$cf.ValLongUri))
    
    mBaseOptions <- sapply(mBaseOptions, function(x) {
      mapped <- ONTOLOGIES.DEFS[Node_Long == x, data.table::first(Definition)] %>% as.character
      if(length(mapped) == 0)
        x
      else
        mapped
    }) %>% paste0(collapse = '; ')
    mValOptions <- sapply(mValOptions, function(x) {
      mapped <- ONTOLOGIES.DEFS[Node_Long == x, data.table::first(Definition)] %>% as.character
      if(length(mapped) == 0)
        x
      else
        mapped
    }) %>% paste0(collapse = '; ')
  } else {
    cond <- mTruthy
    mBaseOptions <- mTruthy[, unique(cf.BaseLongUri)]
    mValOptions <- mTruthy[, unique(cf.ValLongUri)]
  }
  
  lapply(dopeResults$data, function(examples) {
    lapply(examples, function(genes) {
      lapply(1:length(genes$enriched), function(N_genes) {
        if(is.null(genes$enriched[[N_genes]]) || nrow(genes$enriched[[N_genes]]) == 0)
          NULL
        else {
          tmp <- copy(genes$enriched[[N_genes]])[, I := .I]
          tmp2 <- copy(genes$searched[[N_genes]])[, I := .I]
          
          matched <- tmp[cf.Cat %in% cond$cf.Cat &
                            as.character(cf.BaseLongUri) %in% mBaseOptions &
                            as.character(cf.ValLongUri) %in% mValOptions]
          matched2 <- tmp2[grepl('DOPE', rn, fixed = T)]
          
          search_indices <- matched2[, I]
          search_indices_frac <- search_indices / nrow(tmp2)
          enrich_indices <- matched[, I]
          enrich_indices_frac <- enrich_indices / nrow(tmp)
          pvals <- matched[, pv.fisher]
          
          if(length(search_indices) > 1) {
            search_indices <- list(list(search_indices))
            search_indices_frac <- list(list(search_indices_frac))
          }
          
          if(length(enrich_indices) > 1) {
            enrich_indices <- list(list(enrich_indices))
            enrich_indices_frac <- list(list(enrich_indices_frac))
            pvals <- list(list(pvals))
          }

          data.table(meanfc = dopeResults$fc$mean,
                     meanpv = dopeResults$pv$mean,
                     dropout = dopeResults$dropout$frequency,
                     geeq = dopeResults$geeq$mean,
                     regularizer = dopeResults$regularizer,
                     n_examples_spiked = genes$N_examples, # Number of examples spiked in
                     n_searched_genes = genes$N_searched[N_genes], # Number of genes searched in this quest
                     n_total_genes = genes$N_genes,
                     n_condition_found = matched$A,
                     n_other_found = matched$B,
                     n_condition_corpus = matched$C,
                     n_other_corpus = matched$D,
                     search_indices = search_indices,
                     search_indices_frac = search_indices_frac,
                     enrich_indices = enrich_indices,
                     enrich_indices_frac = enrich_indices_frac,
                     pval = pvals)
        }
      }) %>% rbindlist
    }) %>% rbindlist
  }) %>% rbindlist
}

options(app.distance_cutoff = 2)
options(app.mfx = T, app.geeq = F, app.meanval = F, app.min.tags = 1.5)
updateOptions()

geneFreqs <- DATA.HOLDER$mouse@gene.meta[, list(list(entrez.ID)), n.DE] %>% setorder(n.DE)
condFreqs <- DATA.HOLDER$mouse@experiment.meta[, .N, .(cf.Cat, cf.CatLongUri, cf.Val, cf.ValLongUri, cf.Baseline, cf.BaseLongUri)] %>% setorder(N)

# ----
# 10 baseline conditions
mCond <- condFreqs[N == 10][2][, !'N'] %>%
  .[, .(cf.Cat = as.character(cf.Cat), cf.CatLongUri = as.character(cf.CatLongUri), cf.Val, cf.Baseline, cf.ValLongUri = as.character(cf.ValLongUri), cf.BaseLongUri = as.character(cf.BaseLongUri))]
mTruthy <- CACHE.BACKGROUND$mouse[rsc.ID ==
                                    CACHE.BACKGROUND$mouse[cf.Cat == mCond$cf.Cat &
                                                             cf.ValLongUri == mCond$cf.Val &
                                                             cf.BaseLongUri == mCond$cf.Baseline,
                                                           data.table::first(as.character(rsc.ID))]] %>%
  .[, .(cf.Cat = as.character(cf.Cat),
        cf.BaseLongUri = as.character(cf.BaseLongUri),
        cf.ValLongUri = as.character(cf.ValLongUri))]

#lapply(intersect(c(1, 10, 100, 1000), geneFreqs$n.DE), function(i) {
#  g <- unlist(geneFreqs[n.DE == i, V1])
#  lapply(c(0, 1, 3, 5, 10, 100), function(n) {
#    dope('mouse', geneIDs = sample(g, min(length(g), 20)), examples = n,
#         fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
#         dropout = list(p_severity = 0, fc_severity = NA, quantile = 0.95, frequency = 0.4), seed = 206134, mConditions = mCond) %>%
#      evalDope(mTruthy) %>% cbind(n.DE = i)
#  }) %>% rbindlist
#}) %>% rbindlist -> doped # TMP

lapply(intersect(1:100, geneFreqs$n.DE), function(i) {
  g <- unlist(geneFreqs[n.DE == i, V1])
  lapply(0:10, function(n) {
    dope('mouse', geneIDs = sample(g, min(length(g), 20)), examples = n,
         fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
         dropout = list(p_severity = 0, fc_severity = NA, quantile = 0.95, frequency = 0.4), seed = 206134, mConditions = mCond) %>%
      evalDope(mTruthy) %>% cbind(n.DE = i)
  }) %>% rbindlist
}) %>% rbindlist -> doped

saveRDS(doped, '/space/scratch/jsicherman/Thesis Work/data/bootstrap_10_0.4.rds')

# ----
lapply(intersect(1:100, geneFreqs$n.DE), function(i) {
  g <- unlist(geneFreqs[n.DE == i, V1])
  lapply(0:10, function(n) {
    dope('mouse', geneIDs = sample(g, min(length(g), 40)), examples = n,
         fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
         dropout = list(p_severity = 0, fc_severity = NA, quantile = 0.95, frequency = 0.4), seed = 206134, mConditions = mCond) %>%
      evalDope(mTruthy) %>% cbind(n.DE = i)
  }) %>% rbindlist
}) %>% rbindlist -> doped

saveRDS(doped, '/space/scratch/jsicherman/Thesis Work/data/bootstrap_10_0.4_gene40.rds')

# ----
# 21 baseline conditions
mCond <- condFreqs[N == 21][1][, !'N'] %>%
  .[, .(cf.Cat = as.character(cf.Cat), cf.CatLongUri = as.character(cf.CatLongUri), cf.Val, cf.Baseline, cf.ValLongUri = as.character(cf.ValLongUri), cf.BaseLongUri = as.character(cf.BaseLongUri))]
mTruthy <- CACHE.BACKGROUND$mouse[rsc.ID ==
                                    CACHE.BACKGROUND$mouse[cf.Cat == mCond$cf.Cat &
                                                             cf.ValLongUri == mCond$cf.Val &
                                                             cf.BaseLongUri == mCond$cf.Baseline,
                                                           data.table::first(as.character(rsc.ID))]] %>%
  .[, .(cf.Cat = as.character(cf.Cat),
        cf.BaseLongUri = as.character(cf.BaseLongUri),
        cf.ValLongUri = as.character(cf.ValLongUri))]

lapply(intersect(1:100, geneFreqs$n.DE), function(i) {
  g <- unlist(geneFreqs[n.DE == i, V1])
  lapply(0:10, function(n) {
    dope('mouse', geneIDs = sample(g, min(length(g), 20)), examples = n,
         fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
         dropout = list(p_severity = 0, fc_severity = NA, quantile = 0.95, frequency = 0.4), seed = 206134, mConditions = mCond) %>%
      evalDope(mTruthy) %>% cbind(n.DE = i)
  }) %>% rbindlist
}) %>% rbindlist -> doped

saveRDS(doped, '/space/scratch/jsicherman/Thesis Work/data/bootstrap_21_0.4.rds')

mCond <- condFreqs[N == 15][1][, !'N'] %>%
  .[, .(cf.Cat = as.character(cf.Cat), cf.CatLongUri = as.character(cf.CatLongUri), cf.Val, cf.Baseline, cf.ValLongUri = as.character(cf.ValLongUri), cf.BaseLongUri = as.character(cf.BaseLongUri))]
mTruthy <- CACHE.BACKGROUND$mouse[rsc.ID ==
                                    CACHE.BACKGROUND$mouse[cf.Cat == mCond$cf.Cat &
                                                             cf.ValLongUri == mCond$cf.Val &
                                                             cf.BaseLongUri == mCond$cf.Baseline,
                                                           data.table::first(as.character(rsc.ID))]] %>%
  .[, .(cf.Cat = as.character(cf.Cat),
        cf.BaseLongUri = as.character(cf.BaseLongUri),
        cf.ValLongUri = as.character(cf.ValLongUri))]

lapply(intersect(1:100, geneFreqs$n.DE), function(i) {
  g <- unlist(geneFreqs[n.DE == i, V1])
  lapply(0:10, function(n) {
    dope('mouse', geneIDs = sample(g, min(length(g), 20)), examples = n,
         fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
         dropout = list(p_severity = 0, fc_severity = NA, quantile = 0.95, frequency = 0.4), seed = 206134, mConditions = mCond) %>%
      evalDope(mTruthy) %>% cbind(n.DE = i)
  }) %>% rbindlist
}) %>% rbindlist -> doped

saveRDS(doped, '/space/scratch/jsicherman/Thesis Work/data/bootstrap_15_0.4.rds')

lapply(intersect(1:100, geneFreqs$n.DE), function(i) {
  g <- unlist(geneFreqs[n.DE == i, V1])
  lapply(0:10, function(n) {
    dope('mouse', geneIDs = sample(g, min(length(g), 10)), examples = n,
         fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
         dropout = list(p_severity = 0, fc_severity = NA, quantile = 0.95, frequency = 0.4), seed = 206134, mConditions = mCond) %>%
      evalDope(mTruthy) %>% cbind(n.DE = i)
  }) %>% rbindlist
}) %>% rbindlist -> doped

saveRDS(doped, '/space/scratch/jsicherman/Thesis Work/data/bootstrap_15_0.4_gene10.rds')

lapply(intersect(1:100, geneFreqs$n.DE), function(i) {
  g <- unlist(geneFreqs[n.DE == i, V1])
  lapply(0:10, function(n) {
    dope('mouse', geneIDs = sample(g, min(length(g), 100)), examples = n,
         fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
         dropout = list(p_severity = 0, fc_severity = NA, quantile = 0.95, frequency = 0.4), seed = 206134, mConditions = mCond) %>%
      evalDope(mTruthy) %>% cbind(n.DE = i)
  }) %>% rbindlist
}) %>% rbindlist -> doped

saveRDS(doped, '/space/scratch/jsicherman/Thesis Work/data/bootstrap_15_0.4_gene100.rds')

# ----
# 1 baseline condition
mCond <- condFreqs[N == 1][7263][, !'N'] %>%
  .[, .(cf.Cat = as.character(cf.Cat), cf.CatLongUri = as.character(cf.CatLongUri), cf.Val, cf.Baseline, cf.ValLongUri = as.character(cf.ValLongUri), cf.BaseLongUri = as.character(cf.BaseLongUri))]
mTruthy <- CACHE.BACKGROUND$mouse[rsc.ID ==
                                    CACHE.BACKGROUND$mouse[cf.Cat == mCond$cf.Cat &
                                                             cf.ValLongUri == mCond$cf.Val &
                                                             cf.BaseLongUri == mCond$cf.Baseline,
                                                           data.table::first(as.character(rsc.ID))]] %>%
  .[, .(cf.Cat = as.character(cf.Cat),
        cf.BaseLongUri = as.character(cf.BaseLongUri),
        cf.ValLongUri = as.character(cf.ValLongUri))]

lapply(intersect(1:100, geneFreqs$n.DE), function(i) {
  g <- unlist(geneFreqs[n.DE == i, V1])
  lapply(0:10, function(n) {
    dope('mouse', geneIDs = sample(g, min(length(g), 20)), examples = n,
         fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
         dropout = list(severity = 0, frequency = 0), seed = 206134, mConditions = mCond) %>%
      evalDope(mTruthy) %>% cbind(n.DE = i)
  }) %>% rbindlist
}) %>% rbindlist -> doped

saveRDS(doped, '/space/scratch/jsicherman/Thesis Work/data/bootstrap_1.rds')

plotData <- function(N, doped = NULL, suffix = N, n = 2, xlim = NULL, ylim = NULL, col = enrich_indices_frac,
                      title = 'Sensitivity Analysis', subtitle = waiver()) {
  if(is.null(doped))
    doped <- readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/bootstrap_', suffix, '.rds'))
  
  first <- function(x, n) {
    head(unlist(x), n)
  }
  
  tmp <- substitute(col)
  
  doped[, .(frac = first(eval(tmp), n)),
                 .(n.DE, n_examples_spiked)] %>%
    .[is.na(frac), frac := if_else(length(grep('frac', as.character(tmp))) > 0, 1, max(frac, na.rm = T))] %>%
    ggplot(aes(n_examples_spiked, frac,
               group = n_examples_spiked)) +
    geom_jitter() +
    theme(legend.justification = c('right', 'top'),
          legend.position = c(0.99, 0.99),
          legend.box.just = 'right',
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.margin = margin(1, 1, 1, 1)) +
    theme_pavlab() +
    scale_y_pavlab(name = as.character(tmp),
                   as_percentages = length(grep('frac', as.character(tmp))) > 0) +
    scale_x_pavlab(name = '# of Examples Spiked In',
                   breaks = doped$n_examples_spiked %>% unique, limits = xlim) +
    #scale_color_fermenter(name = 'Gene DE Count', palette = 'Dark2') +
    #scale_color_pavlab(type = 'brewer', name = 'Gene DE Count', palette = 'Dark2') +
    ggtitle(title, subtitle)
}

doped <- readRDS('/space/scratch/jsicherman/Thesis Work/data/bootstrap_10_0.4.rds')
plotData(doped = doped, col = search_indices_frac)

doped <- readRDS('/space/scratch/jsicherman/Thesis Work/data/bootstrap_15_0.4.rds')
plotData(doped = doped, col = search_indices_frac)

doped <- readRDS('/space/scratch/jsicherman/Thesis Work/data/bootstrap/10_0.01_0.9_20.rds')
plotData(doped = doped, col = enrich_indices_frac)
