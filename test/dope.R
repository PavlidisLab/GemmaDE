dope <- function(taxa,
                 genes,
                 examples,
                 N_conditions = 1,
                 fc = list(mean = 2, sd = 0.2),
                 pv = list(mean = 0.001, sd = 0.05),
                 dropout = list(severity = 0.2, frequency = 0.1),
                 geeq = list(mean = 1, sd = 0)) {
  # Dope the real data by adding `N_examples` new pseudo-experiments
  # that differentially express the same `N_genes` at the parameters:
  # logFC ~ N(fc$mean, fc$sd)
  # adj.pv ~ N(pv$mean, pv$sd)
  # GEEQ ~ N(geeq$mean, geeq$sd)
  # mean expression level ~ N(meanval$mean, meanval$sd) for the selected genes
  # 
  # And are annotated with the same `N_conditions` conditions
  # Then search/enrich the data for the selected genes
  
  mData <- DATA.HOLDER[[taxa]]
  mData@taxon <- 'dope'
  mData@go <- data.table()
  
  # Select genes that we'll work with
  nGenes <- nrow(mData@gene.meta)
  
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
  mConditions <- mData@experiment.meta[, .(cf.Cat, cf.CatLongUri, cf.Val, cf.ValLongUri, cf.Baseline, cf.BaseLongUri)] %>%
    unique %>% .[sample(1:.N, N_conditions)]
  
  lapply(examples, function(N_examples) {
    message(paste0('Simulating ', N_examples, ' examples...'))
    lapply(genes, function(N_genes) {
      message(paste0('Simulating ', N_genes, ' genes...'))
      
      mSelected <- nSelected[1:N_genes]
      
      mData <- DATA.HOLDER[[taxa]]
      mData@taxon <- 'dope'
      mData@go <- data.table()
      
      # Add `N_examples` rows to experiment.meta
      # Add `N_examples` columns to data$fc, data$meanval, data$adj.pv, data$zscore, data$pvz
      mFC <- mMean <- mZ <- mPVZ <- matrix(0, nGenes, N_examples, dimnames = list(NULL, paste0('RSCID.DOPE.', 1:N_examples)))
      mPV <- matrix(1, nGenes, N_examples, dimnames = list(NULL, paste0('RSCID.DOPE.', 1:N_examples)))
      mFC[mSelected, ] <- matrix(rnorm(N_genes * N_examples, fc$mean, fc$sd), N_genes, N_examples)
      mPV[mSelected, ] <- matrix(rnorm(N_genes * N_examples, pv$mean, pv$sd) +
                                   sample(c(0, dropout$severity), N_genes * N_examples, T, c(1 - dropout$frequency, dropout$frequency)),
                                 N_genes, N_examples) %>%
        pmax(0) %>% pmin(1)
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
                          rsc.ID = paste0('RSCID.DOPE.', 1:N_examples),
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
                          mConditions,
                          n.DE = sum(mPV < 0.05, na.rm = T),
                          mean.fc = mean(mMean, na.rm = T))
      
      message('Precomputing cache...')
      DATA.HOLDER$dope <<- mData
      DATA.HOLDER$dope@experiment.meta <<- mMeta
      CACHE.BACKGROUND$dope <<- precomputeTags('dope') %>% reorderTags
      
      mData@experiment.meta <- rbind(mData@experiment.meta, mMeta)
      mData@data$fc <- cbind(mData@data$fc, mFC)
      mData@data$meanval <- cbind(mData@data$meanval, mMean)
      mData@data$adj.pv <- cbind(mData@data$adj.pv, mPV)
      mData@data$zscore <- cbind(mData@data$zscore, mZ)
      mData@data$pvz <- cbind(mData@data$pvz, mPVZ)
      
      DATA.HOLDER$dope <<- mData
      rm(mData, mMeta, mFC, mMean, mPV, mZ, mPVZ)
      CACHE.BACKGROUND$dope <<- rbind(CACHE.BACKGROUND[[taxa]], CACHE.BACKGROUND$dope) %>% reorderTags
      
      message('Searching...')
      searched <- search(mGenes$entrez.ID, 'dope', verbose = F)
      message('Enriching...')
      enriched <- enrich(searched, 'dope', verbose = F)
      
      list(searched = searched,
           enriched = enriched,
           genes = mGenes,
           conditions = mConditions,
           N_examples = N_examples,
           N_genes = N_genes)
    })
  })
}

options(app.distance_cutoff = Inf)
updateOptions()
doped <- dope('mouse', seq(1, 20, by = 3), seq(1, 10, by = 2), 1)
