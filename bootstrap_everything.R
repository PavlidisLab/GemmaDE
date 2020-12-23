library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(htmlwidgets)
library(sparkline)
library(DT)
library(heatmaply)
library(shinyHeatmaply)
library(shinypanels)
library(RColorBrewer)
library(stringr)

library(jsonlite)
library(async)

library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)

library(gemmaAPI)
library(ermineR)
library(mygene)

if(!exists('DATA.HOLDER'))
  source('dependencies.R')

DATA.HOLDER$human <- NULL
DATA.HOLDER$rat <- NULL
CACHE.BACKGROUND$human <- NULL
CACHE.BACKGROUND$rat <- NULL

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

lapply(c(25, 20, 31, 10, 41, 5, 10, 50, 80), function(condN) {
  while(TRUE) {
    mCond <- condFreqs[N == condN] %>% .[sample(1:nrow(.), 1)] %>% .[, !'N'] %>%
      .[, .(cf.Cat = as.character(cf.Cat),
            cf.CatLongUri = as.character(cf.CatLongUri),
            cf.Val, cf.Baseline,
            cf.ValLongUri = as.character(cf.ValLongUri), cf.BaseLongUri = as.character(cf.BaseLongUri))]
    
    mID <- CACHE.BACKGROUND$mouse[cf.Cat == mCond$cf.Cat &
                                    cf.ValLongUri == mCond$cf.Val &
                                    cf.BaseLongUri == mCond$cf.Baseline,
                                  data.table::first(as.character(rsc.ID))]
    if(length(mID) != 0) {
      mTruthy <- CACHE.BACKGROUND$mouse[rsc.ID == mID] %>%
        .[, .(cf.Cat = as.character(cf.Cat),
              cf.BaseLongUri = as.character(cf.BaseLongUri),
              cf.ValLongUri = as.character(cf.ValLongUri))]
      break
    }
  }
  
  lapply(c(20, 10, 5, 2, 50, 100), function(geneN) {
    lapply(c(0.9, 0.95, 0.99), function(fc) {
      lapply(c(0.001, 0.01, 0.05), function(pv) {
        lapply(c(0.2, 0.3, 0.4, 0.5, 0.75, 0.8, 0.9), function(freq) {
          lapply(intersect(c(1, 2, 5, 30, 50, 100, 300, 1000), geneFreqs$n.DE), function(i) {
            g <- unlist(geneFreqs[n.DE == i, V1])
            lapply(0:10, function(exN) {
              dope('mouse', geneIDs = sample(g, min(length(g), geneN)), examples = exN,
                   fc = list(mean = 1, sd = 0), pv = list(mean = pv, sd = 0),
                   dropout = list(p_severity = 0, fc_severity = NA, quantile = fc, frequency = freq), seed = 206134, mConditions = mCond) %>%
                evalDope(mTruthy) %>% cbind(n.DE = i)
            }) %>% rbindlist
          }) %>% rbindlist -> doped
          
          saveRDS(doped, paste0('/space/scratch/jsicherman/Thesis Work/data/bootstrap/condN_',
                                condN, '_pv_', pv, '_freq_', freq, '_quantile_', fc, '_gene_', geneN, '.rds'))
          rm(doped)
          gc()
          NULL
        })
      })
    })
  })
})
