source('test/setupDope.R')

lapply(c(25, 20, 31, 10, 41, 5, 10, 50, 80), function(condN) {
  while(TRUE) {
    mCond <- condFreqs[N == condN & distance <= getOption('app.distance_cutoff')] %>% .[sample(1:nrow(.), 1)] %>% .[, !'N']
    
    mRSC <- CACHE.BACKGROUND$mouse[mCond$cf.Cat == cf.Cat &
                                     mCond$cf.ValLongUri == cf.ValLongUri &
                                     mCond$cf.BaseLongUri == cf.BaseLongUri &
                                     mCond$distance == distance, as.character(data.table::first(rsc.ID))]
    
    if(length(mRSC) != 0) {
      mCond <- mCond %>% setnames(c('cf.BaseLongUri', 'cf.ValLongUri'), c('cf.Baseline', 'cf.Val')) %>%
        cbind(cf.BaseLongUri = data.table::first(ONTOLOGIES.DEFS[Definition == mCond$cf.Baseline, Node_Long]) %>% as.character) %>%
        cbind(cf.ValLongUri = data.table::first(ONTOLOGIES.DEFS[Definition == mCond$cf.Val, Node_Long]) %>% as.character) %>%
        cbind(cf.CatLongUri = data.table::first(DATA.HOLDER$mouse@experiment.meta[rsc.ID == mRSC, cf.CatLongUri]) %>% as.character) %>%
        .[, .(cf.Cat = as.character(cf.Cat), cf.CatLongUri, cf.Val = as.character(cf.Val), cf.Baseline = as.character(cf.Baseline),
              cf.ValLongUri, cf.BaseLongUri)]
      
      mTruthy <- CACHE.BACKGROUND$mouse[rsc.ID == mRSC & distance <= getOption('app.distance_cutoff'),
                                        .(cf.Cat = as.character(cf.Cat),
                                          cf.BaseLongUri = as.character(cf.BaseLongUri),
                                          cf.ValLongUri = as.character(cf.ValLongUri))]
      
      break
    }
  }
  
  lapply(c(20, 10, 5, 2, 50, 100), function(geneN) {
    lapply(c(1e6, 100, 10, 1, 0.1, 0.01, -0.01, -0.1, -1, -10, -100, -1e6), function(regularizer) {
          lapply(intersect(c(1, 2, 5, 30, 50, 100, 300, 1000, 2000, 2500, 3214, 3695), geneFreqs$n.DE), function(i) {
            g <- unlist(geneFreqs[n.DE == i, V1])
            lapply(c(0:15, 20, 25, 30, 35, 40, 45, 50, 100), function(exN) {
              lapply(206134:206140, function(seed) {
                message(paste0('Examples: ', exN, ' n.DE: ', i, ' reg: ', regularizer, ' genes: ', geneN))
                dopeScores('mouse', geneIDs = sample(g, min(length(g), geneN)), examples = exN,
                     seed = seed, mConditions = mCond, regularizer = regularizer) %>%
                  evalDope(mTruthy) %>% cbind(n.DE = i)
              }) %>% rbindlist
            }) %>% rbindlist
          }) %>% rbindlist -> doped
          
          saveRDS(list(cond = mCond,
                       doped = doped),
                  paste0('/space/scratch/jsicherman/Thesis Work/data/bootstrap_scores_all/condN_',
                         condN, '_reg_', regularizer, '_gene_', geneN, '.rds'))
          rm(doped)
          gc()
          NULL
    })
  })
})
