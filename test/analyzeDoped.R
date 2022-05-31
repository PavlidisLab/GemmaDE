tmp <- readRDS(paste(DATADIR, 'artificial/bootstrapped_scores.rds', sep='/'))

mSimpleCache <- CACHE.BACKGROUND$human[, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, ID = paste0(cf.Cat, cf.BaseLongUri, cf.ValLongUri))]
dope_scores <- lapply(1:length(tmp), function(i) {
  if(is.null(tmp[[i]])) return(NULL)
  
  mFetch <- tmp[[i]]$enrich %>% merge(tmp[[i]]$diff, by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all = T) %>%
    merge(
      mSimpleCache[ID %in% .[, paste0(cf.Cat, cf.BaseLongUri, cf.ValLongUri)], .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, N)],
      by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all.x = T) %>%
    .[is.na(N), N := 0]
  
  if(nrow(mFetch) == 0)
    mFetch <- data.table(N = NA, f = NA, distance = NA, score = NA, scoreD = NA, fD = NA)
  
  data.table(group = i,
             insert = tmp[[i]]$insert,
             n_exp = tmp[[i]]$n_exp,
             n_genes = length(tmp[[i]]$genes),
             n_baseline = mFetch$N,
             f = mFetch$f,
             N = mFetch$N,
             distance = mFetch$distance,
             score = mFetch$score,
             scoreD = mFetch$scoreDelta,
             fD = mFetch$fDelta)
}) %>% rbindlist(fill = T)
rm(tmp)
saveRDS(dope_scores, paste(DATADIR, 'artificial/bootstrapped_scores_processed.rds', sep='/'))
