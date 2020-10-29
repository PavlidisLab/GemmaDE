library(dplyr)
library(data.table)
library(fitdistrplus)
library(pscl)
library(countreg)

lapply(c('human', 'mouse', 'rat'), function(taxa) {
  taxa <- 'human'
  BOOTSTRAP <- readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/bootstrap_', taxa, '.rds'))
  
  lapply(BOOTSTRAP, function(x) data.table(rsc.ID = rownames(x$search), score = x$search$score)) %>%
    rbindlist %>% .[!is.na(score)] -> scores
  
  # Fit experiment scores to a normal distribution
  exp.stats <- scores[, .(mean = mean(score), median = Rfast::med(score), stdev = Rfast::Var(score, T)), rsc.ID]
  saveRDS(exp.stats, paste0('/space/scratch/jsicherman/Thesis Work/data/stats/experiments_', taxa, '.rds'))
  
  lapply(BOOTSTRAP, function(x) melt(x$search[, !c('score', 'direction')], value.name = 'score', variable.name = 'gene.Name')) %>%
    rbindlist %>% .[!is.na(score)] -> scores
  
  # Fit gene scores to an exponential distribution (although may be better as a zinbinom)
  gene.stats <- scores[, .(rate = tryCatch(fitdist(score, 'exp')$estimate, error = function(e) { integer(0) })), gene.Name]
  saveRDS(gene.stats, paste0('/space/scratch/jsicherman/Thesis Work/data/stats/genes_', taxa, '.rds'))
  
  lapply(BOOTSTRAP, function(x) x$enrich[, .(cf.BaseLongUri, cf.ValLongUri, pv.fisher)]) %>%
    rbindlist %>% .[!is.na(pv.fisher)] -> BOOTSTRAP
  
  fitTerms <- function(data) {
    obs <- round(1e5 * (1 - data$pv.fisher))
    data.table(mu = mean(obs),
               theta = tryCatch(zeroinfl(obs ~ 1 | 1, dist = 'negbin')$theta, error = function(e) { integer(0) }),
               pi = sum(obs == 0) / length(obs))
  }
  
  # Fit transformed term scores to a zinbinom distribution. Test statistic is 1e5 * (1 - pv.fisher)
  term.stats <- BOOTSTRAP[, fitTerms(.SD), .(cf.BaseLongUri, cf.ValLongUri)]
  saveRDS(term.stats, paste0('/space/scratch/jsicherman/Thesis Work/data/stats/terms_', taxa, '.rds'))
  rm(exp.stats, gene.stats, term.stats, BOOTSTRAP)
})
