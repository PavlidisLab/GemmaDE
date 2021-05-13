source('/home/jsicherman/Thesis Work/requirements.R')

source('/home/jsicherman/Thesis Work/dependencies.R')

library(parallel)

contrasts <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/contrast_aff.rds')

DATA.HOLDER[c('human', 'mouse', 'rat')] <- NULL
NULLS.EXP[c('human', 'mouse', 'rat')] <- NULL

options(mc.cores = 15)
mclapply(1:5000, function(n) {
  system(paste0('echo ', n))
  lapply(1:10, function(i) {
    GOIs <- sample(200 * floor(n/200) + (1:200), i)
    contrasts[entrez.ID %in% GOIs] %>%
      .[, .(probability = sum(probability)), contrast] %>%
      merge(
        merge(DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
              CACHE.BACKGROUND$artificial[, .(rsc.ID, cf.BaseLongUri, cf.ValLongUri)],
              by = 'rsc.ID') %>% .[, !'rsc.ID'] %>% unique,
        by.x = 'contrast', by.y = 'cf.ID', all = T) %>%
      unique %>%
      setorder(-probability) %>%
      .[, rank := rank(-probability, ties.method = 'min')] %>%
      .[, rank := as.integer(as.character(factor(rank, labels = 1:length(unique(rank)))))] -> tmp
    
    enrich(GOIs, getConfig(taxa = 'artificial'), keepopen = 1) %>%
      copy %>% .[, I := .I] %>%
      merge(tmp[, .(cf.BaseLongUri, cf.ValLongUri, rank, probability)],
            by = c('cf.BaseLongUri', 'cf.ValLongUri')) %>%
      .[, .(stat = max(stat), I = min(I),
            zscore = max(zscore),
            rank = min(rank),
            probability = max(probability)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
      setorder(I) %>% {
        data.table(genes = i,
                   cor = cor(.$stat, .$probability, method = 'spearman'),
                   n = .[rank <= 5, sum(I <= 50)])
      }
  }) %>% rbindlist
}) %>% rbindlist -> new.rankings2

rankings %>% copy %>%
  .[, n := pmin(n, 10)] %>%
  .[, meanline := mean(n)] %>%
  ggplot(aes(n)) +
  geom_histogram(fill = 'black') +
  geom_vline(aes(xintercept = meanline), color = 'red', linetype = 'dashed') +
  scale_x_continuous(name = 'Number of top 5 contrasts on first page',
                     n.breaks = 4, expand = c(0, 0)) +
  scale_y_continuous(name = 'Count', expand = c(0, 0)) +
  theme_classic(base_size = 20)

rankings %>%
  ggplot(aes(cor)) +
  geom_histogram(color = 'black', fill = 'white') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 20) +
  labs(x = 'Correlation between Pr(DE | contrast) and Enrichment Rank', y = 'Count')
