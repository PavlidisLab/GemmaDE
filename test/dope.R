source('test/setupDope.R')

mCond <- condFreqs[N == 56 & distance <= getOption('app.distance_cutoff')] %>% .[sample(1:nrow(.), 1)] %>% .[, !'N']

mRSC <- CACHE.BACKGROUND$mouse[mCond$cf.Cat == cf.Cat &
                                 mCond$cf.ValLongUri == cf.ValLongUri &
                                 mCond$cf.BaseLongUri == cf.BaseLongUri &
                                 mCond$distance == distance, as.character(data.table::first(rsc.ID))]

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

lapply(intersect(c(1, 5, 10, 100, 1000), geneFreqs$n.DE), function(i) {
  g <- unlist(geneFreqs[n.DE == i, V1])
  lapply(0:5, function(n) {
    gid <- sample(g, min(length(g), 20))
    lapply(206134:206140, function(seed) {
      lapply(seq(0.4, 0.9, 0.1), function(freq) {
        dope('mouse', geneIDs = gid, examples = n,
             fc = list(mean = 2, sd = 0), pv = list(mean = 0.01, sd = 0),
             dropout = list(p_severity = NA, fc_severity = NA, quantile = 0.95, frequency = freq), seed = seed, mConditions = mCond) %>%
          evalDope(mTruthy) %>% cbind(n.DE = i) %>% cbind(freq = freq)
      }) %>% rbindlist
    }) %>% rbindlist
  }) %>% rbindlist
}) %>% rbindlist -> doped

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

plotData2 <- function(doped = NULL, n = Inf, cond = NA,
                      xlim = NULL, ylim = NULL, col = enrich_indices_frac,
                      title = 'Sensitivity Analysis', subtitle = waiver()) {
  first <- function(x, n) {
    if(length(grep('search', as.character(tmp))) > 0)
      head(unlist(x), n)
    else
      x
  }
  
  if(is.null(cond))
    cond <- doped$cond[, c(cf.Baseline, cf.Val)]
  doped <- doped$doped
  
  if(!is.na(cond))
    doped <- doped[baseline == cond[1] & contrast == cond[2]]
  
  tmp <- substitute(col)
  
  doped[n.DE %in% c(100, 300, 1000), .(frac = first(eval(tmp), n)), .(n_examples_spiked, n.DE)] %>%
    .[is.na(frac), frac := if_else(length(grep('frac', as.character(tmp))) > 0, 1, max(frac, na.rm = T))] %>%
    ggplot(aes(n_examples_spiked, frac, group = n_examples_spiked)) +
    geom_boxplot() +
    # stat_summary(fun = median, geom = 'line') +
    theme(legend.justification = c('right', 'top'),
          legend.position = c(0.99, 0.99),
          legend.box.just = 'right',
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.margin = margin(1, 1, 1, 1)) +
    theme_pavlab(axis_x_angle = 90) +
    scale_color_discrete() +
    scale_y_pavlab(name = as.character(tmp),
                   as_percentages = length(grep('frac', as.character(tmp))) > 0) +
    scale_x_pavlab(name = '# of Examples Spiked In',
                   breaks = doped$n_examples_spiked %>% unique, limits = xlim) +
    facet_wrap(~n.DE, scales = 'free_y', ncol = 1) +
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

doped <- readRDS('/space/scratch/jsicherman/Thesis Work/data/bootstrap/condN_25_pv_0.001_freq_0.2_quantile_0.9_gene_20.rds')
plotData(doped = doped, col = enrich_indices)

doped <- readRDS('/space/scratch/jsicherman/Thesis Work/data/bootstrap/condN_25_pv_0.001_freq_0.5_quantile_0.95_gene_20.rds')
plotData(doped = doped, col = enrich_indices_frac, n = 1)

doped <- readRDS('/space/scratch/jsicherman/Thesis Work/data/bootstrap_everything_all/condN_25_pv_0.001_freq_0.01_quantile_0.9_gene_20.rds')
plotData2(doped, cond = NULL, col = enrich_indices_frac)

doped <- readRDS('/space/scratch/jsicherman/Thesis Work/data/bootstrap_everything_all/condN_25_pv_0.05_freq_0.5_quantile_0.9_gene_20.rds')

