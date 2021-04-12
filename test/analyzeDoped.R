cMap <- DATA.HOLDER$human@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
  unique %>% .[, cf.ID := 1:nrow(.)] %>%
  merge(DATA.HOLDER$human@experiment.meta[, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
        by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'))

cMap <- CACHE.BACKGROUND$human[, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
  merge(cMap[, .(rsc.ID, cf.ID)], by = 'rsc.ID') %>%
  .[, !'rsc.ID']

rm(DATA.HOLDER, CACHE.BACKGROUND)

tmp <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/bootstrap_scores.rds')

dope_scores <- lapply(1:length(tmp), function(i) {
  data.table(group = i,
             insert_sd = tmp[[i]]$mean_sd,
             insert_mean_index = tmp[[i]]$mean_index,
             insert_groups = tmp[[i]]$groups[tmp[[i]]$associations],
             search_indices = tmp[[i]]$indices,
             search_scores = tmp[[i]]$scores)
}) %>% rbindlist(fill = T)

dope_scores %>% na.omit %>%
  ggplot(aes(insert_mean_index, search_indices)) +
  geom_point(aes(color = factor(insert_mean_index))) +
  stat_smooth(method = 'lm', se = F) +
  theme_classic(base_size = 20) + theme(legend.position = 'none') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Mean insertion index', y = 'Index of occurrence')

library(parallel)
options(mc.cores = 3)
dope_enrich <- mclapply(1:length(tmp), function(i) {
  print(i)
  if(!is.null(tmp[[i]]$enrich) && nrow(tmp[[i]]$enrich) > 0)
    enriched <- tmp[[i]]$enrich[, f := I/max(I)] %>%
      merge(cMap, by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri')) %>%
      .[cf.ID %in% unique(tmp[[i]]$groups), .(stat, distance, I, f, cf.ID)]
  else
    enriched <- data.table(stat = NA, distance = NA, I = NA, f = NA,
                           cf.ID = tmp[[i]]$groups[tmp[[i]]$associations])
  
  data.table(group = i,
             insert_sd = tmp[[i]]$mean_sd,
             n_exp = tmp[[i]]$n_exp,
             n_groups = length(tmp[[i]]$groups),
             n_genes = length(tmp[[i]]$genes),
             insert_mean_index = tmp[[i]]$mean_index,
             enriched)
}) %>% rbindlist(fill = T)
rm(tmp)

dope_enrich[distance == 0] %>%
  .[, c('ms', 'mi', 'mf') := list(mean(stat, na.rm = T),
                                  mean(I, na.rm = T),
                                  mean(f, na.rm = T)), .(group, contrast)] %>%
  .[, .(group, contrast, insert_mean_index,
        n_exp = paste0(n_exp, ' Experiments'),
        n_groups,
        n_genes = paste0(n_genes, ' Genes'), ms = log10(ms), mi, mf)] -> dope_enrich_plottable

library(scattermore)
dope_enrich %>%
  ggplot(aes(insert_mean_index, log10(1 + stat))) +
  geom_scattermore(color = '#002145') +
  stat_smooth(method = 'lm', color = 'red', formula = y ~ poly(x, 2)) +
  theme_classic(base_size = 20) +
  scale_x_continuous(expand = expansion(add = 0.01),
                     breaks = c(0, 0.5, 1),
                     labels = c('Start', 'Middle', 'End'),
                     name = 'Spike in position') +
  scale_y_continuous(expand = c(0, 0),
                     name = expression(log[10](1+stat)),
                     limits = c(-2, 2.5))

dope_enrich %>% copy %>%
  .[, x := factor(plyr::round_any(100 * (1-insert_mean_index), 10))] %>%
  .[, x_f := factor(plyr::round_any(n_exp, 10))] %>%
  .[x_f %in% c(0, 10, 40, 50)] %>%
  #.[, x := paste0(x, '^{"th"}')] %>%
  #.[, x := factor(x, levels = paste0(seq(0, 100, by = 5), '^{"th"}'), ordered = T)] %>%
  ggplot(aes(x, f, fill = x_f)) +
  geom_boxplot(outlier.alpha = 0.3, lwd = 0.4) +
  theme_classic(base_size = 20) +
  theme(plot.margin = unit(c(12, 33, 0, 0), 'points'),
        legend.position = c(1, 1), legend.justification = c(1, 1),
        legend.background = element_rect('white', 'black')) +
  scale_x_discrete(breaks = c(0, 50, 100), labels = c('Unrelated', 'Somewhat related', 'Highly related'), expand = c(0, 0)) +# paste0(seq(0, 100, by = 50), '^{"th"}'), labels = rlang::parse_exprs) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c('First', 'Mid', 'Last'), expand = expansion(add = c(0.005, 0))) +
  scale_fill_brewer(name = '# of Spike-ins', palette = 'Dark2', labels = c('[1, 5]', '(5, 15)', '[35, 45]', '(45, 50]')) + # '[15, 25]', '(25, 35)',
  labs(x = 'True relatedness', y = 'Calculated rank')

dope_enrich[distance == 0] %>%
  ggplot(aes(f, log10(stat))) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  theme_classic(base_size = 20) +
  #ylim(0, 30) +
  labs(x = 'Mean occurrence of contrast', y = 'Mean test statistic') +
  scale_x_continuous(labels = scales::percent)

dope_enrich[distance == 0] %>%
  ggplot(aes(n_genes, log10(stat))) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  theme_classic(base_size = 20) +
  labs(x = 'Number of genes searched', y = 'Mean test statistic')
