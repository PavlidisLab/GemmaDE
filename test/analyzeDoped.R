tmp <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/bootstrapped_scores.rds')

mSimpleCache <- CACHE.BACKGROUND$human[, .N, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, ID = paste0(cf.Cat, cf.BaseLongUri, cf.ValLongUri))]
dope_scores <- lapply(1:length(tmp), function(i) {
  if(is.null(tmp[[i]])) return(NULL)
  
  mFetch <- tmp[[i]]$enrich %>% merge(tmp[[i]]$diff, by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all = T) %>%
    merge(
      mSimpleCache[ID %in% .[, paste0(cf.Cat, cf.BaseLongUri, cf.ValLongUri)], .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, N)],
      by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all.x = T) %>%
    .[is.na(N), N := 0]# %>%
    #.[is.na(fDelta), fDelta := 0] %>%
    #.[is.na(scoreDelta), scoreDelta := 0]
  
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
saveRDS(dope_scores, '/space/scratch/jsicherman/Thesis Work/data/artificial/bootstrapped_scores_processed.rds')

library(gridExtra)
dope_scores %>%
  .[distance == 0] %>%
  #.[, bin := as.factor(plyr::round_any(insert, 0.2))] %>%
  .[, exp_bin := factor(plyr::round_any(n_exp, 10), labels = c('1-5', mapply(paste, seq(6, 45, by = 10), seq(15, 50, by = 10), sep = '-'), '46-50'))] %>%
  ggplot(aes(insert, f, group = exp_bin, color = exp_bin)) +
  stat_smooth(method = 'loess', linetype = 0, alpha = 0.15) +
  stat_smooth(method = 'loess', se = F) + theme_bw(20) +
  scale_x_continuous(expand = expansion(add = 0.001), breaks = c(0, 0.5, 1),
                     labels = c('100th', '50th', '0th')) +
  scale_y_continuous(expand = c(0, 0), name = 'Contrast rank', breaks = c(0, 0.5),
                     labels = c('Best', 'Mid')) +
  scale_color_brewer(palette = 'Dark2', name = '# of Spike ins') -> panelA

dope_scores %>%
  .[distance <= 1] %>%
  .[, bin := as.factor(plyr::round_any(insert, 0.05))] %>%
  .[, exp_bin := factor(plyr::round_any(n_exp, 10), labels = c('1-5', mapply(paste, seq(6, 45, by = 10), seq(15, 50, by = 10), sep = '-'), '46-50'))] %>%
  ggplot(aes(bin, fD, color = exp_bin)) +
  stat_summary() +
  #scale_x_continuous(expand = expansion(add = 0.01), breaks = seq(0, 1, by = 0.25), labels = paste0(seq(100, 0, by = -25), 'th')) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_color_brewer(palette = 'Purples', name = '# of Spike ins') +
  scale_size_discrete(name = '# of Spike ins') +
  #scale_y_continuous(limits = c(0, 0.3)) +
  scale_x_discrete(breaks = seq(0, 1, by = 0.25), labels = paste0(seq(100, 0, by = -25), 'th')) +
  theme_bw(base_size = 20) +
  labs(x = 'Spike in specificity (percentile)', y = expression(paste(Delta, 'Contrast rank'))) -> panelB

dope_scores %>%
  .[distance <= 1] %>%
  .[, bin := as.factor(plyr::round_any(insert, 0.2))] %>%
  .[, exp_bin := factor(plyr::round_any(n_exp, 10), labels = c('1-5', mapply(paste, seq(6, 45, by = 10), seq(15, 50, by = 10), sep = '-'), '46-50'))] %>%
  .[, rowMax := ifelse(insert >= 0.5, 3, 12)] %>%
  ggplot(aes(f, color = exp_bin)) + geom_density(lwd = 1.1) +
  geom_blank(aes(y = rowMax)) +
  facet_wrap(~bin, scales = 'free_y', labeller = labeller(.default = function(x) {
    c('90th - 100th','70th - 90th', '50th - 70th', '30th - 50th', '10th - 30th', '0th - 10th')
  })) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c('Best', 'Mid', 'Worst'), expand = expansion(add = 0.005)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_brewer('# of Spike ins', palette = 'Purples') +
  scale_fill_brewer('# of Spike ins', palette = 'Purples') +
  theme_bw(20) + labs(x = 'Contrast Ranking', y = 'Density') +
  theme(strip.background = element_blank(), plot.margin = unit(c(5.5, 25, 5.5, 5.5), 'pt'), legend.position = 'none', strip.text = element_text(hjust = 0, face = 'bold')) -> panelC

cowplot::plot_grid(panelC,
                   cowplot::plot_grid(panelA + labs(x = element_blank()), panelB, align = 'v', labels = c('B', 'C'), ncol = 1),
                   labels = c('A', '', ''), ncol = 2)



panelA <- dope_scores %>%
  .[distance <= 1] %>%
  #.[is.na(stat), stat := 0] %>%
  .[, bin := as.factor(plyr::round_any(insert, 0.05))] %>%
  .[, exp_bin := factor(plyr::round_any(n_exp, 10), labels = c('1-5', mapply(paste, seq(6, 45, by = 10), seq(15, 50, by = 10), sep = '-'), '46-50'))] %>%
  ggplot(aes(insert, fD, color = exp_bin, lwd = exp_bin), alpha = 0.2) +
  #geom_point() +
  geom_line(stat = 'smooth', method = 'gam', formula = y ~ s(x, bs = 'cs'),
            alpha = 0.9, lineend = 'round') +
  scale_x_continuous(expand = expansion(add = 0.01), breaks = seq(0, 1, by = 0.25), labels = paste0(seq(100, 0, by = -25), 'th')) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_color_brewer(palette = 'Purples', name = '# of Spike ins') +
  scale_size_discrete(name = '# of Spike ins') +
  #scale_y_continuous(limits = c(0, 0.3)) +
  #scale_x_discrete(labels = c('<1', '1', '1.11', '1.25', '1.43', '1.67', '2', '2.5', '3.33', '5', '10')) +
  theme_classic(base_size = 20) +
  labs(x = 'Spike in specificity (percentile)', y = expression(paste(Delta, 'Contrast rank')))


panelB <- dope_scores %>%
  .[distance <= 1] %>%
  #.[is.na(stat), stat := 0] %>%
  .[, bin := as.factor(plyr::round_any(insert, 0.05))] %>%
  .[, exp_bin := factor(plyr::round_any(n_exp, 10), labels = c('1-5', mapply(paste, seq(6, 45, by = 10), seq(15, 50, by = 10), sep = '-'), '46-50'))] %>%
  ggplot(aes(insert, scoreD, color = exp_bin, lwd = exp_bin), alpha = 0.2) +
  #geom_point() +
  geom_line(stat = 'smooth', method = 'gam', formula = y ~ s(x, bs = 'cs'),
            alpha = 0.9, lineend = 'round') +
  scale_x_continuous(expand = expansion(add = 0.01), breaks = seq(0, 1, by = 0.25), labels = paste0(seq(100, 0, by = -25), 'th')) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_brewer(palette = 'Purples', name = '# of Spike ins') +
  scale_size_discrete(name = '# of Spike ins') +
  #scale_y_continuous(limits = c(0, 0.3)) +
  #scale_x_discrete(labels = c('<1', '1', '1.11', '1.25', '1.43', '1.67', '2', '2.5', '3.33', '5', '10')) +
  theme_classic(base_size = 20) +
  labs(x = 'Spike in specificity (percentile)', y = expression(paste(Delta, 'Test statistic')))

library(ggpubr)
ggarrange(panelA + labs(tag = 'A') + xlab('') + scale_x_continuous(expand = expansion(add = 0.01), breaks = seq(0, 1, by = 0.25), labels = rep('', 5)), panelB + labs(tag = 'B'), ncol = 1, common.legend = T, legend = 'right')

dope_scores %>%
  .[, bin := as.factor(plyr::round_any(insert, 0.05))] %>%
  .[distance <= 1 & bin %in% c(0, 0.05, 0.1, 0.9, 0.95, 1)] %>%
  #.[is.na(stat), stat := 0] %>%
  .[, bin := as.factor(plyr::round_any(insert, 0.05))] %>%
  #.[, exp_bin := factor(plyr::round_any(n_exp, 10), labels = c('1-5', mapply(paste, seq(6, 45, by = 10), seq(15, 50, by = 10), sep = '-'), '46-50'))] %>%
  ggplot(aes(n_exp, f, color = factor(bin, levels = c(0, 0.05, 0.1, 0.9, 0.95, 1), labels = c('100th', '95th', '90th', '10th', '5th', '0th'), ordered = T))) +
  stat_summary()
  stat_smooth(method = 'loess', se = F, lwd = 2) +
  scale_color_brewer(palette = 'Dark2', name = 'Spike in specificity (percentile)') +
  scale_y_continuous(breaks = c(0.1, 0.5), labels = c('Top 10%', 'Top 50%')) +
  theme_classic(base_size = 20) +
  labs(x = '# of Spike ins', y = 'Contrast rank')

dope_scores %>%
  .[is.na(stat), stat := 0] %>%
  .[, bin := as.numeric(as.character(as.factor(1/plyr::round_any(index / dropout, 0.1))))] %>%
  .[is.finite(bin)] %>%
  .[, bin := round(bin, 2)] %>%
  .[bin < 1, bin := 0.9] %>%
  #.[, .(index = mean(index, na.rm = T),
  #      stat = mean(stat, na.rm = T),
  #      stat_sd = sd(stat, na.rm = T),
  #      ranking = mean(ranking, na.rm = T),
  #      ranking_sd = sd(ranking, na.rm = T)), bin] %>%
  ggplot(aes(factor(bin), ranking/28327)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic(base_size = 20) +
  #geom_errorbar(aes(ymin = stat - stat_sd,
  #                  ymax = stat + stat_sd))
  scale_y_continuous(breaks = c(1, 0.5, 0), labels = c('End', 'Middle', 'Start')) +
  scale_x_discrete(labels = c('<1', '1', '1.11', '1.25', '1.43', '1.67', '2', '2.5', '3.33', '5', '10')) +
  labs(x = 'Experiment "Goodness"', y = 'Contrast Position')

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
