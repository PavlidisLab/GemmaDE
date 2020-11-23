lapply(sample(1:nrow(DATA.HOLDER$artificial@experiment.meta), 2000), function(i) {
  contrast <- DATA.HOLDER$artificial@experiment.meta[i, .(cf.BaseLongUri, cf.ValLongUri)]
  mData <- DATA.HOLDER$artificial@experiment.meta[cf.BaseLongUri == contrast$cf.BaseLongUri &
                                                    cf.ValLongUri == contrast$cf.ValLongUri]
  
  genes <- paste0('g', mData[, cf.GeneDrive]) %>% unique
  
  list(genes = genes,
       contrast = paste0(contrast$cf.BaseLongUri, ' vs. ', contrast$cf.ValLongUri),
       results = lapply(c('mvsm', 'zscore'), function(method) {
         options(app.search_method = method)
         options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                                        fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                                        geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff'),
                                        max.rows = getOption('max.rows'), method = getOption('app.search_method')))
         
         dt <- Sys.time()
         search(genes, 'artificial', verbose = F) -> searched
         timed <- Sys.time() - dt
         
         tmp <- which(searched$rn %in% mData$rsc.ID)
         tmp1 <- tmp / length(searched$rn)
         enrich(searched, 'artificial', verbose = F) -> enriched
         
         mA <- first(ONTOLOGIES.DEFS[Node_Long == contrast$cf.BaseLongUri, Definition]) %>% {
           switch((length(.) == 0) + 1, ., as.character(contrast$cf.BaseLongUri))
         }
         
         mB <- first(ONTOLOGIES.DEFS[Node_Long == contrast$cf.ValLongUri, Definition]) %>% {
           switch((length(.) == 0) + 1, ., as.character(contrast$cf.ValLongUri))
         }
         
         if(mB %in% CACHE.BACKGROUND$artificial$cf.BaseLongUri) {
           BB <- mA
           mA <- mB
           mB <- BB
         }
         
         list(method = method,
              indices = tmp,
              percs = tmp1,
              scores = searched$score,
              timediff = timed,
              pv = enriched[cf.BaseLongUri == mA & cf.ValLongUri == mB,
                            pv.fisher],
              otherpv = enriched[sample(1:.N, 1), pv.fisher])
       }))
}) -> bootstrap
saveRDS(bootstrap, '/space/scratch/jsicherman/Thesis Work/data/compareMethods.rds')
rm(bootstrap)

rbindlist(lapply(1:length(bootstrap), function(x) {
  data.table(run = x,
             mvsm = bootstrap[[x]]$results[[1]]$timediff / length(bootstrap[[x]]$genes),
             zscore = bootstrap[[x]]$results[[2]]$timediff / length(bootstrap[[x]]$genes))
})) -> data.to.plot

data.to.plot %>% melt(measure.vars = c('mvsm', 'zscore')) %>%
  .[, variable := ifelse(variable == 'mvsm', 'M-VSM', 'Weighted Score')] %>% as.data.frame() %>%
  ggplot(aes(1, log10(as.double(value)), color = variable, fill = variable)) +
  geom_violin(alpha = 0.1, position = 'identity', trim = F) + theme_cowplot(font_size = 19) +
  xlab(element_blank()) + ylab(expression(log[10]~time~per~gene)) +
  scale_color_discrete(name = 'Method') +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = c(0.82, 0.9)) +
  scale_fill_discrete(name = 'Method') + coord_flip() +
  ggtitle('Comparison of experiment ranking methods')

rbindlist(lapply(1:length(bootstrap), function(x) {
  data.table(run = x,
             mvsm = bootstrap[[x]]$results[[1]]$scores,
             zscore = bootstrap[[x]]$results[[2]]$scores)
})) -> data.to.plot

data.to.plot %>% .[sample(1:nrow(data.to.plot), min(nrow(data.to.plot), 100000)), ] %>%
  melt(id.vars = 'run') %>%
  ggplot(aes(value, fill = variable)) + geom_density()

data.to.plot %>% .[sample(1:nrow(data.to.plot), min(nrow(data.to.plot), 100000)), ] %>%
  ggplot(aes(mvsm, zscore)) + geom_point(alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, colour = 'red') +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 5)) + scale_y_continuous(expand = c(0, 0), limits = c(1, 15)) +
  theme_cowplot(font_size = 19) +
  ggtitle('Comparison of experiment ranking methods') +
  xlab('M-VSM') + ylab('Weighted Score')

rbindlist(lapply(1:length(bootstrap), function(x) {
  data.table(run = x,
             mvsm = bootstrap[[x]]$results[[1]]$percs,
             zscore = bootstrap[[x]]$results[[2]]$percs)
})) -> data.to.plot

data.to.plot %>% melt(measure.vars = c('mvsm', 'zscore')) %>%
  .[, variable := ifelse(variable == 'mvsm', 'M-VSM', 'Weighted Score')] %>% as.data.frame() %>%
  ggplot(aes(1, -log10(value), fill = variable, color = variable)) +
  geom_violin(alpha = 0.1, position = 'identity') + theme_cowplot(font_size = 19) +
  xlab(element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x.top = element_text(angle = 45, hjust = 0),
        legend.position = c(0.8, 0.9)) +
  scale_y_continuous(expand = c(0, 0), sec.axis =
                       sec_axis(~10^(-.),
                                name = expression(index[fraction]),
                                breaks = c(1, 0.75, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005),
                                labels = c('100%', '75%', '50%', '40%', '30%', '20%', '10%', '5%', '1%', '0.5%', '0.1%', '0.05%'))) +
  ylab(expression(-log[10]~index[fraction])) +
  coord_flip() +
  scale_color_discrete(name = 'Method') +
  scale_fill_discrete(name = 'Method') +
  ggtitle('Comparison of experiment ranking methods')

rbindlist(lapply(1:length(bootstrap), function(x) {
  data.table(run = x,
             mvsm = bootstrap[[x]]$results[[1]]$pv,
             zscore = bootstrap[[x]]$results[[2]]$pv)
})) -> data.to.plot

data.to.plot %>% melt(measure.vars = c('mvsm', 'zscore')) %>% .[!is.na(value)] %>%
  .[, variable := ifelse(variable == 'mvsm', 'M-VSM', 'Weighted Score')] %>% as.data.frame() %>%
  ggplot(aes(1, -log10(value), fill = variable, color = variable)) +
  geom_violin(alpha = 0.1, position = 'identity', trim = F) + theme_cowplot(font_size = 19) +
  xlab(element_blank()) + stat_summary(fun = median, shape = 4, size = 2) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        legend.position = c(0.8, 0.9)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab(expression(-log[10]~'p-value')) +
  scale_color_discrete(name = 'Method') +
  scale_fill_discrete(name = 'Method') +
  ggtitle('Comparison of experiment ranking methods')

###########

lapply(c('mvsm', 'zscore'), function(method) {
  options(app.search_method = method)
  options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                                 fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                                 geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff'),
                                 max.rows = getOption('max.rows'), min.tags = getOption('app.min.tags'),
                                 method = getOption('app.search_method')))
  
  DATA.HOLDER$human@gene.meta[gene.Name %in% c('KDM5D', 'XIST', 'RPS4Y1'), entrez.ID] -> IDs
  
  IDs %>% search('human', verbose = F) -> tmp
  assign(paste0('human.', method), enrich(tmp, 'human', verbose = F), envir = globalenv())
  
  tmp %>% merge(DATA.HOLDER$human@experiment.meta[, .(rsc.ID, cf.Baseline, cf.Val)],
                by.x = 'rn', by.y = 'rsc.ID', sort = F) %>%
    .[, isProper := grepl('(male|sex)', cf.Baseline) | grepl('male', cf.Val)] %>%
    .[, method := method]
}) -> human.plots

lapply(c('mvsm', 'zscore'), function(method) {
  options(app.search_method = method)
  options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                                 fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                                 geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff'),
                                 max.rows = getOption('max.rows'), min.tags = getOption('app.min.tags'),
                                 method = getOption('app.search_method')))
  
  DATA.HOLDER$mouse@gene.meta[gene.Name %in% c('Kdm5d', 'Xist'), entrez.ID] -> IDs
  
  IDs %>% search('mouse', verbose = F) -> tmp
  assign(paste0('mouse.', method), enrich(tmp, 'mouse', verbose = F), envir = globalenv())
  
  tmp %>% merge(DATA.HOLDER$mouse@experiment.meta[, .(rsc.ID, cf.Baseline, cf.Val)],
                by.x = 'rn', by.y = 'rsc.ID', sort = F) %>%
    .[, isProper := grepl('(male|sex)', cf.Baseline) | grepl('male', cf.Val)] %>%
    .[, method := method]
}) -> mouse.plots

rbindlist(mouse.plots) %>% .[, score, .(isProper, method)] %>%
  .[, isProper := ifelse(isProper, 'Male/female', 'Other')] %>%
  .[, method := ifelse(method == 'mvsm', 'M-VSM', 'Weighted Score')] %>%
  ggplot(aes(isProper, score, color = isProper)) + geom_jitter() +
  coord_flip() + stat_summary(fun = median, shape = 3, size = 2) +
  facet_wrap(~method, scales = 'free_x') + scale_y_continuous(expand = c(0, 0)) +
  theme_cowplot(font_size = 20) + theme(legend.position = 'none') +
  xlab(element_blank()) + ylab('Experiment Score')

(rbindlist(human.plots) %>% .[, species := 'Human'] %>% .[, .(score, isProper, method, species)] %>% rbind(
  rbindlist(mouse.plots) %>% .[, species := 'Mouse'] %>%  .[, .(score, isProper, method, species)] 
) %>% .[, method := ifelse(method == 'mvsm', 'M-VSM', 'Weighted Score')] %>%
  ggplot(aes(d = isProper, m = score, color = method)) + geom_roc(n.cuts = 0) + style_roc() +
  theme(text = element_text(size = 19)) +
  facet_wrap(~species) + scale_color_discrete(name = 'Method') +
  ggtitle('Male/female contrasts based on KDM5D, XIST and RPS4Y1')) %>% calc_auc()

lapply(c('mvsm', 'zscore'), function(method) {
  options(app.search_method = method)
  options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                                 fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                                 geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff'),
                                 max.rows = getOption('max.rows'), method = getOption('app.search_method')))
  
  DATA.HOLDER$mouse@gene.meta[gene.Name %in% c('Olig1', 'Olig2', 'Olig3'), entrez.ID] -> IDs
  
  IDs %>% search('mouse', verbose = F) -> tmp
  list(search = tmp,
       enrich = enrich(tmp, 'mouse', verbose = F))
}) -> mouse.olig2

lapply(c('mvsm', 'zscore'), function(method) {
  options(app.search_method = method)
  options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                                 fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                                 geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff'),
                                 max.rows = getOption('max.rows'), method = getOption('app.search_method')))
  
  DATA.HOLDER$mouse@gene.meta[gene.Name %in% c('Gfap', 'Slc1a2', 'Aldh1l1'), entrez.ID] -> IDs
  
  IDs %>% search('mouse', verbose = F) -> tmp
  list(search = tmp,
       enrich = enrich(tmp, 'mouse', verbose = F))
}) -> mouse.gfap
