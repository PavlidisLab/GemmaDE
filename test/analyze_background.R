library(parallel)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

set.seed(18232)
options('mc.cores' = 15)

N.replicates <- 10000
N.sample <- 100

lapply(unique(ONTOLOGIES$OntologyScope), function(ontology) {
  boot <- mclapply(1:N.replicates, function(n) {
    experiments <- search(sample(DATA.HOLDER$human@gene.meta$entrez.ID, N.sample))
    tags <- enrich(experiments, scope = ontology)
  
    list(experiments = experiments, tags = tags)
  })
  
  list(
    terms = do.call(rbind, mclapply(boot, function(iteration) {
      if(is.list(iteration))
        data.table(tag = iteration$tags$Definition,
                   score = iteration$tags$`ranked N`)
      })) %>% as.data.table %>% .[, c('score.mean', 'score.median', 'score.min', 'score.max', 'score.var') :=
                                    list(mean(score, na.rm = T), median(score, na.rm = T), min(score, na.rm = T),
                                         max(score, na.rm = T), var(score, na.rm = T)), tag] %>% .[, score := NULL] %>% unique,
    experiments = do.call(rbind, mclapply(boot, function(iteration) {
      if(is.list(iteration))
        data.table(experiment = names(iteration$experiments),
                   score = iteration$experiments)
      })) %>% as.data.table %>% .[, c('score.mean', 'score.median', 'score.min', 'score.max', 'score.var') :=
                                    list(mean(score, na.rm = T), median(score, na.rm = T), min(score, na.rm = T),
                                         max(score, na.rm = T), var(score, na.rm = T)), experiment] %>% .[, score := NULL] %>% unique)
}) -> SUMMARY
names(SUMMARY) <- unique(ONTOLOGIES$OntologyScope)

saveRDS(SUMMARY, 'data/bootstrap_summary.rds')

# all.experiments <- do.call(rbind, mclapply(bootstrapped, function(entry) {
#   data.frame(experiment = names(entry$experiments),
#              score = entry$experiments)
# })) %>% as.data.table
# 
# all.experiments.binned.scored <- transform(all.experiments, group = cut(score, breaks = seq(0, 1, 0.1)))
# all.experiments.summarized.scored <- all.experiments.binned.scored.post[, experiment, group]
# 
# all.experiments.summarized.N <- all.experiments[, .N, experiment]
# all.experiments.binned.N <- transform(all.experiments.summarized.N,
#                                       group = cut(N, breaks = c(0, 10, 100, 250, 500, 750, 1000)))
# 
# all.experiments.binned.N <- all.experiments.binned.N %>%
#   merge(DATA.HOLDER$human@experiment.meta[, .(rsc.ID, n.DE)], by.x = 'experiment', by.y = 'rsc.ID')
# 
# all.experiments.binned.N[, .N, group] %>% ggplot(aes(group, N)) +
#   geom_bar(stat = 'identity', fill = 'white', color = 'black') + xlab('# of Times Hit') + ylab('# of Experiments') +
#   ggtitle(paste0('Distribution of Experiments (', N.replicates, ' replicates of ', N.sample, ')')) +
#   scale_y_continuous(expand = c(0, 0))
# 
# all.experiments.binned.scored[, .N, .(experiment, group)] %>% .[, sum(N), group] %>% ggplot(aes(group, V1)) +
#   geom_bar(stat = 'identity', fill = 'white', color = 'black') + xlab('Score') + ylab('# of Experiments') +
#   ggtitle(paste0('Distribution of Experiments (', N.replicates, ' replicates of ', N.sample, ')')) +
#   scale_y_continuous(expand = c(0, 0))
