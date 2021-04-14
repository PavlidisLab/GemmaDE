source('/home/jsicherman/Thesis Work/requirements.R')

source('dependencies.R')

library(lhs)
library(parallel)

DATA.HOLDER[c('human', 'mouse', 'rat')] <- NULL
CACHE.BACKGROUND[c('human', 'mouse', 'rat')] <- NULL
NULLS.EXP[c('human', 'mouse', 'rat')] <- NULL

analyzeArtificial <- function(N_GROUPS, N_GENES, COI = NULL, GOI = NULL, seed = NULL,
                              max.distance = Inf, best.index = 1) {
  if(!is.null(seed))
    set.seed(seed)
  
  if(is.null(COI))
    COI <- DATA.HOLDER$artificial@experiment.meta$cf.ID %>% unique %>% sample(N_GROUPS)
  if(is.null(GOI)) {
    mContrasts <- contrasts[contrast %in% COI] %>% setorder(-probability)
    
    if(best.index + N_GENES > nrow(mContrasts))
      best.index <- nrow(mContrasts) - N_GENES
    
    GOI <- mContrasts[best.index:(best.index + N_GENES - 1), entrez.ID]
  }
  
  tmp <- search(as.character(GOI), getConfig(taxa = 'artificial'))
  
  mSearchData <- tmp %>%
    merge(DATA.HOLDER$artificial@experiment.meta[, .(expected = cf.ID, rn = rsc.ID)], by = 'rn', sort = F) %>%
    merge(dcast(contrasts[entrez.ID %in% GOI, !'effect'],
                contrast ~ entrez.ID, value.var = 'probability'), # No g is probabilities
          by.x = 'expected', by.y = 'contrast', sort = F) %>%
    .[, probability := rowSums(.SD[, -(1:which(names(.) == 'ee.q'))])] %>%
    .[, I := .I] %>%
    .[, .(rn, expected, probability, I)]
  
  tmp2 <- enrich(tmp, getConfig(taxa = 'artificial')) %>%
    .[distance <= max.distance] %>% .[, I := .I]
  
  mEnrichData <- tmp2 %>% merge(
    merge(CACHE.BACKGROUND$artificial[, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
          DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
          by = 'rsc.ID', sort = F) %>% .[, !'rsc.ID'] %>% unique,
    by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), sort = F
  ) %>% merge(
    mSearchData[, .(experiment_count = .N, probability = unique(probability)), expected],
    by.x = 'cf.ID', by.y = 'expected', sort = F
  ) %>% .[, .(distance = mean(distance), A = mean(A), B = mean(B), C = mean(C),
              stat = max(stat), I = mean(I), experiment_count = sum(experiment_count),
              probability = mean(probability)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
  
  list(search = mSearchData,
       enrich = mEnrichData %>% unique,
       genes = GOI,
       contrasts = COI,
       best.index = best.index)
}

if(!exists('hypercube')) {
  print('Making hypercube')
  hypercube <- improvedLHS(5000, 3)
  hypercube[, 1] <- 1L + as.integer(hypercube[, 1] * rexp(nrow(hypercube), 0.8)) # Groups
  hypercube[, 2] <- 1L + as.integer(hypercube[, 2] * 99) # Genes
  hypercube[, 3] <- 1L + as.integer(hypercube[, 3] * 500) # Best index
}

contrasts <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/contrast_aff.rds')

options(mc.cores = 3)
mclapply(1:nrow(hypercube), function(iter) {
  system(paste0('echo ', iter))
  analyzeArtificial(hypercube[iter, 1],
                    hypercube[iter, 2],
                    best.index = hypercube[iter, 3])
}) %>% saveRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/bootstrapped_artificial2.rds')

if(F) {
  art_boot <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/bootstrapped_artificial.rds')
  tmp2 <- lapply(1:length(tmp), function(i) {
    data.table(groups = nrow(tmp[[i]]$diagnostic),
               genes = length(tmp[[i]]$genes),
               insert_index = tmp[[i]]$mean_index,
               insert_sd = tmp[[i]]$mean_sd,
               tmp[[i]]$search %>% merge(
                 tmp[[i]]$diagnostic,
                 by.x = 'group', by.y = 'Group'
               ) %>% .[, .(probability, I, score, N_Genes, N_Exp, Mean_DEGs)])
  }) %>% rbindlist
  tmp3 <- lapply(1:length(tmp), function(i) {
    data.table(groups = nrow(tmp[[i]]$diagnostic),
               genes = length(tmp[[i]]$genes),
               insert_index = tmp[[i]]$mean_index,
               insert_sd = tmp[[i]]$mean_sd,
               tmp[[i]]$enrich[distance == 0, .(cf.ID, stat, I, f)] %>% merge(
                 tmp[[i]]$diagnostic,
                 by.x = 'cf.ID', by.y = 'Group'
               ),
               tmp[[i]]$search[, .(probability = mean(probability),
                                   I = median(I)), group] %>% merge(
                 tmp[[i]]$diagnostic,
                 by.x = 'group', by.y = 'Group'
               ) %>% .[, .(probability, I_search = I)])
  }) %>% rbindlist

  tmp2 %>%
    ggplot(aes(factor(insert_index), score)) +
    geom_boxplot() +
    #facet_wrap(~N_Exp) +
    theme_classic(base_size = 20) + theme(legend.position = 'none') +
    labs(x = 'Mean insertion index', y = 'Search score')
  
  tmp3 %>% na.omit %>%
    ggplot(aes(insert_index, stat)) +
    geom_point() +
    theme_classic(base_size = 20) +
    geom_function(fun = function(x) 38 - 40*(x - 0.05)^(1/17), color = 'red') +
    theme(legend.position = 'none') +
    labs(x = 'Mean insertion index', y = 'Test statistic')
  
  tmp3 %>%
    ggplot(aes(factor(round(insert_index, 1)), stat)) +
    geom_boxplot() +
    theme_classic(base_size = 20) +
    labs(x = 'Mean insertion index', y = 'Test statistic')
  
  plotA <- tmp3 %>% copy %>%
    .[, x := factor(plyr::round_any(100 * (1-insert_index), 5))] %>%
    .[, x := paste0(x, '^{"th"}')] %>%
    .[, x := factor(x, levels = paste0(seq(0, 100, by = 5), '^{"th"}'), ordered = T)] %>%
    ggplot(aes(x, f)) +
    geom_boxplot(fill = '#b3bcc7') +
    theme_bw(base_size = 20) +
    theme(panel.grid = element_blank()) +
    scale_x_discrete(breaks = paste0(seq(0, 100, by = 10), '^{"th"}'), labels = rlang::parse_exprs) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c('First', 'Mid', 'Last')) +
    labs(x = 'Mean gene percentile',
         y = 'Rank of contrast', tag = 'A')
  
  plotB <- tmp3 %>% copy %>%
    .[, ppg := probability / max(probability) / genes] %>%
    .[, ppg := plyr::round_any((ppg - min(ppg)) / (max(ppg) - min(ppg)), 0.05)] %>%
    ggplot(aes(factor(ppg), f)) +
    geom_boxplot(fill = '#b3bcc7') +
    theme_bw(base_size = 20) +
    scale_x_discrete(breaks = seq(0, 1, by = 0.1), labels = function(x) {
      scales::percent(as.numeric(x), accuracy = 1) %>%
        gsub('^0%', '< 5%', .)
    }) +
    scale_y_continuous(breaks = NULL) +
    theme(panel.grid = element_blank()) +
    labs(x = 'Average probability of differential expression (per gene)',
         y = element_blank(), tag = 'B')
  
  gridExtra::grid.arrange(plotA, plotB, ncol = 2)
  
  tmp3 %>% na.omit %>%
    ggplot(aes(insert_index, -log10(f))) +
    geom_point(alpha = 0.2) +
    geom_function(fun = function(x) 17 - 17*(x - 0.05)^(1/41), color = 'red', size = 2) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = expansion(add = c(0, 0.1))) +
    theme_classic(base_size = 20) +
    labs(x = 'Mean insertion index', y = 'Enrichment fraction')
}
