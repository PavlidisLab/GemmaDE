source('/home/jsicherman/Thesis Work/requirements.R')

source('dependencies.R')

library(lhs)
library(parallel)

DATA.HOLDER[c('human', 'mouse', 'rat')] <- NULL
CACHE.BACKGROUND[c('human', 'mouse', 'rat')] <- NULL
NULLS[c('human', 'mouse', 'rat')] <- NULL

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
  ) %>% .[, .(distance = mean(distance), score = max(score),
              stat = max(stat), I = mean(I), experiment_count = sum(experiment_count),
              probability = mean(probability)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
  
  list(search = mSearchData,
       enrich = mEnrichData %>% unique,
       genes = GOI,
       contrasts = COI,
       best.index = best.index)
}

if(Sys.getenv('RSTUDIO') == '1') {

  art_boot <- readRDS(paste(DATADIR, 'artificial/bootstrapped_artificial2.rds', sep='/'))

  
  lapply(art_boot, '[[', 'enrich') %>%
    sapply(function(x) cor(x$probability, -x$I, use = 'complete')) %>%
      data.frame(mean = mean(.)) %>%
    ggplot(aes(.)) + geom_histogram(fill = 'white', color = 'black') +
    geom_vline(aes(xintercept = mean), lwd = 2, color = 'red') +
    theme_classic(base_size = 20) +
    scale_x_continuous(expand = c(0, 0), name = 'Pearson Correlation') +
    scale_y_continuous(expand = c(0, 0), name = 'Count')
  
  contrasts <- CACHE.BACKGROUND$artificial %>%
    merge(DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)], by = 'rsc.ID')
  
  lapply(art_boot, function(x) {
    accepting <- contrasts[distance == 0 & cf.ID %in% x$contrasts, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
      unique %>% .[, paste0(as.character(cf.Cat), as.character(cf.BaseLongUri), as.character(cf.ValLongUri))]
    list(best.index = x$best.index,
         enrich = x$enrich[paste0(cf.Cat, cf.BaseLongUri, cf.ValLongUri) %in% accepting],
         contrasts = x$contrasts)
  }) -> art_boot2
  
  lapply(art_boot, function(x) {
    accepting <- contrasts[distance != 0 & cf.ID %in% x$contrasts, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
      unique %>% .[, paste0(as.character(cf.Cat), as.character(cf.BaseLongUri), as.character(cf.ValLongUri))]
    list(best.index = x$best.index,
         enrich = x$enrich[paste0(cf.Cat, cf.BaseLongUri, cf.ValLongUri) %in% accepting],
         contrasts = x$contrasts)
  }) -> art_boot3
  
  lapply(1:length(art_boot), function(i) {
    if(nrow(art_boot2[[i]]$enrich) != 0)
      head(sort(art_boot3[[i]]$enrich[, I]), length(art_boot3[[i]]$contrasts))
  }) %>% unlist %>% data.table(rep = 1:length(.)) %>% cbind(type = 'Inferred') %>% rbind(
    lapply(art_boot2, function(x) {
      if(nrow(x$enrich) != 0)
        head(sort(x$enrich[, I]), length(x$contrasts))
    }) %>% unlist %>% data.table(rep = 1:length(.)) %>% cbind(type = 'Exact')
  ) %>% .[, .N, .(., type)] %>% .[, f := N / sum(N)] %>%
    .[. %in% 1:10] %>% ggplot(aes(., f, fill = type)) + geom_bar(stat = 'identity', color = 'white') +
    scale_x_continuous(expand = c(0, 0), name = 'Rank', breaks = 1:20) +
    scale_fill_brewer(palette = 'Dark2', name = element_blank()) +
    scale_y_continuous(expand = c(0, 0), name = 'Percent of searches', labels = scales::percent) +
    theme_classic(base_size = 20) + theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.key.size = unit(30, 'points'))
    
  
  plotA <- lapply(art_boot2, function(x) {
    if(nrow(x$enrich) != 0)
      data.table(i = x$best.index, j = head(sort(x$enrich[, I]), length(x$contrasts)))
  }) %>% rbindlist %>% { print(mean(.$i)); .$j } %>% data.table(rep = 1:length(.)) %>% .[, .N, .] %>% setorder(.) %>% {
    print(head(.$N/sum(.$N), 10) %>% sum)
    .
  } %>%
    ggplot(aes(., N/sum(N))) + geom_line() +
    scale_x_continuous(expand = c(0, 0), limits = c(1, 10), name = c('Computed Rank'), breaks = 1:10) +
    scale_y_continuous(expand = c(0, 0), labels = scales::percent, name = 'Percent of searches') +
    theme_classic(base_size = 20)
  
  plotB <- lapply(1:length(art_boot), function(i) {
    if(nrow(art_boot2[[i]]$enrich) != 0)
      head(sort(art_boot3[[i]]$enrich[, I]), length(art_boot3[[i]]$contrasts))
  }) %>% unlist %>% data.table(rep = 1:length(.)) %>% .[, .N, .] %>% setorder(.) %>% {
    print(head(.$N/sum(.$N), 10) %>% sum)
    .
  } %>%
    ggplot(aes(., N/sum(N))) + geom_line() +
    scale_x_continuous(expand = c(0, 0), limits = c(1, 10), name = c('Computed Rank'), breaks = 1:10) +
    scale_y_continuous(expand = c(0, 0), labels = scales::percent, name = 'Percent of searches') +
    theme_classic(base_size = 20)
  
  figure <- ggarrange(plotA + rremove("ylab") + rremove("xlab"),
                      plotB + rremove("ylab") + rremove("xlab"), # remove axis labels from plots
                      labels = LETTERS[1:2],
                      align = "hv", 
                      font.label = list(size = 20, color = "black", face = "bold", family = NULL, position = "top"))
  
  annotate_figure(figure, left = textGrob('Percent of rankings', rot = 90, vjust = 1, gp = gpar(fontsize = 20)),
                  bottom = textGrob('Computed rank', gp = gpar(fontsize = 20)))
} else {
  if(!exists('hypercube')) {
    print('Making hypercube')
    hypercube <- improvedLHS(5000, 3)
    hypercube[, 1] <- 1L + as.integer(hypercube[, 1] * rexp(nrow(hypercube), 0.8)) # Groups
    hypercube[, 2] <- 1L + as.integer(hypercube[, 2] * 99) # Genes
    hypercube[, 3] <- 1L + as.integer(hypercube[, 3] * 500) # Best index
  }
  
  if(!exists('contrasts'))

    contrasts <- readRDS(paste(DATADIR, 'artificial/contrast_aff.rds', sep='/'))

  
  options(mc.cores = 3)
  mclapply(1:nrow(hypercube), function(iter) {
    system(paste0('echo ', iter))
    analyzeArtificial(hypercube[iter, 1],
                      hypercube[iter, 2],
                      best.index = hypercube[iter, 3])

  }) %>% saveRDS(paste(DATADIR, 'artificial/bootstrapped_artificial2.rds', sep='/'))

}
