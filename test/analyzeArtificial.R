analyzeArtificial <- function(N_GROUPS, N_GENES, COI = NULL, GOI = NULL, seed = NULL,
                              max.distance = Inf, plot = T) {
  if(!is.null(seed))
    set.seed(seed)
  
  if(is.null(COI))
    COI <- DATA.HOLDER$artificial@experiment.meta$cf.ID %>% unique %>% sample(N_GROUPS)
  if(is.null(GOI))
    GOI <- contrasts[contrast %in% COI] %>% setorder(-probability) %>% head(N_GENES) %>% .$entrez.ID
  SOI <- contrasts[contrast %in% COI] %>% setorder(-probability) %>% head(N_GENES) %>% .[, .N, contrast]
  
  if(nrow(SOI) < N_GROUPS)
    SOI <- rbind(SOI, rbindlist(lapply(setdiff(COI, SOI$contrast), function(x) data.table(contrast = x, N = 0))))
  
  tmp <- search(as.character(GOI), getConfig(taxa = 'artificial'))
  
  mSearchData <- tmp %>%
    merge(DATA.HOLDER$artificial@experiment.meta[, .(expected = cf.ID, rn = rsc.ID)], by = 'rn', sort = F) %>%
    merge(dcast(contrasts[entrez.ID %in% GOI, !'effect'],
                contrast ~ entrez.ID, value.var = 'probability'), # No g is probabilities
          by.x = 'expected', by.y = 'contrast', sort = F) %>%
    .[, probability := rowSums(.SD[, -(1:which(names(.) == 'ee.q'))])] %>%
    .[, I := .I] %>%
    .[, group := factor(expected,
                        ordered = T,
                        levels = COI) %>% {
                          a <- as.character(.)
                          a[is.na(a)] <- 'Other'
                          factor(a, levels = c(levels(.), 'Other'), ordered = T)
                        }]
  
  tmp2 <- enrich(tmp, getConfig(taxa = 'artificial')) %>%
    .[distance <= max.distance]
  
  tmp3 <- DATA.HOLDER$human@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
    unique %>% .[COI] %>%
    .[, cf.Cat := as.character(cf.Cat)] %>%
    .[, cf.BaseLongUri := as.character(cf.BaseLongUri)] %>%
    .[, cf.ValLongUri := as.character(cf.ValLongUri)]

  tmp4 <- lapply(1:nrow(tmp3), function(i) {
    list(baseline = ONTOLOGIES.DEFS[grepl(
      paste0('(', paste0(strsplit(tmp3$cf.BaseLongUri[i], '; ', T)[[1]], collapse = '|'), ')'), Node_Long),
      as.character(Definition)],
      contrast = ONTOLOGIES.DEFS[grepl(
        paste0('(', paste0(strsplit(tmp3$cf.ValLongUri[i], '; ', T)[[1]], collapse = '|'), ')'), Node_Long),
        as.character(Definition)])
  }) %>% `names<-`(COI)
  
  mEnrichData <- tmp2 %>% .[, I := .I] %>%
    merge(
      merge(CACHE.BACKGROUND$artificial[, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
            DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
            by = 'rsc.ID') %>% .[, !'rsc.ID'],
      by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri')
    )
  
  mEnrichData2 <- tmp2 %>% .[, I := .I] %>%
    merge(
      merge(CACHE.BACKGROUND$artificial[, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
            DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
            by = 'rsc.ID') %>% .[, !'rsc.ID'],
      by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri')
    )
  
  print(COI)
  print(SOI)
  print(GOI)
  
  mTableData <- data.frame(Group = COI,
                           `N_Genes` = SOI[match(COI, SOI$contrast), N],
                           `N_Inf` = merge(CACHE.BACKGROUND$artificial[distance <= max.distance, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                                           DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
                                           by = 'rsc.ID') %>%
                             .[, .(rsc.ID, cf.ID)] %>% unique %>%
                             .[cf.ID %in% COI, .N, cf.ID] %>%
                             .[match(COI, cf.ID), N],
                           `N_Exp` = merge(CACHE.BACKGROUND$artificial[distance == 0, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                                           DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
                                           by = 'rsc.ID') %>%
                             .[, .(rsc.ID, cf.ID)] %>% unique %>%
                             .[cf.ID %in% COI, .N, cf.ID] %>%
                             .[match(COI, cf.ID), N],
                           `Mean_DEGs` = merge(CACHE.BACKGROUND$artificial[distance == 0, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                                               DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, n.DE, cf.ID)],
                                               by = 'rsc.ID') %>%
                             .[, .(rsc.ID, n.DE, cf.ID)] %>% unique %>%
                             .[cf.ID %in% COI, .(D = round(sum(n.DE) / .N, 2)), cf.ID] %>%
                             .[match(COI, cf.ID), D])
  
  if(plot) {
    mSearchData %>%
      ggplot(aes(probability, score, color = group)) + geom_point() +
      theme_classic(base_size = 20) + scale_color_brewer(palette = 'Dark2') +
      labs(x = 'Total weight', y = 'Score', title = 'Search Results', colour = 'Group', tag = 'A') -> searchPlot
    
    mEnrichData %>%
      .[cf.ID %in% COI] %>%
      .[, cf.ID := factor(cf.ID, levels = COI, ordered = T)] %>%
      ggplot(aes(cf.ID, stat, color = cf.ID)) + geom_jitter() +
      stat_summary(fun = median, color = 'black', shape = 3) +
      theme_classic(base_size = 20) + scale_color_brewer(palette = 'Dark2') +
      labs(x = 'Group', y = 'SD from mean', title = 'Enrich Results', tag = 'B') +
      theme(legend.position = 'none') +
      coord_flip() -> enrichPlot
    
    mEnrichData2 %>%
      .[, f := I / max(I)] %>%
      .[cf.ID %in% COI] %>%
      .[, cf.ID := factor(cf.ID, levels = COI, ordered = T)] %>%
      ggplot(aes(cf.ID, f, color = cf.ID)) + geom_boxplot() +
      theme_classic(base_size = 20) + scale_color_brewer(palette = 'Dark2') +
      scale_y_continuous(labels = scales::percent) +
      labs(x = 'Group', y = 'Fractional index', tag = 'C') +
      theme(legend.position = 'none') +
      coord_flip() -> enrichPlot2
    
    mTable <- tableGrob(mTableData, theme = ttheme_minimal(base_size = 17))
    mContrast <- tableGrob(data.frame(
      Group = COI,
      Contrast = merge(CACHE.BACKGROUND$artificial[distance == 0, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                       DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
                       by = 'rsc.ID') %>%
        .[cf.ID %in% COI, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, cf.ID)] %>% unique %>%
        .[, .(P = paste0(cf.Cat, ': ', cf.BaseLongUri, ' vs. ', cf.ValLongUri)), cf.ID] %>%
        .[match(COI, cf.ID), P]
    ), theme = ttheme_minimal(base_size = 10))
    
    grid.arrange(grobs = list(searchPlot, enrichPlot, enrichPlot2, mTable, mContrast),
                 layout_matrix = rbind(c(1, 2),
                                       c(1, 2),
                                       c(4, 3),
                                       c(5, 3)))
  }
  
  invisible(list(search = mSearchData[, .(expected, rn, score, probability, I, group)],
                 enrich = mEnrichData %>% unique,
                 diagnostic = mTableData,
                 genes = GOI))
}

lapply(1:10, function(groups) {
  lapply(1:10, function(genes) {
    lapply(1, function(rep) {
      tmp <- analyzeArtificial(groups, genes, plot = F)
      list(search = tmp$search[group %in% tmp$diagnostic$Group],
           enrich = tmp$enrich[cf.ID %in% tmp$diagnostic$Group],
           diagnostic = tmp$diagnostic,
           genes = tmp$genes)
    })
  })
}) -> tmp

lapply(1:length(tmp), function(groups) {
  lapply(1:length(tmp[[groups]]), function(genes) {
    lapply(1:length(tmp[[groups]][[genes]]), function(rep) {
      data.table(groups = groups,
                 genes = genes,
                 rep = rep,
                 tmp[[groups]][[genes]][[rep]]$search %>%
                   .[, group := as.integer(as.character(group))] %>% merge(
                   tmp[[groups]][[genes]][[rep]]$diagnostic,
                   by.x = 'group', by.y = 'Group'
                 ) %>% .[, .(probability, I, score, N_Genes, N_Exp, Mean_DEGs)])
    }) %>% rbindlist
  }) %>% rbindlist
}) %>% rbindlist %>%
  .[, groups := factor(paste0(groups, ' Groups'), ordered = T, levels = paste0(1:length(tmp), ' Groups'))] %>%
  ggplot(aes(factor(N_Genes), score)) +
  geom_boxplot() + facet_wrap(~groups) +
  theme_classic() +
  labs(x = 'Number of genes', y = 'Search score')

lapply(1:length(tmp), function(groups) {
  lapply(1:length(tmp[[groups]]), function(genes) {
    lapply(1:length(tmp[[groups]][[genes]]), function(rep) {
      data.table(groups = groups,
                 genes = genes,
                 rep = rep,
                 tmp[[groups]][[genes]][[rep]]$enrich[distance == 0, .(stat, I)])
    }) %>% rbindlist
  }) %>% rbindlist
}) %>% rbindlist %>% ggplot(aes(factor(groups), I)) + geom_boxplot() + facet_wrap(~genes)
