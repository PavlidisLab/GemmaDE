devtools::load_all()
source(here::here("main/requirements.R"))
source(here::here("main/dependencies.R"))
# GEO holdings ----

experiments = readRDS('report/experiments.rds')

GEO <- fread('generate//GEO_holdings_2021.csv')
GEO %>% melt(measure.vars = 'Count') %>%
  rbind(GEO %>% melt(measure.vars = 'Count') %>%
          .[Organism %in% c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'), .(value = sum(value)), .(Year, Type)] %>% {
            data.table(Year = .$Year,
                       Type = .$Type,
                       Organism = 'Mammalian',
                       variable = 'Count',
                       value = .$value)
          }) %>%
  .[, .(value = sum(value)), .(Year, Organism)] %>% {
    data.table(Year = .$Year,
               Type = 'Aggregate',
               Organism = .$Organism,
               variable = 'Count',
               value = .$value)
  } %>%
  rbind(GEO %>% melt(measure.vars = 'Count') %>%
          rbind(GEO %>% melt(measure.vars = 'Count') %>%
                  .[Organism %in% c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus'), .(value = sum(value)), .(Year, Type)] %>% {
                    data.table(Year = .$Year,
                               Type = .$Type,
                               Organism = 'Mammalian',
                               variable = 'Count',
                               value = .$value)
                  })) %>%
  #.[Organism == 'Mammalian'] %>%
  #.[, ingemma := experiments[taxon.Name %in% c('human', 'mouse', 'rat'), length(unique(ee.ID))]] %>%
  .[, Type := factor(Type, levels = c('Microarray', 'Sequencing', 'Aggregate'), labels = c('Microarray', 'Sequencing', 'Total'), ordered = T)] %>%
  .[Organism %in% c('Homo sapiens', 'Mus musculus', 'Rattus norvegicus')] %>%
  .[, ingemma := case_when(Organism == 'Homo sapiens' ~ experiments[experiment.Database == 'GEO' & taxon.Name == 'human' &
                                                                      (Type == 'Total' | (platform.Type %in% c('GENELIST', 'SEQUENCING')) == (Type == 'Sequencing')), length(unique(experiment.ID))],
                        Organism == 'Mus musculus' ~ experiments[experiment.Database == 'GEO' & taxon.Name == 'mouse' &
                                                                   (Type == 'Total' | (platform.Type %in% c('GENELIST', 'SEQUENCING')) == (Type == 'Sequencing')), length(unique(experiment.ID))],
                        Organism == 'Rattus norvegicus' ~ experiments[experiment.Database == 'GEO' & taxon.Name == 'rat' &
                                                                        (Type == 'Total' | (platform.Type %in% c('GENELIST', 'SEQUENCING')) == (Type == 'Sequencing')), length(unique(experiment.ID))]), Type] %>%
  ggplot(aes(Year, value, color = Type)) +
  geom_point() + geom_line(lwd = 1.1) +
  #geom_hline(aes(yintercept = ingemma), color = 'black', linetype = 'dashed', lwd = 1.5) +
  geom_hline(aes(color = Type, yintercept = ingemma), lwd = 1.4, linetype = 'dashed') +
  facet_wrap(~Organism, scales = 'free_y', ncol = 1) +
  theme_bw(20) + xlab(element_blank()) + ylab('# of Data Series') +
  theme(strip.background = element_blank(), strip.text = element_text(face = 'bold.italic', hjust = 0)) +
  scale_color_brewer(palette = 'Dark2') +
  scale_x_continuous(expand = c(0, 0))

# Number of data series exceeds half million in 20 years
GEO %>%
  mutate(YearQuarter = Year + (Quarter - 1) / 4,
         RootSeries = sqrt(2/3 * Series)) %>% {
           lm(RootSeries ~ YearQuarter, .)
         } %>% predict(data.frame(YearQuarter = 2022:2050)) %>% `^`(2)

DATA.HOLDER$human@experiment.meta %>% nrow
DATA.HOLDER$human@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique %>% nrow

DATA.HOLDER$mouse@experiment.meta %>% nrow
DATA.HOLDER$mouse@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique %>% nrow

DATA.HOLDER$rat@experiment.meta %>% nrow
DATA.HOLDER$rat@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique %>% nrow

rbindlist(list(DATA.HOLDER$human@experiment.meta,
               DATA.HOLDER$mouse@experiment.meta,
               DATA.HOLDER$rat@experiment.meta))[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique %>% nrow

CACHE.BACKGROUND$human %>% nrow
CACHE.BACKGROUND$human[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique %>% nrow

CACHE.BACKGROUND$mouse %>% nrow
CACHE.BACKGROUND$mouse[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique %>% nrow

CACHE.BACKGROUND$rat %>% nrow
CACHE.BACKGROUND$rat[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique %>% nrow

rbindlist(list(CACHE.BACKGROUND$human,
               CACHE.BACKGROUND$mouse,
               CACHE.BACKGROUND$rat))[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique %>% nrow

# Overlaps ----
rbindlist(list(DATA.HOLDER$human@experiment.meta[, .(N = length(unique(ee.ID)), species = 'Human'),
                                                 .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
               DATA.HOLDER$rat@experiment.meta[, .(N = length(unique(ee.ID)), species = 'Rat'),
                                               .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
               DATA.HOLDER$mouse@experiment.meta[, .(N = length(unique(ee.ID)), species = 'Mouse'),
                                                 .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)])) %>%
  .[, whichOne := 'Pre-inference'] %>%
  .[, mean := mean(N), species] %>%
  .[, median := median(N), species] %>%
  rbind(
    rbindlist(list(CACHE.BACKGROUND$human[, .(N = length(unique(ee.ID)), species = 'Human'),
                                          .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                   CACHE.BACKGROUND$rat[, .(N = length(unique(ee.ID)), species = 'Rat'),
                                        .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                   CACHE.BACKGROUND$mouse[, .(N = length(unique(ee.ID)), species = 'Mouse'),
                                          .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)])) %>%
      .[, whichOne := 'Post-inference'] %>%
      .[, mean := mean(N), species] %>%
      .[, median := median(N), species]
  ) %>%
  .[, whichOne := factor(whichOne, levels = c('Pre-inference', 'Post-inference'))] %>% {
    print(.[, .(whichOne, mean, median)] %>% unique)
    .
  } %>%
  ggplot(aes(log10(N), fill = species)) +
  geom_density() +
  geom_vline(aes(xintercept = log10(mean)), color = 'red') +
  geom_vline(aes(xintercept = log10(median)), color = 'black') +
  facet_grid(vars(whichOne), vars(species), scales = 'free_y', switch = 'y') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2.25), breaks = c(0, log10(c(2, 5)), 1, log10(c(25, 50)), 2), labels = function(x) round(10^x)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = 'Dark2') +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none',
        text = element_text(colour = 'black'), panel.grid = element_blank(),
        strip.placement = 'outside',
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = 'bold')) +
  xlab(expression("#"~Overlapping~Experiments)) +
  ylab('Density')

# Fraction of ontology terms ----
rbindlist(list(DATA.HOLDER$human@experiment.meta,
               DATA.HOLDER$mouse@experiment.meta,
               DATA.HOLDER$rat@experiment.meta)) %>% {
                 n <- nrow(.[grepl('http', cf.BaseLongUri) | grepl('http', cf.ValLongUri)])
                 n / nrow(.)
               }

# Tag summary ----
rbindlist(lapply(names(CACHE.BACKGROUND.FULL), function(i) data.table(taxon = i, CACHE.BACKGROUND.FULL[[i]]))) %>%
  .[, rsc.ID := as.integer(as.factor(rsc.ID))] -> mCache

lapply(c('nervous', 'reproductive', 'digestive', 'respiratory', 'hemolymphoid', 'endocrine', 'exocrine', 'cardiovascular', 'hepatobilliary',
         'brain', 'blood', 'liver', 'lung', 'muscle', 'intestine', 'heart', 'spleen', 'kidney', 'spinal cord',
         'hematopoietic cell', 'leukocyte', 'epithelial cell', 'glial cell', 'fibroblast', 'macrophage', 'T cell', 'B cell', 'embryonic stem cell', 'microglial cell'),
       function(topic) {
         data.table(Topic = topic, mCache[grepl(topic, cf.BaseLongUri, T) | grepl(topic, cf.ValLongUri, T), .(EE = length(unique(ee.ID)), RSC = length(unique(rsc.ID))), taxon])
       }) %>% rbindlist -> tagSummary

tagSummary %>%
  .[, subject := c(rep('System', 25), rep('Organ/Tissue', 30), rep('Cell Type', 29))] %>%
  melt(measure.vars = c('EE', 'RSC')) %>%
  .[!is.na(value) & variable == 'RSC'] %>%
  .[, stat := sum(value), Topic] %>%
  setorder(-stat) %>%
  .[, Topic := factor(Topic, levels = unique(Topic), ordered = T)] %>%
  ggplot(aes(Topic, value, fill = taxon)) +
  geom_bar(stat = 'identity') + facet_wrap(~subject, scales = 'free') + coord_flip() +
  theme_classic(20) + scale_y_continuous(expand = c(0, 0)) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = 'bold'), legend.text.align = 0, legend.justification = 'top') +
  xlab(element_blank()) + ylab('Condition Comparisons') +
  scale_fill_brewer(palette = 'Dark2', name = 'Taxon', labels = c(expression(italic('H. sapiens')), expression(italic('M. musculus')), expression(italic('R. norvegicus'))))

# Purple bar for astro genes ----
DATA.HOLDER$mouse@gene.meta[gene.Name %in% strsplit('Cyp4f15, Grin2c, Col4a5, Heph, Celsr1, Egfr, Lgi4, Slc38a3, Aass, Fkbp10, Slc7a10, Cyp2d22, Cd38, Cyp4f14, Cbs, Slc1a2, Emp2, Axl, Slc14a1, Btbd17', ', ')[[1]], entrez.ID] %>%
  search %>% enrich %>% head(50) %>%
  .[, contrast := paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri) %>%
      factor(levels = unique(.), ordered = T)] %>%
  ggplot(aes(contrast, score - 9, fill = distance)) +
  geom_bar(stat = 'identity') + coord_flip() +
  scale_y_continuous(expand = c(0, 0), labels = function(x) x + 9) + theme_classic(20) +
  labs(x = element_blank(), y = 'Test Statistic') +
  scale_fill_fermenter(palette = 'Purples', name = 'Ontology Steps') +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1))

# Brain specific ----
list.files('/home/qinkaiwu/data/GTEx/tissueMean/') %>% .[grepl('^Brain', .)] %>%
  lapply(function(file) {
    fread(paste0('/home/qinkaiwu/data/GTEx/tissueMean/', file))
  }) -> mFiles

list.files('/home/qinkaiwu/data/GTEx/tissueMean/') %>% .[!grepl('^Brain', .)] %>%
  lapply(function(file) {
    fread(paste0('/home/qinkaiwu/data/GTEx/tissueMean/', file))
  }) -> mOtherFiles

lapply(mOtherFiles, function(x) {
  setorder(x, V2) %>% .[1:floor(nrow(.)/4), V1]
}) %>% unlist %>% unique

lapply(mFiles, function(x) {
  setorder(x, -V2) %>% .[1:500, V1]
}) %>% unlist %>% unique %>%
  intersect(lapply(mOtherFiles, function(x) {
    setorder(x, V2) %>% .[1:floor(nrow(.)/4), V1]
  }) %>% unlist %>% unique
  ) %>% {
    DATA.HOLDER$human@gene.meta[ensembl.ID %in% gsub('\\..*$', '', .), gene.Name]
  } %>% paste0(collapse = ', ') %>% search %>% enrich

# VAD studies plot ----
mGenes <- strsplit('BCL6, BIRC3, CEBPD, ERRFI1, FBXL16, FKBP5, GADD45B, IRS2, KLF9, PDK4, PER1, RGCC, RGS2, SEC14L2, SLC16A12, TFCP2L1, TSC22D3', ', ')[[1]]
synchronise({
  geneExpression(CACHE.BACKGROUND$human[grepl('ventricular assist', cf.ValLongUri), ee.ID][1:2], CACHE.BACKGROUND$human[grepl('ventricular assist', cf.ValLongUri), rsc.ID][1:2], 'human',
                 DATA.HOLDER$human@gene.meta[gene.Name %in% mGenes, entrez.ID])
}) %>% {
  generateResultsPlot(rownames(.$expr), 'reference subject role vs. ventricular assist device', ., getConfig(), 'human', rownames(.$expr),
                      'reference subject role vs. ventricular assist device', 'Heatmap', 'Gene Expression')
}

# Spike-ins ----
library(gridExtra)
dope_scores %>%
  .[distance == 0] %>%
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
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  scale_color_brewer(palette = 'Purples', name = '# of Spike ins') +
  scale_size_discrete(name = '# of Spike ins') +
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

# Artificial histogram ----
search_grid2 %>% copy %>% .[, n := pmin(n, 5)] %>%
  .[, sum_n := sum(n == N), .(N = n, genes, random)] %>%
  .[genes == 1, c('random', 'sum_n') := list(F, sum_n/2)] %>%
  .[, random := ifelse(random, 'Yes', 'No')] %>%
  .[, random := factor(random, levels = c('Yes', 'No'), ordered = T)] %>%
  ggplot(aes(factor(n), sum_n, fill = random)) +
  geom_bar(stat = 'identity', width = 0.3, position = position_dodge(width = 0.4)) +
  facet_wrap(~genes, labeller = labeller(.default = function(x) paste(x, 'genes') %>% { gsub('1 genes', '1 gene', .) })) +
  theme_bw(base_size = 20) + labs(x = 'Number of top 5 contrasts in top 50', y = 'Count') +
  theme(strip.background = element_blank(), strip.text = element_text(face = 'bold', hjust = 0)) +
  scale_fill_brewer(palette = 'Dark2', name = 'Random') +
  scale_y_continuous(expand = c(0, 0))

search_grid2 %>% ggplot(aes(cor.score, fill = factor(genes))) +
  geom_histogram() +
  scale_fill_manual(name = 'Number of genes', values = rev(colorRampPalette(brewer.pal(8, 'Purples'))(length(unique(search_grid2$genes)) * 2))[1:length(unique(search_grid2$genes))] %>% rev) +
  theme_classic(20) + facet_wrap(~random) +
  labs(x = 'Correlation', y = 'Count') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = 'none') -> mPlot

(search_grid2 %>% ggplot(aes(cor.stat, cor.score, fill = genes)) +
    geom_tile() +
    scale_fill_fermenter(name = '# of genes', palette = 'Purples', direction = 1) +
    theme_classic(20)) %>% cowplot::get_legend() -> mLegend

cowplot::ggdraw(cowplot::plot_grid(mPlot, mLegend, rel_widths = c(1, 0.2)))

# Saturation ----
mConstrasts <- DATA.HOLDER$human@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
  unique %>% .[, paste0(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
mViewed <- list(0, 0, 0, 0, 0)

for(j in 1:5) {
  mSeen <- c()
  mNonVisited <- 1:nrow(DATA.HOLDER$human@experiment.meta) %>% sample
  
  for(mWhich in mNonVisited) {
    mViewing <- DATA.HOLDER$human@experiment.meta %>% .[mWhich, paste0(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
    mViewed[[j]] <- c(mViewed[[j]], mViewing %in% mSeen)
    mSeen <- unique(c(mSeen, mViewing))
  }
}

data.table(A = mViewed[[1]], B = mViewed[[2]], C = mViewed[[3]], D = mViewed[[4]], E = mViewed[[5]]) %>%
  .[, lapply(.SD, function(x) cumsum(!x))] %>% .[, x := .I] %>% melt(measure.vars = 1:j) %>%
  ggplot(aes(x, value, color = variable)) + geom_line(lwd = 1.5, alpha = 0.7) + theme_bw(20) + theme(legend.position = 'none') +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = 'Index', y = '# Unique Condition Comparisons')
