source('/home/jsicherman/Thesis Work/requirements.R')
source('/home/jsicherman/Thesis Work/dependencies.R')

library(parallel)

if(!exists('contrasts'))
  contrasts <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/contrast_aff.rds')
uContrasts <- contrasts[, unique(contrast)]

DATA.HOLDER[c('human', 'mouse', 'rat')] <- NULL
NULLS[c('human', 'mouse', 'rat')] <- NULL

options(mc.cores = 30)
mclapply(1:250, function(n) {
  lapply(c(T, F), function(k) {
    lapply(c(1:10, 20, 50), function(i) {
      system(paste0('echo ', n, k, i))
      if(k)
        GOIs <- sample(1:nrow(DATA.HOLDER$artificial@gene.meta), i)
      else
        GOIs <- contrasts[contrast == sample(uContrasts, 1)] %>% setorder(-probability) %>%
          .[, head(entrez.ID, i)]
      
      contrasts[entrez.ID %in% GOIs] %>%
        .[, .(probability = sum(probability)), contrast] %>%
        merge(
          merge(DATA.HOLDER$artificial@experiment.meta[, .(rsc.ID, cf.ID)],
                CACHE.BACKGROUND$artificial[, .(rsc.ID, cf.Cat, cf.BaseLongUri, cf.ValLongUri)],
                by = 'rsc.ID') %>% .[, !'rsc.ID'] %>% unique,
          by.x = 'contrast', by.y = 'cf.ID', all = T) %>%
        unique %>%
        setorder(-probability) %>%
        .[, rank := rank(-probability, ties.method = 'min')] %>%
        .[, rank := as.integer(as.character(factor(rank, labels = 1:length(unique(rank)))))] -> tmp
      
      mSearch <- search(GOIs, getConfig(taxa = 'artificial'))
      if(is.null(mSearch))
        return(NULL)
      
      mSearch %>% enrich(getConfig(taxa = 'artificial', categories = getConfig('categories')$choices %>% unlist %>% unname)) %>%
        .[, I := .I] %>%
        merge(tmp[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, rank, probability)],
              by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), all = T) %>%
        .[, .(stat = max(stat), I = min(I),
              rank = min(rank),
              probability = max(probability),
              score = max(score)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
        .[is.na(score), score := 0] %>%
        .[is.na(stat), stat := 0] %>%
        setorder(I) %>% {
          data.table(genes = i,
                     cor.stat = cor(.$stat, .$probability, method = 'spearman'),
                     cor.score = cor(.$score, .$probability, method = 'spearman'),
                     n = .[rank <= 5, sum(I <= 50, na.rm = T)],
                     random = k)
        }
    }) %>% rbindlist
  }) %>% rbindlist
}) %>% rbindlist -> new.rankings2
saveRDS(new.rankings2, '/space/scratch/jsicherman/Thesis Work/data/artificial/search_grid2.rds')

search_grid %>% copy %>% .[, n := pmin(n, 5)] %>%
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
