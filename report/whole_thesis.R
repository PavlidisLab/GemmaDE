# this is an attempt to have an up to date working version of all plots and claims
# included in Jordan's thesis, in order of appearance -Ogan

devtools::load_all()
source(here::here("main/requirements.R"))
source(here::here("main/dependencies.R"))

# Chapter 1, figure 1.2 GEO and Gemma Contents----------
experiments = readRDS('report/experiments.rds')

# did not regenerate this file from scratch. generate/quantifyGEO.R deals with it
# and needs to be tested
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


# Chapter 2 Figure 2.2 --------


# Chapter 3 Figure 3.1 Overlaps -----------
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
