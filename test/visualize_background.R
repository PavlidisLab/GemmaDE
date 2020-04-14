library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(parallel)

options('mc.cores' = 6)
theme_set(theme_cowplot())

SUMMARY <- readRDS('data/bootstrap/bootstrap.OBI.rds')
SUMMARY <- do.call(rbind, lapply(SUMMARY, function(x) data.table(rsc.ID = names(x$experiments), score = x$experiments))) %>%
  as.data.table

SUMMARY[rsc.ID %in% (DATA.HOLDER$human@experiment.meta[, .(n.DE, rsc.ID)] %>% setorder(-n.DE) %>% .[1:10, rsc.ID])] %>%
  as.data.frame %>% ggplot(aes(score, color = rsc.ID)) + geom_density(fill = NA) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  xlab('Experiment Score') + ylab('Density') + theme(legend.position = 'none')

boot_long[Definition %in% c('male', 'heart', 'brain', 'liver', 'adipose tissue')] %>%
  ggplot(aes(100 * N.x / (N.y + N.x), color = Definition)) + geom_density(fill = NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + scale_y_continuous(expand = c(0, 0)) +
  xlab('Percentage of Tags Detected (%)') + ylab('Density') + theme(legend.position = c(0.7, 0.7))
