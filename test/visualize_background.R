library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
library(parallel)

options('mc.cores' = 6)
theme_set(theme_cowplot())

SUMMARY$OBI[rsc.ID %in% unique(SUMMARY$OBI$rsc.ID)[10:20]] %>%
  ggplot(aes(score, color = rsc.ID)) + geom_density(fill = NA) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  xlab('Experiment Score') + ylab('Density') + theme(legend.position = 'none')
