library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(htmlwidgets)
library(sparkline)
library(DT)
library(heatmaply)
library(shinyHeatmaply)
library(shinypanels)
library(RColorBrewer)
library(stringr)
library(bit)

library(jsonlite)
library(async)

library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)

library(gemmaAPI)
library(ermineR)
library(mygene)
library(homologene)

source('/home/jsicherman/Thesis Work/dependencies.R')

library(parallel)
options(mc.cores = 30)

lapply(getConfig(key = 'taxa')$core, function(x) {
  lapply(1:20, function(i) {
    mclapply(0:9999, function(j) {
      if(j %% 200 == 0)
        message(paste0(Sys.time(), ' ... ', x, ' ', i, ' ... ', round(j / 100, 2), '%'))
      
      mGenes <- sample(1:nrow(DATA.HOLDER[[x]]@gene.meta), i)

      ret <- try({
        tmp <- search(DATA.HOLDER[[x]]@gene.meta[mGenes, entrez.ID], getConfig(taxa = x), verbose = F)
        
        if(is.null(tmp) || class(tmp) == 'try-error') {
          message(paste0('No rankings on genes ', paste0(mGenes, collapse = ', ')))
          NULL
        } else {
          enrich(tmp, getConfig(taxa = x), verbose = F, inprod = F) %>%
            .[, .(cf.Cat = cf.Cat, cf.BaseLongUri = cf.BaseLongUri, cf.ValLongUri = cf.ValLongUri,
                  stat = A / B, index = .I)]
        }
      })
      
      if(j == 9999)
        message('Saving...')
      
      if(is.null(ret) || class(ret) == 'try-error') {
        message(paste0('Error on genes ', paste0(mGenes, collapse = ', ')))
        NULL
      } else ret
    }) %>% rbindlist %>%
      .[, .(st.mean = mean(stat), st.sd = sd(stat), st.md = median(stat),
            in.mean = mean(index), in.sd = sd(index), in.md = median(index)), .(cf.BaseLongUri, cf.ValLongUri)] %>%
      saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/nulls/', x, '_', i, '.rds'))
  })
})
