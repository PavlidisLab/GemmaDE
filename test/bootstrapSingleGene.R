# Visualizations
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders) # From jsicherman/shinycssloaders, NOT daattali
library(htmlwidgets)
library(DT)
library(heatmaply)
library(shinyHeatmaply)
library(shinypanels) # From jsicherman/shinypanels, NOT datasketch
library(circlepackeR)
library(d3wordcloud)
library(data.tree)
# library(sparkline)
library(RColorBrewer)
library(sass)

library(async)
library(memoise)

# Data drivers
library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)
library(stringr)
library(bit)

# Parsing helpers
library(gemmaAPI, lib.loc = '/home/omancarci/R/x86_64-redhat-linux-gnu-library/3.6/')
library(ermineR)
library(mygene)
library(homologene)
library(jsonlite)
library(XML)

# Concurrent users
library(promises)
library(future)
plan(multicore, workers = 5)
options(future.globals.maxSize = 30 * 1000 * 1024^2)

library(parallel)

source('/home/jsicherman/Thesis Work/dependencies.R')
rm(NULLS)
DATA.HOLDER[c('artificial', 'mouse', 'rat')] <- NULL

options(mc.cores = 3)

mRange <- 1:length(DATA.HOLDER$human@gene.meta$entrez.ID) %>% split(ceiling(seq_along(.[]) / 500))
mclapply(names(mRange), function(mID) {
  mBlock <- lapply(1:length(mRange[[mID]]), function(indx) {
    i <- mRange[[mID]][indx]
    if(indx %% 25 == 0)
      message(paste0(mID, ': ', (100 * indx/length(mRange[[mID]])), '%'))
    
    tmp <- search(DATA.HOLDER$human@gene.meta$entrez.ID[i])
    if(is.null(tmp)) NULL
    else {
      tmp %>% enrich %>%
      .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, distance, A, B, C, D)]
    }
  })
  
  saveRDS(mBlock, paste0('/space/scratch/jsicherman/Thesis Work/data/singlegene/human_', mID, '.rds'))
  rm(mBlock)
  gc()
  NULL
})
