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
options(mc.cores = 14)

ITERS <- 10000

OPTIONS <- c('human', 'artificial', 'mouse', 'rat')

for(x in OPTIONS) {
  if(!exists('DATA.HOLDER') || names(DATA.HOLDER) != x) {
    rm(DATA.HOLDER)
    source('/home/jsicherman/Thesis Work/dependencies.R')
    rm(NULLS)
    DATA.HOLDER[OPTIONS[-which(OPTIONS == x)]] <- NULL
  }
  
  for(i in 1:20) {
    message(paste0(Sys.time(), ' ... Starting ', x, ' ', i))
    if(file.exists(paste0('/space/scratch/jsicherman/Thesis Work/data/nulls2/', x, '_', i, '.rds'))) {
      message(paste0('File for ', x, '_', i, ' already exists... Skipping.'))
    } else {
      mclapply(1:ITERS, function(j) {
        if(j %% (ITERS / 20) == 0)
          message(paste0(Sys.time(), ' ... ', x, ' ', i, ' ... ', round(100 * j / ITERS, 2), '%'))
        
        tmp <- DATA.HOLDER[[x]]@gene.meta[sample(1:nrow(DATA.HOLDER[[x]]@gene.meta), i), entrez.ID] %>%
          search(getConfig(taxa = x))
        
        if(is.null(tmp)) {
          message('No rankings on genes')
          NULL
        } else {
          tmp[, .(score, rn)]
          #enrich(tmp, getConfig(taxa = x), inprod = F) %>% .[, index := .I] %>%
          #  .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, A, B, index, stat)]
        }
      }) %>% rbindlist %>% {
        message('Saving...')
        
        #.[, .(A.mean = mean(A, na.rm = T),
        #      B.mean = mean(B, na.rm = T),
        #      st.mean = mean(stat, na.rm = T),
        #      st.sd = sd(stat, na.rm = T),
        #      in.mean = mean(index, na.rm = T)),
        #  .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
        #  .[st.mean != 0 | st.sd != 0] %>%
        
        .[, .(score.mean = mean(score, na.rm = T),
              score.median = median(score, na.rm = T),
              score.sd = sd(score, na.rm = T)), rn] %>%
          saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/nulls2/', x, '_', i, '.rds'))
        
        # Save full
        if(F) {#i %in% c(1, 5, 10, 20)) {
          .[, list(stat = list(stat),
                   index = list(index)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
            saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/full_nulls/', x, '_', i, '.rds'))
        }
        
        gc()
        NULL
      }
    }
  }
}
