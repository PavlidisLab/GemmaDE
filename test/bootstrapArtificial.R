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

source('/home/jsicherman/Thesis Work/dependencies.R')
rm(DRUGBANK)
CACHE.BACKGROUND$human <- NULL
CACHE.BACKGROUND$mouse <- NULL
CACHE.BACKGROUND$rat <- NULL
DATA.HOLDER$human <- NULL
DATA.HOLDER$mouse <- NULL
DATA.HOLDER$rat <- NULL
TAGS$human <- NULL
TAGS$mouse <- NULL
TAGS$rat <- NULL
NULLS$human <- NULL
NULLS$mouse <- NULL
NULLS$rat <- NULL
NULLS$any <- NULL

library(parallel)
options(mc.cores = 19)

eContrasts <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/experiment_contrasts.rds')

whichOnes <- sample(2:nrow(DATA.HOLDER$artificial@experiment.meta), 1000) - 1
lapply(1:length(whichOnes), function(i) {
  message(paste0('Progress: ', round(100 * i/length(whichOnes), 2), '%'))
  N <- whichOnes[i]
  # Pick an experiment. We'll pick genes that are DE in it and search those
  
  mclapply(2:20, function(N_GENES) {
    # Pick a number of "true" genes to search
    
    lapply(0:min(20 - N_GENES, N_GENES), function(MIX_IN) {
      message(paste0(N_GENES, '/', MIX_IN))
      # Pick a number of "false" genes to search
      
      genes.A <- DATA.HOLDER$artificial@gene.meta[entrez.ID %in% as.character(which(DATA.HOLDER$artificial@data$adj.pv[, N] < 0.05)), .(entrez.ID, mfx.Rank)] %>% setorder(mfx.Rank) %>%
        .[1:N_GENES, .(mfx.Rank, entrez.ID)]
      
      if(MIX_IN == 0) {
        genes.B <- data.table()
      } else {
        genes.B <- DATA.HOLDER$artificial@gene.meta[entrez.ID %in% as.character(which(DATA.HOLDER$artificial@data$adj.pv[, N + 1] < 0.05)), .(entrez.ID, mfx.Rank)] %>% setorder(mfx.Rank) %>%
          .[1:MIX_IN, .(mfx.Rank, entrez.ID)]
      }
      
      rbind(genes.A, genes.B) %>% .[, unique(entrez.ID)] %>% search(getConfig(taxa = 'artificial')) -> tmp
      
      whichExperiments <- switch((MIX_IN == 0) + 1, N:(N+1), N)
      search_results <- tmp %>% .[, I := .I] %>% .[rn %in% eContrasts[whichExperiments, rsc.ID], .(rn, I)] %>%
        .[, mfx.A := median(genes.A$mfx.Rank)] %>% {
          if(nrow(genes.B) == 0)
            .[, mfx.B := NA_real_]
          else
            .[, mfx.B := median(genes.B$mfx.Rank)]
        }
      
      tmp2 <- enrich(tmp, getConfig(taxa = 'artificial'))
      enrich_results <- TAGS$artificial[rsc.ID %in% eContrasts[whichExperiments, rsc.ID]] %>%
        .[distance == 0, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
        merge(tmp2[, c('I', 'frac') := list(.I, .I / max(.I))],
              by = c('cf.Cat', 'cf.BaseLongUri', 'cf.ValLongUri'), sort = F) %>%
        .[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, distance, stat, I, frac)]
      
      list(search = search_results[, mixin := MIX_IN] %>%
             .[, entrez.A := list(as.integer(genes.A$entrez.ID))] %>% {
               if(nrow(genes.B) == 0)
                 .[, entrez.B := list(NA_integer_)]
               else
                 .[, entrez.B := list(As.integer(genes.B$entrez.ID))]
             },
           enrich = enrich_results[, mixin := MIX_IN])
    }) %>% {
      list(search = rbindlist(lapply(., '[[', 'search'))[, genes := N_GENES],
           enrich = rbindlist(lapply(., '[[', 'enrich'))[, genes := N_GENES])
    }
  }) %>% {
    list(search = rbindlist(lapply(., '[[', 'search')),
         enrich = rbindlist(lapply(., '[[', 'enrich')))
  }
}) %>% {
  list(search = rbindlist(lapply(., '[[', 'search')),
       enrich = rbindlist(lapply(., '[[', 'enrich')))
} %>% saveRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/bootstrapped.rds')
