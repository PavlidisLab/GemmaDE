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
options(mc.cores = 24)

ITERS <- 1000

rm(DRUGBANK, NULLS)

# Copied from server for "any" impl
tidyGenes <- function(genes, taxa) {
  if(taxa == 'any') {
    taxOptions <- getConfig(key = 'taxa')$core
    taxIDs <- getConfig(key = 'taxa')$mapping
    
    orthologs <- lapply(taxOptions, function(i) {
      lapply(taxOptions, function(j) {
        homologene(genes, inTax = unname(taxIDs[i]), outTax = unname(taxIDs[j])) %>% {
          if(nrow(.) > 0) {
            if(i == j) {
              .[, !duplicated(colnames(.))] %>% .[, grepl('_ID', colnames(.))] %>% {
                data.frame(.) %>% `colnames<-`(paste0(unname(taxIDs[i]), '_ID')) %>%
                  rbind(
                    data.frame(
                      DATA.HOLDER[[i]]@gene.meta[gene.Name %in% genes | entrez.ID %in% genes, entrez.ID]
                    ) %>% `colnames<-`(paste0(unname(taxIDs[i]), '_ID'))
                  )
              }
            } else
              .[, grepl('_ID', colnames(.))]
          }
        }
      }) %>% rbindlist(fill = T) %>% {
        if(nrow(.) > 0) .
      }
    }) %>% rbindlist(fill = T) %>% .[, names(.) := lapply(.SD, as.character)] %>%
      .[, key := fcoalesce(.)] %>%
      group_by(key) %>% summarise_all(~na.omit(unique(.))[1]) %>%
      as.data.table
    
    return(orthologs)
  }
  
  # Clean numerics (interpreted as entrez IDs) and remove them from further processing.
  cleanGenes <- suppressWarnings(Filter(function(x) !is.na(as.integer(x)), genes))
  genes <- genes[!(genes %in% cleanGenes)]
  
  # If it matches (ENSG|ENSMUS|ENSRNO)\d{11}, it's an Ensembl ID (for human, mouse or rat).
  if(length(genes) > 0) {
    ensembl <- grep('(ENSG|ENSMUS|ENSRNO)\\d{11}', genes, value = T)
    
    if(length(ensembl) != 0) {
      # Extract genes with a matching Ensembl ID and clean them too.
      ensembls <- DATA.HOLDER[[taxa]]@gene.meta[ensembl.ID %in% ensembl, .(entrez.ID, ensembl.ID)]
      cleanGenes <- c(cleanGenes, ensembls[, entrez.ID])
      genes <- genes[!(genes %in% ensembls[, ensembl.ID])]
    }
  }
  
  if(length(genes > 0)) {
    go <- grep('(GO:)\\d{7}', genes, value = T)
    
    if(length(go) != 0) {
      gos <- DATA.HOLDER[[taxa]]@go[id %in% go, entrez.ID]
      cleanGenes <- c(cleanGenes, gos)
      genes <- genes[!(genes %in% go)]
    }
  }
  
  # Try to match to gene names and descriptions.
  if(length(genes) > 0) {
    descriptors <- DATA.HOLDER[[taxa]]@gene.meta[gene.Name %in% genes | gene.Desc %in% genes,
                                                 .(entrez.ID, gene.Name, gene.Desc)]
    if(nrow(descriptors) != 0) {
      cleanGenes <- c(cleanGenes, descriptors[, entrez.ID])
      genes <- genes[!(genes %in% descriptors[, c(gene.Name, gene.Desc)])]
    }
  }
  
  # If anything is left, try to match it to gene aliases.
  if(length(genes) > 0) {
    aliases <- DATA.HOLDER[[taxa]]@gene.meta[, parseListEntry(alias.Name), entrez.ID] %>%
      .[grepl(paste0(genes, collapse = '|'), V1)]
    if(nrow(aliases) > 0)
      cleanGenes <- c(cleanGenes, aliases[, unique(entrez.ID)])
  }
  
  cleanGenes
}

# While doing artificial
DATA.HOLDER$human <- NULL
DATA.HOLDER$mouse <- NULL
DATA.HOLDER$rat <- NULL

lapply('artificial', function(x) { # c(getConfig(key = 'taxa')$core, 'artificial')
  lapply(1:20, function(i) {
    #if(file.exists(paste0('/space/scratch/jsicherman/Thesis Work/data/nulls/', x, '_', i, '.rds'))) {
    #  message(paste0('File for ', x, '_', i, ' already exists... Skipping.'))
    #} else {
      mclapply(1:ITERS, function(j) {
        if(j %% (ITERS / 20) == 0)
          message(paste0(Sys.time(), ' ... ', x, ' ', i, ' ... ', round(100 * j / ITERS, 2), '%'))
        
        if(x == 'any')
          mGenes <- sample(1:nrow(DATA.HOLDER$human@gene.meta), i)
        else
          mGenes <- sample(1:nrow(DATA.HOLDER[[x]]@gene.meta), i)
        
        ret <- try({
          if(x == 'any')
            searchGenes <- tidyGenes(DATA.HOLDER$human@gene.meta[mGenes, entrez.ID], x)
          else
            searchGenes <- DATA.HOLDER[[x]]@gene.meta[mGenes, entrez.ID]
          
          tmp <- search(searchGenes, getConfig(taxa = x), verbose = F)
          
          if(is.null(tmp) || class(tmp) == 'try-error') {
            message(paste0('No rankings on genes ', paste0(mGenes, collapse = ', ')))
            NULL
          } else
            enrich(tmp, getConfig(taxa = x), verbose = F, inprod = F) %>% .[, index := .I]
        })
        
        if(is.null(ret) || class(ret) == 'try-error') {
          message(paste0('Error on genes ', paste0(mGenes, collapse = ', ')))
          NULL
        } else ret
      }) %>% rbindlist %>% {
        message('Saving...')
        
        .[, .(st.mean = mean(stat, na.rm = T), st.sd = sd(stat, na.rm = T), in.mean = mean(index, na.rm = T)),
          .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
          .[st.mean != 0 | st.sd != 0] %>%
          saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/nulls/', x, '_', i, '.rds'))
        
        # Save full
        if(F) {
          .[, list(stat = list(stat), index = list(index)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
            saveRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/full_nulls/', x, '_', i, '.rds'))
        }
        
        gc()
        NULL
      #}
    }
  })
})
