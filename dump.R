library(gemmaAPI, lib.loc = '~/R/x86_64-redhat-linux-gnu-library/3.6/')
library(data.table)
library(gemmaAPI)
library(async)
library(glue)
library(dplyr)
library(jsonlite)

MAX_EXP <- 20000
BATCH_SIZE <- 1000

synchronise({
  async_map(split(as.character(1:MAX_EXP), ceiling(seq_along(as.character(1:MAX_EXP)) / BATCH_SIZE)),
            async(function(x) {
              gemmaAPI::getDatasets(x, limit = 0, async = T)
            }))
}) %>% rbindlist(fill = T) -> experiments

synchronise({
  async_map(experiments[, ee.ID], async(function(x) {
    gemmaAPI::getDatasetDEA(x, async = T)$then(function(anns) {
      if(length(anns) != 0 && anns[, any(is.na(ee.ID))])
        anns[, ee.ID := as.integer(ee.ID)] %>% .[, ee.ID := x]
      anns
    })
  }), .limit = floor(BATCH_SIZE / 5))
}) %>% rbindlist -> annotations

gemmaAPI::getPlatforms(limit = 0) %>%
  .[, .(ad.ID = platform.ID, ad.ShortName = platform.ShortName,
        ad.Troubled = platform.Troubled, ad.Sequencing = technology.Type == 'SEQUENCING',
        ad.Type = technology.Type, ad.Company = sapply(strsplit(platform.Name, ' '), '[[', 1))] -> platforms

synchronise({
  async_map(experiments[, ee.ID], async(function(x) {
    gemmaAPI::getDatasetAnnotations(x, async = T)$then(function(tags) {
      data.table(ee.ID = x, tags)
    })
  }), .limit = floor(BATCH_SIZE / 5))
}) %>% rbindlist(fill = T) -> tags

lapply(unique(annotations[, ad.ID]), function(ads) {
  platforms[ad.ID %in% ads,
            .(ad.ID = paste0(ads, collapse = ', '),
              ad.Name = paste0(ad.ShortName, collapse = ', '),
              ad.Company = paste0(ad.Company, collapse = ', '),
              ad.Sequencing = all(ad.Sequencing),
              ad.Troubled = any(ad.Troubled))]
}) %>% rbindlist %>%
  merge(annotations %>% .[, .(ad.ID = paste0(unlist(ad.ID), collapse = ', ')), setdiff(names(annotations), 'ad.ID')], by = 'ad.ID') %>%
  merge(experiments, by = 'ee.ID') %>%
  merge(tags[class.Type == 'ExperimentTag',
             .(ee.Cat = paste0(class.Name, collapse = ', '),
               ee.CatLongUri = paste0(class.URL, collapse = ', '),
               ee.Tag = paste0(term.Name, collapse = ', '),
               ee.TagLongUri = paste0(term.URL, collapse = ', ')), ee.ID],
        by = 'ee.ID') %>%
  .[, .(ee.ID, rsc.ID, ee.Public, ee.Troubled, ee.Taxon = taxon.Name, ee.Database,
        ee.Description, ee.ShortName, ee.Accession, ee.Name,
        ee.BatchCorrected = geeq.batchCorrected, ee.Batch.P, ee.Batch.PC,
        ee.qScore = geeq.qScore, ee.sScore = geeq.sScore,
        ee.NumSamples = ee.Samples, ee.Cat, ee.CatLongUri, ee.Tag, ee.TagLongUri,
        ee.numDE = stats.DE, ee.numDown = stats.Down, ee.numUp = stats.Up,
        ee.numProbes = probes.Analyzed, ee.numGenes = genes.Analyzed,
        ad.Name, ad.Sequencing, ad.Troubled, ad.Type = technology.Type,
        sf.Subset, sf.Cat, sf.CatLongUri, sf.Val, sf.ValLongUri,
        cf.Cat, cf.CatLongUri, cf.Val, cf.ValLongUri, cf.Baseline, cf.BaseLongUri
        )] %>% setorder(ee.Taxon, rsc.ID) -> experiment.meta

saveRDS(experiment.meta, '/space/scratch/jsicherman/Gemma/meta_exp.rds')

rm(experiments, platforms, tags, annotations)

# This bootlegged auto-retry is pretty gross, but should work

# Initialize with RSC IDs to fetch
# TODO Update later to only redownload updated ones.
for(taxon in na.omit(experiment.meta[, unique(ee.Taxon)])) {
  assign(paste0('QUEUED_', taxon), experiment.meta[ee.ID %in% experiment.meta[ee.Taxon == taxon, unique(ee.ID)]] %>%
           .[, rsc.ID])
}

rm(list = ls(pattern = 'diffEx_'))

# Keep trying to pull data down as long as we have some left
while(length(ls(pattern = 'QUEUED_')) > 0) {
  for(taxon in na.omit(experiment.meta[, unique(ee.Taxon)])) {
    # Finish a taxon
    if(!exists(paste0('QUEUED_', taxon)) || length(get(paste0('QUEUED_', taxon))) == 0) {
      rm(list = paste0('QUEUED_', taxon))
      next
    }

    message(taxon)

    anns <- experiment.meta[ee.ID %in% experiment.meta[ee.Taxon == taxon, unique(ee.ID)]] %>%
      .[rsc.ID %in% get(paste0('QUEUED_', taxon))] %>%
      .[, result.ID := as.integer(sapply(strsplit(rsc.ID, '.', T), '[[', 2))] %>%
      .[, .(result.ID, ee.ID, rsc.ID, cf.Cat, cf.Baseline, cf.Val)]

    synchronise({
      async_map(1:nrow(anns), async(function(i) {
        message(paste0(i, '/', nrow(anns), '...'))
        passthrough <- anns[i, .(rsc.ID, cf.Cat, cf.Baseline, cf.Val)]

        # Pull the DE results
        gemmaAPI::getDatasetDE(anns[i, ee.ID], keepNonSpecific = F, consolidate = 'average', limit = 0, diffExSet = anns[i, result.ID], async = T)$then(function(de) {
          # Test if our HTTP request worked...
          # Sometimes we get 500 errors if the server times out on us,
          # so we'll have to come back to those
          if('data.table' %in% class(de)) {
            # Unqueue safely if there are no associated vectors
            if(nrow(de[, unlist(expr, F)]) == 0)
              list(unqueue = passthrough[, rsc.ID])
            else { # Compute fold changes otherwise
              sampleInfo <- gemmaAPI::getDatasetSamples(de[, ee.ID]) %>%
                .[, .(sample.Name, sample.Characteristics, sample.FactorValues)]

              # Decide if each sample is in the baseline or the contrast
              baselines <- sampleInfo[sampleInfo[, sapply(sample.FactorValues, function(x) {
                sum(sapply(x[, name], grepl, passthrough[, cf.Baseline])) > sum(sapply(x[, name], grepl, passthrough[, cf.Val]))
              })], sample.Name]

              list(unqueue = passthrough[, rsc.ID],
                   data = de[, unlist(expr, F)] %>% .[, setdiff(sampleInfo[, sample.Name], baselines), with = F] %>%
                     rowMeans(na.rm = T) %>% `/`(de[, unlist(expr, F)] %>% .[, baselines, with = F] %>% rowMeans(na.rm = T)) %>% {
                       data.table(entrez.ID = de[, unlist(expr, F)][, gene.ID], fc = log2(.))
                     } %>% setnames(2, passthrough[, rsc.ID])
              )
            }
          }
        })
      }), .limit = floor(BATCH_SIZE / 30))
    }) %>% {
      # Unqueue non-failures
      assign(paste0('QUEUED_', taxon),
             Filter(function(x) !(x %in% unlist(lapply(., '[[', 'unqueue'))), get(paste0('QUEUED_', taxon))), envir = globalenv())

      # Recompute the differential expression set if there are any successes
      if(any(sapply(., length) == 2)) {
        diffEx <- Reduce(function(...) merge(..., all = T, by = 'entrez.ID'), lapply(.[lengths(.) == 2], '[[', 'data'))

        if(exists(paste0('diffEx_', taxon)))
          assign(paste0('diffEx_', taxon), merge(get(paste0('diffEx_', taxon)), diffEx, by = 'entrez.ID', all = T), envir = globalenv())
        else
          assign(paste0('diffEx_', taxon), diffEx, envir = globalenv())
        saveRDS(get(paste0('diffEx_', taxon)), paste0('/space/scratch/jsicherman/Gemma/diffEx_', taxon, '.rds'))
      }

      invisible(NULL)
    }
  }
}
