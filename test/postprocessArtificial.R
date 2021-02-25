artificial <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial4/experiment_dat.rds')
contrasts <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial4/contrast_aff.rds')
contrastMap <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial4/contrast_map.rds')

tmp <- rbindlist(lapply(1:length(artificial), function(i) {
  if(class(artificial[[i]]) == 'try-error' || is.null(artificial[[i]]))
    data.table(experiment = i, entrez.ID = NA_integer_, contrast = NA_integer_,
               fc = NA_real_, adj.pv = NA_real_)
  else {
    eID <- artificial[[i]]$entrez.ID
    if(is.null(eID))
      eID <- 1:length(artificial[[i]]$fc)
    
    data.table(experiment = i, entrez.ID = eID, contrast = artificial[[i]]$contrast,
               fc = artificial[[i]]$fc, adj.pv = artificial[[i]]$adj.pv)
  }
}))

fc <- dcast(tmp[, .(entrez.ID, experiment, fc)], entrez.ID ~ experiment, value.var = 'fc')
pv <- dcast(tmp[, .(entrez.ID, experiment, adj.pv)], entrez.ID ~ experiment, value.var = 'adj.pv')
eContrasts <- tmp[, unique(contrast), experiment] %>% .[, V1]

rm(tmp, artificial)

fc <- as.data.frame(fc) %>% .[!is.na(.$entrez.ID), ]
pv <- as.data.frame(pv) %>% .[!is.na(.$entrez.ID), ]

rownames(fc) <- paste0('g', fc[, 1])
rownames(pv) <- paste0('g', pv[, 1])

fc <- fc[, -1] %>% as.matrix
pv <- pv[, -1] %>% as.matrix

N <- ncol(fc)

letterWrap <- function(n, depth = 1) {
  x <- do.call(paste0,
               do.call(expand.grid, args = list(lapply(1:depth, function(x) return(LETTERS)), stringsAsFactors = F)) %>%
                 .[, rev(names(.[])), drop = F])
  
  if(n <= length(x)) return(x[1:n])
  
  return(c(x, letterWrap(n - length(x), depth = depth + 1)))
}

EXPERIMENTS <- letterWrap(N)
colnames(fc) <- EXPERIMENTS
colnames(pv) <- EXPERIMENTS
names(eContrasts) <- EXPERIMENTS

#goodExperiments <- which(!is.na(eContrasts)) %>% unname
#fc <- fc[, goodExperiments]
#pv <- pv[, goodExperiments]
#eContrasts <- eContrasts[goodExperiments]
#EXPERIMENTS <- EXPERIMENTS[goodExperiments]

eMeta <- DATA.HOLDER$human@experiment.meta %>% copy# %>% .[goodExperiments]
eMeta[, c('cf.Cat', 'cf.Baseline', 'cf.Val', 'cf.BaseLongUri', 'cf.ValLongUri') := NULL]

eMeta[, cf.ID := eContrasts]
eMeta <- merge(eMeta, contrastMap[, cf.ID := .I], by = 'cf.ID', sort = F)
#eMeta[, cf.Baseline := '']
#eMeta[, cf.Val := ''] # TODO
#eMeta[, ee.qScore := ...]
#eMeta[, ee.sScore := ...]
eMeta[, n.detect := nrow(fc)]
#eMeta[, mean.fc := ...]
eMeta[, ee.Name := EXPERIMENTS]
eMeta[, rsc.ID := EXPERIMENTS]

#eMeta[, n.DE := colSums2(pv < 0.05, na.rm = T)]
eMeta[, mean.fc := colMeans2(fc, na.rm = T)]

gMeta <- data.table(entrez.ID = as.character(1:nrow(fc)),
                    gene.ID = NA_character_,
                    ensembl.ID = NA_character_,
                    gene.Name = paste0('g', 1:nrow(fc)),
                    alias.Name = NA_character_,
                    gene.Desc = NA_character_,
                    mfx.Rank = contrasts[, sum(probability), entrez.ID][, sqrt(V1 / max(V1))])

DATA.HOLDER$artificial <- new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
                              experiment.meta = eMeta, gene.meta = gMeta,
                              go = data.table(entrez.ID = NA, category = NA, id = NA, term = NA))
rm(fc, pv, N, EXPERIMENTS, eMeta, gMeta, eContrasts, contrastMap)

DATA.HOLDER$artificial@gene.meta <- DATA.HOLDER$artificial@gene.meta[, c('n.DE', 'dist.Mean', 'dist.SD') :=
                                                                       list(rowSums2(DATA.HOLDER$artificial@data$adj.pv < 0.05, na.rm = T),
                                                                            rowMeans2(DATA.HOLDER$artificial@data$fc, na.rm = T),
                                                                            Rfast::rowVars(DATA.HOLDER$artificial@data$fc, std = T, na.rm = T))]

DATA.HOLDER$artificial@data$zscore <- (DATA.HOLDER$artificial@data$fc - DATA.HOLDER$artificial@gene.meta$dist.Mean) / DATA.HOLDER$artificial@gene.meta$dist.SD

CACHE.BACKGROUND$artificial <- precomputeTags('artificial')
