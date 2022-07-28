artificial <- readRDS(paste(DATADIR, 'artificial/experiment_dat.rds', sep='/'))
contrasts <- readRDS(paste(DATADIR, 'artificial/contrast_aff.rds', sep='/'))

tmp <- data.table::rbindlist(lapply(1:length(artificial), function(i) {
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

fc <- data.table::dcast(tmp[, .(entrez.ID, experiment, fc)], entrez.ID ~ experiment, value.var = 'fc')
pv <- data.table::dcast(tmp[, .(entrez.ID, experiment, adj.pv)], entrez.ID ~ experiment, value.var = 'adj.pv')
eContrasts <- tmp[, unique(contrast), experiment] %>% .[, V1]

rm(tmp, artificial)

fc <- as.data.frame(fc)
pv <- as.data.frame(pv)

rownames(fc) <- paste0('g', fc[, 1])
rownames(pv) <- paste0('g', pv[, 1])

fc <- fc[, -1] %>% as.matrix
pv <- pv[, -1] %>% as.matrix

N <- ncol(fc)

EXPERIMENTS <- letterWrap(N)
colnames(fc) <- EXPERIMENTS
colnames(pv) <- EXPERIMENTS
names(eContrasts) <- EXPERIMENTS

source(paste(PROJDIR, 'main/dependencies.R', sep='/'))

eMeta <- DATA.HOLDER$human@experiment.meta %>% data.table::copy()
eMeta[, c('cf.Cat', 'cf.Baseline', 'cf.Val', 'cf.BaseLongUri', 'cf.ValLongUri') := NULL]

eMeta[, cf.ID := eContrasts]
eMeta <- merge(eMeta, unique(DATA.HOLDER$human@experiment.meta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]) %>%
                 .[, cf.ID := .I], by = 'cf.ID', sort = F)
eMeta[, ee.Name := EXPERIMENTS]
eMeta[, rsc.ID := EXPERIMENTS]

eMeta[, n.DE := matrixStats::colSums2(pv < 0.05, na.rm = T)]
eMeta[, mean.fc := matrixStats::colMeans2(fc, na.rm = T)]

gMeta <- data.table::data.table(entrez.ID = as.character(1:nrow(fc)),
                    gene.ID = NA_character_,
                    ensembl.ID = NA_character_,
                    gene.Name = paste0('g', 1:nrow(fc)),
                    alias.Name = NA_character_,
                    gene.Desc = NA_character_,
                    mfx.Rank = contrasts[, sum(probability), entrez.ID][, sqrt(V1 / max(V1))])

DATA.HOLDER$artificial <- new('EData', taxon = 'artificial', data = list(adj.pv = pv),
                              experiment.meta = eMeta, gene.meta = gMeta,
                              go = data.table(entrez.ID = NA, category = NA, id = NA, term = NA))
rm(pv, N, EXPERIMENTS, eMeta, gMeta, eContrasts)

DATA.HOLDER$artificial@gene.meta <- DATA.HOLDER$artificial@gene.meta[, c('n.DE', 'dist.Mean', 'dist.SD') :=
                                                                       list(matrixStats::rowSums2(DATA.HOLDER$artificial@data$adj.pv < 0.05, na.rm = T),
                                                                            matrixStats::rowMeans2(fc, na.rm = T),
                                                                            Rfast::rowVars(fc, std = T, na.rm = T))]

DATA.HOLDER$artificial@data$zscore <- (fc - DATA.HOLDER$artificial@gene.meta$dist.Mean) / DATA.HOLDER$artificial@gene.meta$dist.SD

CACHE.BACKGROUND$artificial <- precomputeTags('artificial')

saveRDS(CACHE.BACKGROUND, paste(DATADIR, 'CACHE.BACKGROUND.rds', sep='/'))
saveRDS(DATA.HOLDER, paste(DATADIR, 'DATA.HOLDER.rds', sep='/'))
