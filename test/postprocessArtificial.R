library(dplyr)
library(data.table)
library(matrixStats)

artificial <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial.rds')
artificial <- lapply(artificial, function(x) x[, entrez.ID := as.integer(substring(entrez.ID, 2))])
artificial <- lapply(1:length(artificial), function(i) data.table(artificial[[i]], experiment = i))
tmp <- rbindlist(artificial)
fc <- dcast(tmp[, .(entrez.ID, experiment, fc)], entrez.ID ~ experiment, value.var = 'fc')
pv <- dcast(tmp[, .(entrez.ID, experiment, adj.pv)], entrez.ID ~ experiment, value.var = 'adj.pv')
rm(tmp, artificial)

artificial.exp.meta <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial.experiment.meta.rds')
artificial.gene.meta <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial.gene.meta.rds')

fc <- as.data.frame(fc)
pv <- as.data.frame(pv)

rownames(fc) <- paste0('g', fc[, 1])
rownames(pv) <- paste0('g', pv[, 1])

fc <- fc[, -1]
pv <- pv[, -1]

colnames(fc) <- artificial.exp.meta$rsc.ID
colnames(pv) <- artificial.exp.meta$rsc.ID

saveRDS(new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
            experiment.meta = artificial.exp.meta, gene.meta = artificial.gene.meta),
        paste0('/space/scratch/jsicherman/Thesis Work/data/', ifelse(USE_DESEQ, 'DESeq2', 'Limma'), '/artificial.rds'))

source('main/load.R')

DATA.HOLDER$artificial <- new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
                              experiment.meta = artificial.exp.meta, gene.meta = artificial.gene.meta)
DATA.HOLDER$artificial@gene.meta$n.DE <- rowSums2(DATA.HOLDER$artificial@data$adj.pv < 0.05, na.rm = T)
rm(fc, pv, artificial.gene.meta, artificial.exp.meta)

CACHE.BACKGROUND$artificial <- precomputeTags('artificial')
