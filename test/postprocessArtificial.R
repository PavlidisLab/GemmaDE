artificial <- readRDS('/space/scratch/jsicherman/Thesis Work/data/Limma/artificial.rds')

tmp <- rbindlist(lapply(1:length(artificial), function(i) {
  data.table(experiment = i, entrez.ID = artificial[[i]]$entrez.ID,
             fc = artificial[[i]]$fc, adj.pv = artificial[[i]]$adj.pv)
  }))

fc <- dcast(tmp[, .(entrez.ID, experiment, fc)], entrez.ID ~ experiment, value.var = 'fc')
pv <- dcast(tmp[, .(entrez.ID, experiment, adj.pv)], entrez.ID ~ experiment, value.var = 'adj.pv')

contrasts <- rbindlist(lapply(artificial, '[[', 'contrast'))

rm(tmp, artificial)

fc <- as.data.frame(fc)
pv <- as.data.frame(pv)

rownames(fc) <- paste0('g', fc[, 1])
rownames(pv) <- paste0('g', pv[, 1])

fc <- fc[, -1]
pv <- pv[, -1]

artificial.gene.associations <- readRDS('/space/scratch/jsicherman/Thesis Work/data/Limma/gene.associations.rds')
artificial.gene.meta <- readRDS('/space/scratch/jsicherman/Thesis Work/data/Limma/artificial.gene.meta.rds')

N <- ncol(fc)
EXPERIMENTS <- letterWrap(N)
colnames(fc) <- EXPERIMENTS
colnames(pv) <- EXPERIMENTS

experiment.meta <- data.table(rsc.ID = EXPERIMENTS,
                              ee.ID = 1:N,
                              ee.Name = EXPERIMENTS,
                              ee.Source = 'Artificial' %>% as.factor,
                              ee.NumSamples = readRDS('/space/scratch/jsicherman/Thesis Work/data/Limma/samples.rds'),
                              # TODO ee.TagLongUri = exp.assoc[, ee.TagLongUri] %>% as.factor,
                              ee.qScore = rnorm(N, 0.2, 0.3) %>% pmin(1), # TODO This should be more meaningful
                              ad.Name = 'TODO' %>% as.factor,
                              ad.Company = 'TODO' %>% as.factor,
                              ad.Sequencing = T %>% as.factor,
                              sf.Subset = F,
                              sf.Cat = 'TODO' %>% as.factor,
                              sf.CatLongUri = 'TODO' %>% as.factor,
                              sf.ValLongUri = 'TODO' %>% as.factor,
                              cf.Cat = contrasts[, cf.Cat] %>% as.factor,
                              cf.CatLongUri = contrast[, cf.CatLongUri] %>% as.factor,
                              cf.ValLongUri = contrasts[, cf.ValLongUri] %>% as.factor,
                              cf.BaseLongUri = contrasts[, cf.BaseLongUri] %>% as.factor,
                              n.DE = colSums2(pv < 0.05, na.rm = T),
                              mean.fc = colMeans2(fc, na.rm = T))

saveRDS(new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
            experiment.meta = experiment.meta, gene.meta = artificial.gene.meta),
        paste0('/space/scratch/jsicherman/Thesis Work/data/', ifelse(USE_DESEQ, 'DESeq2', 'Limma'), '/artificial.rds'))

DATA.HOLDER$artificial <- new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
                              experiment.meta = experiment.meta, gene.meta = artificial.gene.meta)
DATA.HOLDER$artificial@gene.meta$n.DE <- rowSums2(DATA.HOLDER$artificial@data$adj.pv < 0.05, na.rm = T)
rm(fc, pv, artificial.gene.meta, experiment.meta)
