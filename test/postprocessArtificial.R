artificial <- readRDS('/space/scratch/jsicherman/Thesis Work/data/Limma/normalized/artificial_normalized.rds')

tmp <- rbindlist(lapply(1:length(artificial), function(i) {
  data.table(experiment = i, entrez.ID = artificial[[i]]$entrez.ID,
             fc = artificial[[i]]$fc, adj.pv = artificial[[i]]$adj.pv)
  }))

fc <- dcast(tmp[, .(entrez.ID, experiment, fc)], entrez.ID ~ experiment, value.var = 'fc')
pv <- dcast(tmp[, .(entrez.ID, experiment, adj.pv)], entrez.ID ~ experiment, value.var = 'adj.pv')

contrasts <- rbindlist(lapply(artificial, '[[', 'contrast'))
saveRDS(contrasts, '/space/scratch/jsicherman/Thesis Work/data/Limma/normalized/experiment.contrasts_normalized.rds')

rm(tmp, artificial)

fc <- as.data.frame(fc)
pv <- as.data.frame(pv)

rownames(fc) <- paste0('g', fc[, 1])
rownames(pv) <- paste0('g', pv[, 1])

fc <- fc[, -1]
pv <- pv[, -1]

artificial.gene.associations <- readRDS('/space/scratch/jsicherman/Thesis Work/data/Limma/normalized/gene.associations_normalized.rds')
artificial.gene.meta <- readRDS('/space/scratch/jsicherman/Thesis Work/data/Limma/normalized/artificial.gene.meta_normalized.rds')

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
                              cf.GeneDrive = contrasts[, entrez.ID],
                              cf.Cat = contrasts[, cf.Cat] %>% as.factor,
                              cf.CatLongUri = contrasts[, cf.CatLongUri] %>% as.factor,
                              cf.ValLongUri = contrasts[, cf.ValLongUri] %>% as.factor,
                              cf.BaseLongUri = contrasts[, cf.BaseLongUri] %>% as.factor,
                              n.DE = colSums2(pv %>% as.matrix < 0.05, na.rm = T),
                              mean.fc = colMeans2(fc %>% as.matrix, na.rm = T))

saveRDS(new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
            experiment.meta = experiment.meta, gene.meta = artificial.gene.meta),
        '/space/scratch/jsicherman/Thesis Work/data/Limma/normalized/artificial_normalized.rds')

DATA.HOLDER$artificial <- new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
                              experiment.meta = experiment.meta, gene.meta = artificial.gene.meta)
rm(fc, pv, artificial.gene.meta, experiment.meta, N, EXPERIMENTS)

DATA.HOLDER$artificial@gene.meta$n.DE <- rowSums2(DATA.HOLDER$artificial@data$adj.pv < 0.05, na.rm = T)
CACHE.BACKGROUND$artificial <- precomputeTags('artificial')
