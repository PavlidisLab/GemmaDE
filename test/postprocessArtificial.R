VERSION <- '_superuniform-binary'
VERSION2 <- substring(VERSION, 2)

artificial <- readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/Limma/', VERSION2, '/',
                             VERSION, '.rds'))

tmp <- rbindlist(lapply(1:length(artificial), function(i) {
  data.table(experiment = i, entrez.ID = artificial[[i]]$entrez.ID,
             fc = artificial[[i]]$fc, adj.pv = artificial[[i]]$adj.pv)
  }))

fc <- dcast(tmp[, .(entrez.ID, experiment, fc)], entrez.ID ~ experiment, value.var = 'fc')
pv <- dcast(tmp[, .(entrez.ID, experiment, adj.pv)], entrez.ID ~ experiment, value.var = 'adj.pv')

contrasts <- rbindlist(lapply(artificial, '[[', 'contrast'))
saveRDS(contrasts, paste0('/space/scratch/jsicherman/Thesis Work/data/Limma/', VERSION2, '/experiment.contrasts',
                          VERSION, '.rds'))

rm(tmp, artificial)

fc <- as.data.frame(fc)
pv <- as.data.frame(pv)

rownames(fc) <- paste0('g', fc[, 1])
rownames(pv) <- paste0('g', pv[, 1])

fc <- fc[, -1] %>% as.matrix
pv <- pv[, -1] %>% as.matrix

artificial.gene.associations <- readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/Limma/', VERSION2, '/gene.associations',
                                               VERSION, '.rds'))
artificial.gene.meta <- readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/Limma/', VERSION2, '/artificial.gene.meta',
VERSION, '.rds'))

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
                              ee.NumSamples = readRDS(paste0('/space/scratch/jsicherman/Thesis Work/data/Limma/', VERSION2,
                                                             '/samples', VERSION, '.rds')),
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
        paste0('/space/scratch/jsicherman/Thesis Work/data/Limma/', VERSION2, '/artificial', VERSION, '.rds'))

DATA.HOLDER$artificial <- new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
                              experiment.meta = experiment.meta, gene.meta = artificial.gene.meta)
rm(fc, pv, artificial.gene.meta, experiment.meta, N, EXPERIMENTS)

DATA.HOLDER$artificial@gene.meta <- DATA.HOLDER$artificial@gene.meta[, c('n.DE', 'dist.Mean', 'dist.SD') :=
                                                                       list(rowSums2(DATA.HOLDER$artificial@data$adj.pv < 0.05, na.rm = T),
                                                                            rowMeans2(DATA.HOLDER$artificial@data$fc, na.rm = T),
                                                                            Rfast::rowVars(DATA.HOLDER$artificial@data$fc, std = T, na.rm = T))]

DATA.HOLDER$artificial@data$zscore <- (DATA.HOLDER$artificial@data$fc - DATA.HOLDER$artificial@gene.meta$dist.Mean) / DATA.HOLDER$artificial@gene.meta$dist.SD
DATA.HOLDER$artificial@data$pvz <- DATA.HOLDER$artificial@data$zscore %>% {
  tmp <- DATA.HOLDER$artificial@data$adj.pv
  tmp[is.na(tmp)] <- 1
  tmp[tmp < 1e-20] <- 1e-20
  abs(.) * -log(tmp, 100)
}
DATA.HOLDER$artificial@experiment.meta$version <- VERSION

CACHE.BACKGROUND$artificial <- precomputeTags('artificial')
CACHE.BACKGROUND$artificial$version <- VERSION
