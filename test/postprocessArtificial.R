artificial <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/experiment_data.rds')
#contrasts <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/contrast_table.rds')
exp.meta <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/experiment_meta.rds')
gene.meta <- readRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/gene_meta.rds')

tmp <- rbindlist(lapply(1:length(artificial), function(i) {
  if(class(artificial[[i]]) == 'try-error')
    data.table(experiment = i, entrez.ID = NA_integer_,
               fc = NA_real_, adj.pv = NA_real_)
  else
    data.table(experiment = i, entrez.ID = artificial[[i]]$entrez.ID,
               fc = artificial[[i]]$fc, adj.pv = artificial[[i]]$adj.pv)
  }))

fc <- dcast(tmp[, .(entrez.ID, experiment, fc)], entrez.ID ~ experiment, value.var = 'fc')
pv <- dcast(tmp[, .(entrez.ID, experiment, adj.pv)], entrez.ID ~ experiment, value.var = 'adj.pv')

eContrasts <- rbindlist(lapply(artificial[sapply(artificial, class) != 'try-error'], '[[', 'contrast'))

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

exp.meta[, n.DE := colSums2(as.matrix(pv) < 0.05, na.rm = T)]
exp.meta[, mean.fc := colMeans2(as.matrix(fc), na.rm = T)]

DATA.HOLDER$artificial <- new('EData', taxon = 'artificial', data = list(fc = fc, adj.pv = pv),
                              experiment.meta = exp.meta, gene.meta = gene.meta,
                              go = data.table(entrez.ID = NA, category = NA, id = NA, term = NA))
rm(fc, pv, N, EXPERIMENTS, exp.meta, gene.meta)

DATA.HOLDER$artificial@gene.meta <- DATA.HOLDER$artificial@gene.meta[, c('n.DE', 'dist.Mean', 'dist.SD') :=
                                                                       list(rowSums2(DATA.HOLDER$artificial@data$adj.pv < 0.05, na.rm = T),
                                                                            rowMeans2(DATA.HOLDER$artificial@data$fc, na.rm = T),
                                                                            Rfast::rowVars(DATA.HOLDER$artificial@data$fc, std = T, na.rm = T))]

DATA.HOLDER$artificial@data$zscore <- (DATA.HOLDER$artificial@data$fc - DATA.HOLDER$artificial@gene.meta$dist.Mean) / DATA.HOLDER$artificial@gene.meta$dist.SD

CACHE.BACKGROUND$artificial <- precomputeTags('artificial')
TAGS$artificial <- getTags('artificial')
