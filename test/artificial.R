library(data.table)
library(dplyr)
library(matrixStats)

set.seed(18232)

letterWrap <- function(n, depth = 1) {
  x <- do.call(paste0,
               do.call(expand.grid, args = list(lapply(1:depth, function(x) return(LETTERS)), stringsAsFactors = F)) %>%
                 .[, rev(names(.[])), drop = F])

  if(n <= length(x)) return(x[1:n])
  
  return(c(x, letterWrap(n - length(x), depth = depth + 1)))
}

rpdist <- function(n, nna) {
  sample(c(rep(NaN, nna),
           sapply(c(rgamma(floor((n - nna) / 2), runif(1, 0.3, 0.8)),
                    runif(ceiling((n - nna) / 2))), function(x) min(x, 1))))
}

rontdist <- function(n) {
  rdist <- function(n) {
    sample(c('', ONTOLOGIES$ChildNode_Long %>% as.character), n, T,
           c(0.5, rep(1 / length(ONTOLOGIES$ChildNode_Long),
                      length(ONTOLOGIES$ChildNode_Long))))
  }
  
  gsub('(^; |; $)', '', paste(rdist(n), rdist(n), sep = '; ')) %>% as.factor
}

experiments <- letterWrap(10000)
genes <- as.character(1:20000)

# Generate random FCs
fc <- matrix(sample(c(NaN, rnorm(length(experiments) * length(genes) / 2, mean = 0.03, sd = 4)),
                    length(experiments) * length(genes), T,
                    c(0.1, rep(1 / (length(experiments) * length(genes) / 2), length(experiments) * length(genes) / 2))),
             nrow = length(genes))
rownames(fc) <- genes
colnames(fc) <- experiments

# Generate random P-values and match their position with NAs in fc
adj.pv <- do.call(cbind, lapply(1:length(experiments), function(x) {
  ps <- rpdist(length(genes), 0)
  ps[which(is.na(fc[, x]))] <- NaN
  ps
}))
rownames(adj.pv) <- genes
colnames(adj.pv) <- experiments

metaData <- data.table(rsc.ID = experiments,
                       ee.ID = as.character(1:length(experiments)),
                       ee.Name = experiments,
                       ee.Source = 'TODO' %>% as.factor,
                       ee.NumSamples = runif(length(experiments), 1, 100),
                       ee.TagLongUri = rontdist(length(experiments)),
                       ad.Name = 'TODO' %>% as.factor,
                       ad.Company = 'TODO' %>% as.factor,
                       ad.Sequencing = F %>% as.factor,
                       sf.Subset = F,
                       sf.Cat = 'TODO' %>% as.factor,
                       sf.CatLongUri = 'TODO' %>% as.factor,
                       sf.ValLongUri = 'TODO' %>% as.factor,
                       cf.Cat = 'TODO' %>% as.factor,
                       cf.CatLongUri = 'TODO' %>% as.factor,
                       cf.ValLongUri = rontdist(length(experiments)),
                       cf.BaseLongUri = rontdist(length(experiments)))

metaData$n.DE <- colSums2(adj.pv < 0.05, na.rm = T)
metaData$mean.fc <- colMeans2(fc, na.rm = T)

metaGene <- data.table(entrez.ID = genes,
                       gene.ID = as.numeric(genes),
                       ensembl.ID = paste0('ENSJ', genes),
                       gene.Name = paste0('GENE', genes),
                       alias.Name = '',
                       gene.Desc = '',
                       mfx.Rank = runif(length(genes), 0.01, 0.99))

DATA.HOLDER$artificial <- new('EData', taxon = 'artificial',
                              data = list(fc = fc, adj.pv = adj.pv),
                              experiment.meta = metaData, gene.meta = metaGene)

rm(metaGene, metaData, fc, adj.pv, experiments, genes)
