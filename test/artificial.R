library(data.table)
library(dplyr)
library(compcodeR)
library(DESeq2)
library(pbapply)
library(edgeR)
library(matrixStats)
library(dqrng)
library(matrixStats)
library(parallel)

mu.phi.estimates <- system.file("extdata", "Pickrell.Cheung.Mu.Phi.Estimates.rds",
                                package = "compcodeR")
mu.phi.estimates <- readRDS(mu.phi.estimates)

# Copied from compcodeR::generateSyntheticData to optimize
generateSyntheticData <- function(dataset, n.vars, samples.per.cond, n.diffexp, repl.id = 1, 
                                  seqdepth = 1e7, minfact = 0.7, maxfact = 1.4, 
                                  relmeans = "auto", dispersions = "auto", 
                                  fraction.upregulated = 1, between.group.diffdisp = FALSE, 
                                  filter.threshold.total = 1, filter.threshold.mediancpm = 0, 
                                  fraction.non.overdispersed = 0, random.outlier.high.prob = 0, 
                                  random.outlier.low.prob = 0, single.outlier.high.prob = 0, 
                                  single.outlier.low.prob = 0, effect.size = 1.5, 
                                  output.file = NULL) {
  ## Define conditions
  condition <- rep(c(1, 2), each = samples.per.cond)
  S1 <- which(condition == 1)
  S2 <- which(condition == 2)
  
  ## Define sets of upregulated, downregulated and non-differentially regulated genes
  genes.upreg <- which(effect.size > 1)
  genes.downreg <- which(effect.size < 1)
  genes.nonreg <- which(effect.size == 1)
  n.upregulated <- length(genes.upreg)
  n.diffexp <- length(genes.upreg) + length(genes.downreg)
  fraction.upregulated <- n.upregulated/n.diffexp
  
  ### Differentially expressed genes
  differential.expression <- rep(0, n.vars)
  differential.expression[genes.upreg] <- 1
  differential.expression[genes.downreg] <- 1
  upregulation <- rep(0, n.vars)
  upregulation[genes.upreg] <- 1
  downregulation <- rep(0, n.vars)
  downregulation[genes.downreg] <- 1
  
  ### Load mu and phi estimates from real data (Pickrell data set and Cheung data set)
  mu.estimates <- mu.phi.estimates$pickrell.cheung.mu
  phi.estimates <- mu.phi.estimates$pickrell.cheung.phi
  
  ### Sample a mu and a phi for each gene in condition S1
  to.include <- dqsample(1:length(mu.estimates), n.vars,
                         replace = ifelse(n.vars > length(mu.estimates), TRUE, FALSE))
  truedispersions.S1 <<- phi.estimates[to.include]
  truemeans.S1 <- mu.estimates[to.include]
    
  ### Generate sequencing depths (nfacts * Nk)
  seq.depths <<- dqrunif(2 * samples.per.cond, min = minfact, max = maxfact) * seqdepth
  
  ### Find rates of mapping to each gene in each condition
  prob.S1 <- truemeans.S1
  prob.S2 <- effect.size * prob.S1
  
  p.S1 <<- prob.S1 / sum(prob.S1)
  p.S2 <- prob.S2 / sum(prob.S2)
  
  ### Find new dispersions for condition S2, depending on what prob.S2 is. 
  ### From the mu/phi estimates, sample a phi value from 
  ### the pairs where mu is similar to prob.S2.
  truedispersions.S2 <- truedispersions.S1
  
  ### Initialize data matrix
  
  ### Generate data
  Z <- do.call(cbind, lapply(1:(length(S1) + length(S2)), function(j) {
    if(j %in% S1)
      rnbinom(n.vars, 1/truedispersions.S1, mu = p.S1 * seq.depths[j])
    else
      rnbinom(n.vars, 1/truedispersions.S2, mu = p.S2 * seq.depths[j])
  }))
  
  ### Assign variable names to rows
  rownames(Z) <- paste0('g', 1:nrow(Z))
  colnames(Z) <- paste0('sample', 1:ncol(Z))
  
  ### Create a sample annotation data frame
  sample.annotations <- data.frame(condition = condition)
  
  ### Filter the data with respect to total count
  Z <- Z[rowSums2(Z) >= filter.threshold.total, ]
  
  ### Filter the data with respect to median cpm
  Z <- Z[rowMedians(t(t(Z) / colSums2(Z))) >= (filter.threshold.mediancpm / 1e6), ]
  
  ### Generate sample and variable names
  rownames(sample.annotations) <- colnames(Z)
  
  list(count.matrix = Z, sample.annotations = sample.annotations)
}

letterWrap <- function(n, depth = 1) {
  x <- do.call(paste0,
               do.call(expand.grid, args = list(lapply(1:depth, function(x) return(LETTERS)), stringsAsFactors = F)) %>%
                 .[, rev(names(.[])), drop = F])
  
  if(n <= length(x)) return(x[1:n])
  
  return(c(x, letterWrap(n - length(x), depth = depth + 1)))
}

set.seed(18232)

N_EXPERIMENTS <- 11880
N_GENES <- 27151

experiments <- letterWrap(N_EXPERIMENTS)
genes.diff.prob <- rexp(N_GENES)
genes.diff.prob <- genes.diff.prob / (1.2 * max(genes.diff.prob))

genes.diff.prob.up <- runif(N_GENES)

gene.meta <- data.table(entrez.ID = paste0('g', 1:N_GENES),
                        gene.ID = 1:N_GENES,
                        ensembl.ID = paste0('ENSJ', 1:N_GENES),
                        gene.Name = paste0('GENE', 1:N_GENES),
                        alias.Name = '',
                        gene.Desc = '',
                        mfx.Rank = runif(N_GENES, 0.01, 0.99),
                        prob.DE = genes.diff.prob,
                        prob.DE.up = genes.diff.prob.up)

saveRDS(gene.meta, '/space/scratch/jsicherman/Thesis Work/data/artificial.gene.meta.rds')

experiment.samples <- 2 + as.integer(rgamma(N_EXPERIMENTS, 0.5, 8) * 1e3)

options(mc.cores = 30)
experiment.data <- mclapply(1:N_EXPERIMENTS, function(experiment) {
  cat(paste0(Sys.time(), ' ..... ', round(100 * experiment/N_EXPERIMENTS), '%\n'))
  genes.dysregulated <- sapply(1:N_GENES, function(x) {
    sample(-1:1, 1, prob = c(genes.diff.prob[x] * (1 - genes.diff.prob.up[x]),
                             1 - genes.diff.prob[x],
                             genes.diff.prob[x] * genes.diff.prob.up[x]))
  })
  
  effects <- (1.5 + dqrexp(N_GENES)) * (genes.dysregulated != 0)
  
  effects[genes.dysregulated == 0] <- 1
  effects[genes.dysregulated == -1] <- 1 / effects[genes.dysregulated == -1]
  
  tmp <- generateSyntheticData(experiments[experiment],
                               n.vars = N_GENES, samples.per.cond = experiment.samples[experiment],
                               effect.size = effects)
  
  tmp$sample.annotations$condition <- factor(tmp$sample.annotations$condition, labels = LETTERS[1:2])
  
  suppressMessages(DESeqDataSetFromMatrix(countData = tmp$count.matrix,
                                          colData = tmp$sample.annotations,
                                          design = ~ condition) %>% DESeq(quiet = T) %>% results %>%
                     as.data.table(keep.rownames = T) %>% .[, .(entrez.ID = rn, fc = log2FoldChange, adj.pv = padj)])
})

saveRDS(experiment.data, '/space/scratch/jsicherman/Thesis Work/data/artificial.rds')
source('/home/jsicherman/Thesis Work/main/load.R')

geneAssociations <- function(N) {
  # Divide all genes into discrete categories (age, behavior, biological process, biological sex, etc.)
  de.cat <- sample(DATA.HOLDER$human@experiment.meta[, .(sum(!is.na(cf.ValLongUri)),
                                                         sum(!is.na(cf.BaseLongUri))), cf.CatLongUri] %>%
                     .[V1 > 50 & V2 > 50, cf.CatLongUri], N, T)

  # Draw a gene function from any of the possibilities within each category
  data <- cbind(entrez.ID = paste0('g', 1:N), rbindlist(lapply(unique(de.cat), function(category) {
    vals <- unique(DATA.HOLDER$human@experiment.meta[cf.CatLongUri == category, .(cf.ValLongUri, cf.BaseLongUri)])
    indices <- sample(1:nrow(vals), sum(de.cat == category), T)

    data.table(cf.CatLongUri = category, vals[indices, ])
  })))
}

gene.assoc <- geneAssociations(N_GENES)

saveRDS(gene.assoc, '/space/scratch/jsicherman/Thesis Work/data/gene.associations.rds')

pipeFilter <- function(data, FUNC) {
  Filter(FUNC, data)
}

exp.assoc <- rbindlist(mclapply(experiment.data, function(experiment) {
  experiment <- experiment[!is.na(adj.pv)]
  # Each gene has an associated factor its upregulation is associated with (reverse for downregulation)

  # An experiment's tags can be driven by any DEG
  genes <- experiment[adj.pv < 0.05, .(entrez.ID, upregulated = fc > 0)]

  # Select between one and three genes whose tags will be used to generate this experiment's tags
  whichGenes <- genes$entrez.ID[sample(1:nrow(genes), sample(1:min(nrow(genes), 3), 1, prob = c(0.8, 0.15, 0.05)[1:min(nrow(genes), 3)]))]

  mData <- merge(genes[entrez.ID %in% whichGenes], gene.assoc[entrez.ID %in% whichGenes], by = 'entrez.ID')

  possibleTags <- c(as.character(mData$cf.ValLongUri), as.character(mData$cf.BaseLongUri))

  # Generate an experiment tag by picking between 1-3 derived from the chosen gene(s)' contrasts
  if(length(possibleTags) == 0)
    ee.TagLongUri <- NA
  else
    ee.TagLongUri <- DATA.HOLDER$human@experiment.meta[cf.ValLongUri %in% possibleTags | cf.BaseLongUri %in% possibleTags,
                                                       sample(ee.TagLongUri, sample(1:min(length(ee.TagLongUri), 3), 1,
                                                                                    prob = c(0.8, 0.15, 0.05)[1:min(length(ee.TagLongUri), 3)]))] %>%
    na.omit %>% as.character %>% pipeFilter(function(x) x != '') %>%  paste0(collapse = '; ')
  
  cf.ValLongUri <- lapply(whichGenes, function(gene) {
    tmp <- mData[entrez.ID == gene, ifelse(upregulated, as.character(cf.ValLongUri), as.character(cf.BaseLongUri))]
  }) %>% na.omit %>% pipeFilter(function(x) x != '' & x != 'NA') %>% paste0(collapse = '; ')

  cf.BaseLongUri <- lapply(whichGenes, function(gene) {
    tmp <- mData[entrez.ID == gene, ifelse(!upregulated, as.character(cf.ValLongUri), as.character(cf.BaseLongUri))]
  }) %>% na.omit %>% pipeFilter(function(x) x != '' & x != 'NA') %>% paste0(collapse = '; ')
  
  if(ee.TagLongUri == '') ee.TagLongUri <- NA
  if(cf.BaseLongUri == '') cf.BaseLongUri <- NA
  if(cf.ValLongUri == '') cf.ValLongUri <- NA

  data.table(cf.CatLongUri = mData[entrez.ID %in% whichGenes, sample(cf.CatLongUri, 1)],
             ee.TagLongUri = ee.TagLongUri,
             cf.ValLongUri = cf.ValLongUri,
             cf.BaseLongUri = cf.BaseLongUri,
             n.DE = experiment[adj.pv < 0.05, .N],
             mean.fc = experiment[, mean(fc, na.rm = T)])
}))

rm(pipeFilter)

experiment.meta <- data.table(rsc.ID = experiments,
                       ee.ID = 1:N_EXPERIMENTS,
                       ee.Name = experiments,
                       ee.Source = 'Artificial' %>% as.factor,
                       ee.NumSamples = experiment.samples,
                       ee.TagLongUri = exp.assoc[, ee.TagLongUri] %>% as.factor,
                       ee.qScore = rnorm(N_EXPERIMENTS, 0.2, 0.3) %>% pmin(1), # TODO This should be more meaningful
                       ad.Name = 'TODO' %>% as.factor,
                       ad.Company = 'TODO' %>% as.factor,
                       ad.Sequencing = T %>% as.factor,
                       sf.Subset = F,
                       sf.Cat = 'TODO' %>% as.factor,
                       sf.CatLongUri = 'TODO' %>% as.factor,
                       sf.ValLongUri = 'TODO' %>% as.factor,
                       cf.Cat = 'TODO' %>% as.factor,
                       cf.CatLongUri = exp.assoc[, cf.CatLongUri] %>% as.factor,
                       cf.ValLongUri = exp.assoc[, cf.ValLongUri] %>% as.factor,
                       cf.BaseLongUri = exp.assoc[, cf.BaseLongUri] %>% as.factor,
                       n.DE = exp.assoc[, n.DE],
                       mean.fc = exp.assoc[, mean.fc])

saveRDS(experiment.meta, '/space/scratch/jsicherman/Thesis Work/data/artificial.experiment.meta.rds')
