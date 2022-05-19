source('/home/jsicherman/Thesis Work/requirements.R')
source('/home/jsicherman/Thesis Work/dependencies.R')

library(parallel)
library(edgeR)
library(compcodeR)
options(mc.cores = 5)
USE <- 'simulated'

rm(NULLS, CACHE.BACKGROUND, ONTOLOGIES, ONTOLOGIES.DEFS)
DATA.HOLDER[c('artificial', 'mouse', 'rat')] <- NULL

mu.phi.estimates <- system.file("extdata", "Pickrell.Cheung.Mu.Phi.Estimates.rds",
                                package = "compcodeR")
mu.phi.estimates <- readRDS(mu.phi.estimates)

generateSyntheticData <- function (dataset, n.vars, samples.per.cond, n.diffexp, repl.id = 1, 
                                   seqdepth = 1e+07, minfact = 0.7, maxfact = 1.4, relmeans = "auto", 
                                   dispersions = "auto", fraction.upregulated = 1, between.group.diffdisp = FALSE, 
                                   filter.threshold.total = 0, filter.threshold.mediancpm = 0, 
                                   fraction.non.overdispersed = 0, random.outlier.high.prob = 0, 
                                   random.outlier.low.prob = 0, single.outlier.high.prob = 0, 
                                   single.outlier.low.prob = 0, effect.size = 1.5, output.file = NULL) 
{
  if (!is.null(output.file)) {
    if (!(substr(output.file, nchar(output.file) - 3, nchar(output.file)) == 
          ".rds")) {
      stop("output.file must be an .rds file.")
    }
  }
  uID <- paste(sample(c(0:9, letters, LETTERS), 10, replace = TRUE), 
               collapse = "")
  condition <- rep(c(1, 2), each = samples.per.cond)
  S1 <- 1:samples.per.cond
  S2 <- 1:samples.per.cond + samples.per.cond
  if (length(effect.size) == 1) {
    n.upregulated <- floor(fraction.upregulated * n.diffexp)
    if (fraction.upregulated != 0 & n.diffexp != 0) {
      genes.upreg <- 1:n.upregulated
    }
    else {
      genes.upreg <- NULL
    }
    if (fraction.upregulated != 1 & n.diffexp != 0) {
      genes.downreg <- (n.upregulated + 1):n.diffexp
    }
    else {
      genes.downreg <- NULL
    }
    genes.nonreg <- setdiff(1:n.vars, union(genes.upreg, 
                                            genes.downreg))
  }
  else {
    if (length(effect.size) != n.vars) {
      stop("The length of the effect.size vector must be the same as the number of simulated genes.")
    }
    else {
      genes.upreg <- which(effect.size > 1)
      genes.downreg <- which(effect.size < 1)
      genes.nonreg <- which(effect.size == 1)
      n.upregulated <- length(genes.upreg)
      n.diffexp <- length(genes.upreg) + length(genes.downreg)
      fraction.upregulated <- n.upregulated/n.diffexp
    }
  }
  differential.expression <- rep(0, n.vars)
  differential.expression[genes.upreg] <- 1
  differential.expression[genes.downreg] <- 1
  upregulation <- rep(0, n.vars)
  upregulation[genes.upreg] <- 1
  downregulation <- rep(0, n.vars)
  downregulation[genes.downreg] <- 1
  if (is.character(relmeans) | is.character(dispersions)) {
    mu.estimates <- mu.phi.estimates$pickrell.cheung.mu
    phi.estimates <- mu.phi.estimates$pickrell.cheung.phi
    to.include <- sample(1:length(mu.estimates), n.vars, 
                         replace = ifelse(n.vars > length(mu.estimates), TRUE, 
                                          FALSE))
    truedispersions.S1 <- phi.estimates[to.include]
    truemeans.S1 <- mu.estimates[to.include]
  }
  if (!is.character(relmeans)) {
    if (length(relmeans) != n.vars) 
      stop("The length of the relmeans vector must be the same as the number of simulated genes.")
    truemeans.S1 <- c(relmeans)
  }
  if (!is.character(dispersions)) {
    if (nrow(cbind(dispersions)) != n.vars) 
      stop("The number of provided dispersions must be the same as the number of simulated genes.")
    truedispersions.S1 <- cbind(dispersions)[, 1]
    if (ncol(cbind(dispersions)) > 1) {
      truedispersions.S2 <- cbind(dispersions)[, 2]
    }
    else {
      truedispersions.S2 <- truedispersions.S1
    }
  }
  nfacts <- runif(2 * samples.per.cond, min = minfact, max = maxfact)
  seq.depths <- nfacts * seqdepth
  overdispersed <- rep(1, n.vars)
  if (fraction.non.overdispersed > 0) {
    overdispersed[genes.upreg[1:round(fraction.non.overdispersed * 
                                        length(genes.upreg))]] <- 0
    overdispersed[genes.downreg[1:round(fraction.non.overdispersed * 
                                          length(genes.downreg))]] <- 0
    overdispersed[genes.nonreg[1:round(fraction.non.overdispersed * 
                                         length(genes.nonreg))]] <- 0
  }
  prob.S1 <- truemeans.S1
  prob.S2 <- rep(0, length(prob.S1))
  if (length(effect.size) == 1) {
    for (i in 1:n.vars) {
      if (i %in% genes.upreg) {
        prob.S2[i] <- (effect.size + rexp(1, rate = 1)) * 
          prob.S1[i]
      }
      else {
        if (i %in% genes.downreg) {
          prob.S2[i] <- 1/(effect.size + rexp(1, rate = 1)) * 
            prob.S1[i]
        }
        else {
          prob.S2[i] <- prob.S1[i]
        }
      }
    }
  }
  else {
    prob.S2 <- c(effect.size) * prob.S1
  }
  true.log2foldchange <- log2(prob.S2/prob.S1)
  sum.S1 <- sum(prob.S1)
  sum.S2 <- sum(prob.S2)
  if (is.character(dispersions)) {
    truedispersions.S2 <- truedispersions.S1
    if (between.group.diffdisp == TRUE) {
      for (i in 1:length(truedispersions.S2)) {
        sample.base <- phi.estimates[abs(log10(mu.estimates) - 
                                           log10(prob.S2[i])) < 0.05]
        if (length(sample.base) < 50) {
          sample.base <- phi.estimates[order(abs(log10(mu.estimates) - 
                                                   log10(prob.S2[i])))][1:500]
        }
        truedispersions.S2[i] <- sample(sample.base, 
                                        1)
      }
    }
  }
  truedispersions.S1 <- truedispersions.S1 * overdispersed
  truedispersions.S2 <- truedispersions.S2 * overdispersed
  do.call(cbind, lapply(1:(length(S1) + length(S2)), function(j) {
    if(j %in% S1) {
      rnbinom(n.vars, mu = prob.S1/sum.S1 * 
                seq.depths[j], size = 1/truedispersions.S1)
    } else {
      rnbinom(n.vars, mu = prob.S2/sum.S2 * 
                seq.depths[j], size = 1/truedispersions.S2)
    }
  })) -> Z
  
  random.outliers <- matrix(0, nrow(Z), ncol(Z))
  random.outliers.factor <- matrix(1, nrow(Z), ncol(Z))
  if (random.outlier.high.prob != 0 | random.outlier.low.prob != 0) {
    tmp <- matrix(runif(prod(dim(Z))), nrow = nrow(Z))
    for (i in 1:nrow(Z)) {
      for (j in 1:ncol(Z)) {
        if (tmp[i, j] < random.outlier.high.prob) {
          random.outliers[i, j] <- 1
          random.outliers.factor[i, j] <- runif(1, min = 5, 
                                                max = 10)
        }
        else if (tmp[i, j] < random.outlier.low.prob + random.outlier.high.prob) {
          random.outliers[i, j] <- (-1)
          random.outliers.factor[i, j] <- 1/runif(1, 
                                                  min = 5, max = 10)
        }
      }
    }
    Z <- round(random.outliers.factor * Z)
  }
  
  has.single.outlier <- rep(0, n.vars)
  single.outliers <- matrix(0, nrow(Z), ncol(Z))
  single.outliers.factor <- matrix(1, nrow(Z), ncol(Z))
  if (single.outlier.high.prob != 0 | single.outlier.low.prob !=  0) {
    has.single.outlier[sample(1:n.vars, floor((single.outlier.high.prob +
                                                 single.outlier.low.prob) * n.vars))] <- 1
    
    seeds <- runif(nrow(Z))
    samples <- sample(1:ncol(Z), nrow(Z), T)
    factors <- runif(nrow(Z), min = 5, max = 10)
    
    rout <- lapply(1:nrow(Z), function(i) {
      if(has.single.outlier[i] == 1) {
        if(seeds[i] < single.outlier.high.prob / (single.outlier.high.prob + single.outlier.low.prob)) {
          list(c(rep(0, samples[i] - 1), 1, rep(0, ncol(Z) - samples[i])),
               c(rep(0, samples[i] - 1), factors[i], rep(0, ncol(Z) - samples[i])))
        } else {
          list(c(rep(0, samples[i] - 1), -1, rep(0, ncol(Z) - samples[i])),
               c(rep(0, samples[i] - 1), 1 / factors[i], rep(0, ncol(Z) - samples[i])))
        }
      } else
        list(rep(0, ncol(Z)), rep(1, ncol(Z)))
    })
    
    single.outliers <- do.call(rbind, lapply(rout, '[[', 1))
    single.outliers.factor <- do.call(rbind, lapply(rout, '[[', 2))
    
    Z <- round(single.outliers.factor * Z)
  }
  
  rownames(Z) <- 1:n.vars
  n.random.outliers.up.S1 <- apply(random.outliers[, S1] > 0, 1, sum)
  n.random.outliers.up.S2 <- apply(random.outliers[, S2] >  0, 1, sum)
  n.random.outliers.down.S1 <- apply(random.outliers[, S1] < 0, 1, sum)
  n.random.outliers.down.S2 <- apply(random.outliers[, S2] < 0, 1, sum)
  n.single.outliers.up.S1 <- apply(single.outliers[, S1] > 0, 1, sum)
  n.single.outliers.up.S2 <- apply(single.outliers[, S2] > 0, 1, sum)
  n.single.outliers.down.S1 <- apply(single.outliers[, S1] < 0, 1, sum)
  n.single.outliers.down.S2 <- apply(single.outliers[, S2] < 0, 1, sum)
  
  nf <- calcNormFactors(Z)
  norm.factors <- nf * colSums(Z)
  common.libsize <- exp(mean(log(colSums(Z))))
  pseudocounts <- sweep(Z + 0.5, 2, norm.factors, "/") * common.libsize
  log2.pseudocounts <- log2(pseudocounts)
  M.value <- apply(log2.pseudocounts[, S2], 1, mean) - apply(log2.pseudocounts[, S1], 1, mean)
  A.value <- 0.5 * (apply(log2.pseudocounts[, S2], 1, mean) + 
                      apply(log2.pseudocounts[, S1], 1, mean))
  variable.annotations <- data.frame(truedispersions.S1 = truedispersions.S1, 
                                     truedispersions.S2 = truedispersions.S2, truemeans.S1 = prob.S1, 
                                     truemeans.S2 = prob.S2, n.random.outliers.up.S1 = n.random.outliers.up.S1, 
                                     n.random.outliers.up.S2 = n.random.outliers.up.S2, n.random.outliers.down.S1 = n.random.outliers.down.S1, 
                                     n.random.outliers.down.S2 = n.random.outliers.down.S2, 
                                     n.single.outliers.up.S1 = n.single.outliers.up.S1, n.single.outliers.up.S2 = n.single.outliers.up.S2, 
                                     n.single.outliers.down.S1 = n.single.outliers.down.S1, 
                                     n.single.outliers.down.S2 = n.single.outliers.down.S2, 
                                     M.value = M.value, A.value = A.value, truelog2foldchanges = true.log2foldchange, 
                                     upregulation = upregulation, downregulation = downregulation, 
                                     differential.expression = differential.expression)
  rownames(variable.annotations) <- rownames(Z)
  sample.annotations <- data.frame(condition = condition, depth.factor = nfacts)
  info.parameters <- list(n.diffexp = n.diffexp, fraction.upregulated = fraction.upregulated, 
                          between.group.diffdisp = between.group.diffdisp, filter.threshold.total = filter.threshold.total, 
                          filter.threshold.mediancpm = filter.threshold.mediancpm, 
                          fraction.non.overdispersed = fraction.non.overdispersed, 
                          random.outlier.high.prob = random.outlier.high.prob, 
                          random.outlier.low.prob = random.outlier.low.prob, single.outlier.high.prob = single.outlier.high.prob, 
                          single.outlier.low.prob = single.outlier.low.prob, effect.size = effect.size, 
                          samples.per.cond = samples.per.cond, repl.id = repl.id, 
                          dataset = dataset, uID = uID, seqdepth = seqdepth, minfact = minfact, 
                          maxfact = maxfact)
  s <- apply(Z, 1, sum)
  keep.T <- which(s >= filter.threshold.total)
  Z.T <- Z[keep.T, ]
  variable.annotations.T <- variable.annotations[keep.T, ]
  filtering <- paste("total count >=", filter.threshold.total)
  cpm <- sweep(Z.T, 2, apply(Z.T, 2, sum), "/") * 1e+06
  m <- apply(cpm, 1, median)
  keep.C <- which(m >= filter.threshold.mediancpm)
  Z.TC <- Z.T[keep.C, ]
  variable.annotations.TC <- variable.annotations.T[keep.C, ]
  filtering <- paste(filtering, "; ", paste("median cpm >=", 
                                            filter.threshold.mediancpm))
  rownames(Z.TC) <- paste("g", 1:nrow(Z.TC), sep = "")
  colnames(Z.TC) <- paste("sample", 1:ncol(Z.TC), sep = "")
  rownames(sample.annotations) <- colnames(Z.TC)
  rownames(variable.annotations.TC) <- rownames(Z.TC)
  data.object <- compData(count.matrix = Z.TC, variable.annotations = variable.annotations.TC, 
                          sample.annotations = sample.annotations, filtering = filtering, 
                          info.parameters = info.parameters)
  if (!is.null(output.file)) {
    saveRDS(data.object, file = output.file)
  }
  return(invisible(data.object))
}

eMeta <- DATA.HOLDER$human@experiment.meta %>% copy
N_GENES <- nrow(DATA.HOLDER$human@gene.meta)
N_EXP <- nrow(eMeta)
CONTRASTS <- eMeta[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>% unique
EXP_CONTRASTS <- sample(1:nrow(CONTRASTS), N_EXP, T, sort(rexp(nrow(CONTRASTS)), T))

rm(DATA.HOLDER)

r1exp <- function(n, n.1, val.1 = 2, val.1.rate = 1, rate = 1) {
  sample(c(val.1 + rexp(n.1, val.1.rate), val.1 * rexp(n - n.1, rate) %>% `/`(max(.))))
}

CONTRAST_AFFINITY <- lapply(unique(EXP_CONTRASTS), function(contrast) {
  data.table(contrast = contrast,
             entrez.ID = 1:N_GENES,
             probability = r1exp(N_GENES, n.1 = 10)) %>%
    .[, effect := 1.5 + probability / 2] %>%
    .[sample(c(T, F), N_GENES, T, c(0.5012399, 0.4987578)), effect := 1 / effect]
}) %>% rbindlist

saveRDS(CONTRAST_AFFINITY, '/space/scratch/jsicherman/Thesis Work/data/artificial/contrast_aff.rds')

mclapply(1:N_EXP, function(experiment) {
  message(paste0(Sys.time(), ' ... ', round(100 * experiment / N_EXP, 2), '%'))
  
  # Sample some genes to DE
  if(eMeta$n.DE[experiment] > 0) {
    mDat <- CONTRAST_AFFINITY[contrast == EXP_CONTRASTS[experiment]] %>%
      setorder(-probability) %>%
      head(floor(1.5 * eMeta$n.DE[experiment])) %>%
      .[sample(1:nrow(.), eMeta$n.DE[experiment], F, .$probability)] %>%
      .[, .(entrez.ID, effect)]
  } else
    mDat <- data.table(entrez.ID = NA_character_, effect = NA_real_)
  
  # Make the others have effect of 1
  mDat <- data.table(entrez.ID = 1:N_GENES,
                     effect = 1) %>% merge(mDat, by = 'entrez.ID', all.x = T) %>%
    .[!is.na(effect.y), effect.x := effect.y] %>%
    .[, c('effect', 'effect.x', 'effect.y') := list(effect.x, NULL, NULL)]
  
  outliers <- 0.05 * 10^(-eMeta$ee.qScore[experiment])
  outliers[is.na(outliers)] <- 0.5
  
  tmp <- generateSyntheticData('',
                               n.vars = N_GENES, samples.per.cond = max(2, floor(eMeta[experiment, ee.NumSample] / 2)),
                               effect.size = mDat$effect,
                               single.outlier.high.prob = outliers/2,
                               single.outlier.low.prob = outliers/2)
  tmp@sample.annotations$condition <- factor(tmp@sample.annotations$condition, labels = LETTERS[1:2])
  
  if(sum(mDat[, effect != 1]) > 0)
    MISSING <- (1:N_GENES)[-mDat[effect != 1, entrez.ID]]
  else
    MISSING <- 1:N_GENES
  
  tmp@variable.annotations$truelog2foldchanges[
    sample(MISSING, N_GENES - eMeta$ad.NumGenes[experiment])] <- NA_real_
  
  if(USE == 'DESeq2') {
    suppressMessages(DESeqDataSetFromMatrix(countData = tmp@count.matrix,
                                            colData = tmp@sample.annotations,
                                            design = ~ condition) %>% DESeq(quiet = T) %>% results %>%
                       as.data.table(keep.rownames = T) %>% {
                         list(contrast = EXP_CONTRASTS[experiment],
                              DEGs = mDat[effect != 1, entrez.ID],
                              entrez.ID = .[, as.integer(substring(rn, 2))],
                              fc = .[, log2FoldChange],
                              adj.pv = .[, padj])
                       })
  } else if(USE == 'Limma') {
    d0 <- DGEList(tmp@count.matrix) %>% calcNormFactors()
    cutoff <- 1
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    d <- d0[-drop, ]
    mm <- model.matrix(~0 + condition, tmp@sample.annotations)
    yy <- voom(d, mm, plot = F)
    fit <- lmFit(yy, mm)
    contr <- makeContrasts(conditionB - conditionA, levels = colnames(coef(fit)))
    contrasts.fit(fit, contr) %>% eBayes %>% topTable(sort.by = 'P', number = Inf) %>%
      as.data.table(keep.rownames = T) %>% {
        list(contrast = EXP_CONTRASTS[experiment],
             DEGs = mDat[effect != 1, entrez.ID],
             entrez.ID = .[, as.integer(substring(rn, 2))],
             fc = .[, logFC],
             adj.pv = .[, adj.P.Val])
      }
  } else {
    list(contrast = EXP_CONTRASTS[experiment],
         DEGs = mDat[effect != 1, entrez.ID],
         fc = tmp@variable.annotations$truelog2foldchanges + rnorm(N_GENES, sd = 0.05),
         adj.pv = runif(N_GENES) %>% {
           a <- .
           spikes <- as.integer(mDat[, 1.05 * sum(effect != 1)])
           a[sample(1:length(a), spikes)] <- runif(spikes, 1e-32, 0.05)
           a <- p.adjust(a, method = 'BH')
           a[tmp@variable.annotations$differential.expression == 1] <- runif(tmp@info.parameters$n.diffexp, 1e-32, 0.05)
           a[is.na(tmp@variable.annotations$truelog2foldchanges)] <- NA_real_
           a
         })
  }
}) %>% saveRDS('/space/scratch/jsicherman/Thesis Work/data/artificial/experiment_dat.rds')

rm(mu.phi.estimates, generateSyntheticData, eMeta, N_GENES,
   N_EXP, CONTRASTS, EXP_CONTRASTS, r1exp, CONTRAST_AFFINITY, USE)
source('/home/jsicherman/Thesis Work/generate/postprocessArtificial.R')
