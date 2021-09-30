source('main/process.R')
assign('search.default', search)
assign('enrich.default', enrich)

search <- function(genes, options = getConfig(), DATA = NULL) UseMethod('search')
enrich <- function(genes, options = getConfig(), CACHE = NULL) UseMethod('search')

search.genechaser <- function(genes, options = getConfig(), DATA = NULL) {
  if(is.null(DATA))
    DATA <- DATA.HOLDER
  
  mData <- DATA[[options$taxa$value]]
  
  experimentMask <- rep(T, nrow(mData@experiment.meta))
  
  if(!is.null(options$filter$value))
    experimentMask <- experimentMask & options$filter$value
  
  # Only retain GOI
  geneMask <- which(mData@gene.meta$entrez.ID %in% genes)
  
  n.genes <- length(geneMask)
  if(n.genes == 0)
    return(NULL)
  
  # P-values for only the GOI
  pv <- mData@data$adj.pv[geneMask, experimentMask, drop = F]
  fc <- mData@data$fc[geneMask, experimentMask, drop = F]
  
  pv[is.na(pv)] <- 1
  
  data.table(rn = colnames(pv),
             passing = Rfast::colsums(pv <= options$pv$value) == n.genes & colSums2(abs(fc > 1), na.rm = T) == n.genes) %>%
    .[, c('score.up', 'score.down') := list(passing * colMeans2(fc * ifelse(fc > 0, 1, NA_integer_), na.rm = T),
                                            passing * -colMeans2(fc * ifelse(fc < 0, 1, NA_integer_), na.rm = T))] %>%
    `class<-`(c('genechaser', 'data.table', 'data.frame'))
}

search.gxa <- function(genes, options = getConfig(), DATA = NULL) {
  if(is.null(DATA))
    DATA <- DATA.HOLDER
  
  mData <- DATA[[options$taxa$value]]
  
  experimentMask <- rep(T, nrow(mData@experiment.meta))
  
  if(!is.null(options$filter$value))
    experimentMask <- experimentMask & options$filter$value
  
  # Only retain GOI
  geneMask <- which(mData@gene.meta$entrez.ID %in% genes)
  
  n.genes <- length(geneMask)
  if(n.genes == 0)
    return(NULL)
  
  # P-values for only the GOI
  pv <- mData@data$adj.pv[geneMask, experimentMask, drop = F]
  pv[is.na(pv)] <- 1
  
  n.overlap <- Rfast::colsums(pv < 0.05)
  n.DE <- mData@experiment.meta[experimentMask == T] %>% .[, n.DE]
  n.numgenes <- mData@experiment.meta[experimentMask == T] %>% .[, ad.NumGenes]
  n.other <- n.DE - n.overlap
  n.nooverlap <- n.genes - n.overlap
  n.other2 <- n.numgenes - n.DE - n.genes + n.overlap
  
  mNames <- colnames(pv)
  
  mclapply(1:ncol(pv), function(i) {
    data.table(rn = mNames[i],
               pv = tryCatch(fisher.test(matrix(c(n.overlap[i], n.other[i], n.nooverlap[i], n.other2[i]), byrow = T, nrow = 2), alternative = 'greater')$p.value,
                             error = function(e) 1),
               observed = n.overlap[i],
               expected = n.genes * n.DE[i] / n.numgenes[i])
  }) %>% rbindlist %>% .[, padj := p.adjust(pv, 'BH')] %>% .[, effectsize := observed / expected] %>%
    .[is.finite(effectsize)] %>%
    setorder(-effectsize, pv) %>%
    `class<-`(c('gxa', 'data.table', 'data.frame'))
}

enrich.genechaser <- function(rankings, options = getConfig(), CACHE = NULL) {
  terms <- DATA.HOLDER[[options$taxa$value]]@experiment.meta[rsc.ID %in% rankings$rn, .(rsc.ID, cf.Cat, cf.Baseline, cf.Val)]
  
  terms <- rankings %>%
    merge(terms, by.x = 'rn', by.y = 'rsc.ID', sort = F) %>%
    .[cf.Cat %in% options$categories$value]
}

mFilter <- DATA.HOLDER$human@experiment.meta$ad.Type != 'GENELIST'
