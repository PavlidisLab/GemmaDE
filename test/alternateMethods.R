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
  n.nooverlap <- n.numgenes - n.overlap
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

search.genechaser(DATA.HOLDER$human@gene.meta[gene.Name %in% c('NANOG', 'POU5F1', 'SOX2', 'LIN28A', 'LIN28B'), entrez.ID], getConfig(filter = mFilter)) %>%
  enrich.genechaser %>%
  .[score.up > 0 | score.down > 0]

search.genechaser(DATA.HOLDER$mouse@gene.meta[gene.Name %in% c('Nanog', 'Pou5f1', 'Sox2', 'Lin28a', 'Lin28b'), entrez.ID],
                  getConfig(filter = DATA.HOLDER$mouse@experiment.meta$ad.Type != 'GENELIST', taxa = 'mouse')) %>%
  enrich.genechaser(getConfig(taxa = 'mouse')) %>%
  .[score.up > 0 | score.down > 0]

DATA.HOLDER$human@gene.meta[gene.Name %in% strsplit('KLF9, PER1, TSC22D3, ZBTB16', ', ')[[1]], entrez.ID] %>%
  search.genechaser() %>% enrich.genechaser() %>% .[score.up > 0] %>% setorder(-score.up) %>% print

res <- c()
mGenes <- DATA.HOLDER$human@gene.meta[gene.Name %in% (strsplit('ALDOA, ARL3, ATP2A2, CA8, CACNA1C, CACNB2, CHRNA3, CHRNA5, CHRNB4, CKB, CLCN3, CYP17A1, DPP4, DPYD, DRD2, FES, GRIA1, GRIN2A, GRM3, IREB2, KLC1, LRP1, MAN2A2, MEF2C, MMP16, MYO1A, NAB2, NEK1, FURIN, PSMA4, PTPRF, RRAS, SHMT2, STAT6, TAC3, TAF5, TCF4, TLE3, VRK2, XRCC3, CDK2AP1, MAD1L1, CUL3, INA, BAG5, PLCH2, KDM4A, TRANK1, MPHOSPH9, FUT9, NXPH4, R3HDM2, CNKSR2, NT5C2, PDCD11, RIMS1, ABCB9, NGEF, GIGYF2, EPC2, REEP2, ARL6IP4, CNNM2, WBP1L, FANCL, SBNO1, AS3MT, PITPNM2, ZSWIM6, SCAF1, SLC39A8, PJA1, BCL11B, ZFYVE21, GDPD3, OGFOD2, ACTR5, TRIM8, IMMP2L, PCGF6, APOPT1, ESAM, C12orf65, TRMT61A, SFXN2, TYW5, SLC32A1, CNTN4, RILPL2, TSNARE1, STAC3, ASPHD1, SNX19, MIR137', ', ')[[1]]), entrez.ID]
for (i in 1:1000) {
  print(i)
  res <- unique(c(res, sample(mGenes, 10) %>% search.genechaser() %>% enrich.genechaser() %>%
                    .[score.up > 0 | score.down > 0, paste0(c(as.character(cf.Cat), as.character(cf.Baseline), as.character(cf.Val)), collapse = '')]))
}
