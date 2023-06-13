.DATA_PATH <- paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/')


# Load the lite versions if they're already created.
if (!exists("DATA.HOLDER")) {
  if (file.exists(.DATA_PATH)) {
    DATA.HOLDER <- readRDS(.DATA_PATH)
    DATA.HOLDER$artificial <- NULL
  } else {
    stop("Couldn't find DATA.HOLDER, run generate/fbm.R first.")
  }
}

# Read existing FBMs
for (i in names(DATA.HOLDER)) {
  DATA.HOLDER[[i]]@data$zscore <- bigstatsr::big_attach(paste0(paste(DATADIR, 'fbm/', sep='/'), i, "/zscores.rds"))
  attr(DATA.HOLDER[[i]]@data$zscore, ".dimnames") <- readRDS(paste0(paste(DATADIR, 'fbm/', sep='/'), i, "/z.dimnames.rds"))

  DATA.HOLDER[[i]]@data$adj.pv <- bigstatsr::big_attach(paste0(paste(DATADIR, 'fbm/', sep='/'), i, "/adjpvs.rds"))
  attr(DATA.HOLDER[[i]]@data$adj.pv, ".dimnames") <- readRDS(paste0(paste(DATADIR, 'fbm/', sep='/'), i, "/p.dimnames.rds"))
  
  score_file = paste0(file.path(DATADIR, 'fbm/'), i, "/scores.rds")
  if(file.exists(score_file)){
    DATA.HOLDER[[i]]@data$scores <- bigstatsr::big_attach(score_file)
    attr(DATA.HOLDER[[i]]@data$adj.pv, ".dimnames") <- readRDS(paste0(file.path(DATADIR, 'fbm/'), i, "/score.dimnames.rds"))
  }
  
  fc_file = paste0(file.path(DATADIR, 'fbm/'), i, "/fc.rds")
  if(file.exists(score_file)){
    DATA.HOLDER[[i]]@data$fc <- bigstatsr::big_attach(fc_file)
    attr(DATA.HOLDER[[i]]@data$fc, ".dimnames") <- readRDS(paste0(file.path(DATADIR, 'fbm/'), i, "/fc.dimnames.rds"))
  }
  
}

gc()

dimnames.FBM <- function(object, ...) {
  attr(object, ".dimnames")
}


