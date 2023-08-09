.DATA_PATH <- paste(DATADIR, 'DATA.HOLDER.light.rds', sep='/')
fbm_path = file.path(DATADIR,'data_fbm')

# Load the lite versions if they're already created.
if (!exists("DATA.HOLDER")) {
  if (file.exists(.DATA_PATH)) {
    DATA.HOLDER <- readRDS(.DATA_PATH)
    DATA.HOLDER$artificial <- NULL
  } else {
    stop("Couldn't find DATA.HOLDER, run create_fbm.R first.")
  }
}

# Read existing FBMs
for (t in names(DATA.HOLDER)) {
  DATA.HOLDER[[t]]@data <- load_fbms(file.path(fbm_path,t),suffix = '_contrast')
}

gc()


