# source(paste(PROJDIR, 'main/requirements.R', sep='/'))
print('bootstrap data')
devtools::load_all()
library(parallel)
RNGkind("L'Ecuyer-CMRG")

options(mc.cores = 14)

ITERS <- 1000
BLOCK <- 500


# removed artificial from options. the current code never saves it to data 
# holder light version and don't think it's used in downstream for the app -ogan
OPTIONS <- c( 'human', 'mouse', 'rat')

# x = 'human'
dir.create(file.path(DATADIR,'updated_nulls'),showWarnings = FALSE, recursive = TRUE)

for(x in OPTIONS) {
  rm(DATA.HOLDER)
  source(paste(PROJDIR, 'main/dependencies.R', sep='/'))
  DATA.HOLDER[OPTIONS[-which(OPTIONS == x)]] <- NULL
  message(paste0(Sys.time(), ' ... Starting ', x))
  # note the hard skip presumably added to force overwriting at some point - ogan
  if(F && file.exists(file.path(DATADIR, 'updated_nulls',paste0(x, '.rds')))) {
    message(paste0('File for ', x, ' already exists... Skipping.'))
  } else {
    mclapply(1:ITERS, function(j) {
      set.seed(j)
      if(j %% (ITERS / 20) == 0)
        message(paste0(Sys.time(), ' ... ', x, ' ... ', round(100 * j / ITERS, 2), '%'))
      
      opts = getConfig(taxa = x)
      
      DATA.HOLDER[[x]]@gene.meta[sample(1:nrow(DATA.HOLDER[[x]]@gene.meta), BLOCK), entrez.ID] %>%
        vsmSearch(taxa = opts$taxa$value, 
                  confounds = TRUE, 
                  filter = opts$filter$value, 
                  mfx = FALSE, 
                  geeq = opts$geeq$value,
                  p_threshold = 0.05) %>%
        .[, c(1:500, 505)] %>% # TODO BLOCK
        data.table::melt(id.vars = 'rn') %>%
        .[, !'variable'] %>%
        .[, .(score.sum = sum(value, na.rm = T),
              score.sqsum = sum(value^2, na.rm = T)), rn]
    }) %>% data.table::rbindlist() %>% {
      message('Saving...')
      
      .[, .(score.mean = sum(score.sum) / (BLOCK * ITERS),
            score.sd = sqrt(sum(score.sqsum) / (BLOCK * ITERS) - (sum(score.sum) / (BLOCK * ITERS))^2)), rn] %>%
        saveRDS(file.path(DATADIR, 'updated_nulls',paste0(x, '.rds')))
      
      gc()
      NULL
    }
  }
}
