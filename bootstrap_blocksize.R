# source(paste(PROJDIR, 'main/requirements.R', sep='/'))
print('bootstrap data')
devtools::load_all()
library(parallel)
RNGkind("L'Ecuyer-CMRG")

options(mc.cores = 14)

ITERS <- 1000
BLOCK <- c(1,5,10,100,250,500,100)


OPTIONS <- c( 'human', 'mouse', 'rat')


for (b in BLOCK){
  dir.create(file.path(DATADIR,'blocksize_nulls',b),recursive = TRUE,showWarnings = FALSE)
  for(x in OPTIONS) {
    rm(DATA.HOLDER)
    source(paste(PROJDIR, 'main/dependencies.R', sep='/'))
    DATA.HOLDER[OPTIONS[-which(OPTIONS == x)]] <- NULL
    message(paste0(Sys.time(), ' ... Starting ', x))
    mclapply(1:ITERS, function(j) {
      set.seed(j)
      if(j %% (ITERS / 20) == 0)
        message(paste0(Sys.time(), ' ... ', x, ' ... ', round(100 * j / ITERS, 2), '%'))
      
      opts = getConfig(taxa = x)
      
      DATA.HOLDER[[x]]@gene.meta[sample(1:nrow(DATA.HOLDER[[x]]@gene.meta), b), entrez.ID] %>%
        vsmSearch(taxa = opts$taxa$value, 
                  confounds = TRUE, 
                  filter = opts$filter$value, 
                  mfx = FALSE, 
                  geeq = opts$geeq$value,
                  p_threshold = 0.05) %>%
        .[, c(seq_len(b), b+5),with = FALSE] %>% # TODO BLOCK
        data.table::melt(id.vars = 'rn') %>%
        .[, !'variable'] %>%
        .[, .(score.sum = sum(value, na.rm = T),
              score.sqsum = sum(value^2, na.rm = T)), rn]
    }) %>% data.table::rbindlist() %>% {
      message('Saving...')
      
      .[, .(score.mean = sum(score.sum) / (b * ITERS),
            score.sd = sqrt(sum(score.sqsum) / (b * ITERS) - (sum(score.sum) / (b * ITERS))^2)), rn] %>%
        saveRDS(file.path(DATADIR,'blocksize_nulls',b,paste0(x, '.rds')))
      
      gc()
      NULL
    }
  }
}
