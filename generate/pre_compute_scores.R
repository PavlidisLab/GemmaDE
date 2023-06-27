# calculate score matrices twice to be able to access both gene and contrast level data

devtools::load_all()

data.holder = readRDS(file.path(RAWDIR, 'DATA.HOLDER.rds'))


for (tax in names(data.holder)){
  file.remove( file.path(DATADIR, 'fbm', tax, 'scores.bk'))
  file.remove( file.path(DATADIR, 'fbm', tax, 'score.dimnames.rds'))
  file.remove( file.path(DATADIR, 'fbm', tax, 'scores.rds'))
  
  file.remove( file.path(DATADIR, 'fbm', tax, 'scores_contrast.bk'))
  file.remove( file.path(DATADIR, 'fbm', tax, 'score_contrast.dimnames.rds'))
  file.remove( file.path(DATADIR, 'fbm', tax, 'scores_contrast.rds'))
  
  
  sc = abs(data.holder[[tax]]@data$zscore*(1 - log10(Rfast::Pmax(matrix(1e-10, 
                                                                    ncol = ncol(data.holder[[tax]]@data$adj.pv),
                                                                    nrow = nrow(data.holder[[tax]]@data$adj.pv)),
                                                             data.holder[[tax]]@data$adj.pv))))
  
  
  dimnames(sc) %>% saveRDS(file.path(DATADIR, 'fbm', tax, 'score_contrast.dimnames.rds'))
  sc_fb_contrast  = bigstatsr::as_FBM(sc,
                                      backingfile = file.path(DATADIR, 'fbm', tax, 'scores_contrast'),
                                      is_read_only = T)$save()
  
  
  dimnames(sc) %>% 
    rev() %>%
    saveRDS(file.path(DATADIR, 'fbm', tax, 'score.dimnames.rds'))
  
  sc_fb = bigstatsr::as_FBM(sc %>% t,
                    backingfile = file.path(DATADIR, 'fbm', tax, 'scores'),
                    is_read_only = T)$save()
  gc()
}


