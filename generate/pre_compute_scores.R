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
  
  
  
  # temporary until next generation, these also exist within compile.R
  dimnames(data.holder[[taxon]]@data$fc) %>%
    saveRDS(file.path(DATADIR, 'fbm', taxon, 'fc_contrast.dimnames.rds'))
  
  fc_contrast <- bigstatsr::as_FBM(data.holder[[taxon]]@data$fc,
                                   backingfile = file.path(DATADIR, 'fbm', taxon, 'fc_contrast'),
                                   is_read_only = T)$save()
  
  
  
  dimnames(data.holder[[taxon]]@data$zscore) %>% 
    saveRDS(file.path(DATADIR, 'fbm', taxon, 'z_contrast.dimnames.rds'))
  
  zscore_contrast <- bigstatsr::as_FBM(data.holder[[taxon]]@data$zscore,
                              backingfile = file.path(DATADIR, 'fbm', taxon, 'zscores_contrast'),
                              is_read_only = T)$save()
  
  
  dimnames(data.holder[[taxon]]@data$adj.pv) %>% 
    saveRDS(file.path(DATADIR, 'fbm', taxon, 'p_contrast.dimnames.rds'))
  
  adj.pv <- bigstatsr::as_FBM(data.holder[[taxon]]@data$adj.pv,
                              backingfile = file.path(DATADIR, 'fbm', taxon, 'adjpvs_contrast'),
                              is_read_only = T)$save()
  
  
  gc()
}


