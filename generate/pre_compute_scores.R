# calculate score matrices twice to be able to access both gene and contrast level data

devtools::load_all()

data.holder = readRDS(file.path(RAWDIR, 'DATA.HOLDER_fixed_terms.rds'))
fbm_path = file.path(DATADIR,'data_fbm')


for (tax in names(data.holder)){
  
  sc = abs(data.holder[[tax]]@data$zscore*(1 - log10(Rfast::Pmax(matrix(1e-10,
                                                                    ncol = ncol(data.holder[[tax]]@data$adj.pv),
                                                                    nrow = nrow(data.holder[[tax]]@data$adj.pv)),
                                                             data.holder[[tax]]@data$adj.pv))))
  
  create_fbm(sc,path = file.path(fbm_path,tax,'score'))
  
  
  gc()
}


