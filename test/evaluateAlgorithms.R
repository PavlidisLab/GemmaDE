updateAlgorithm <- function(ID, EID) {
  if(ID == 1) {
    options(app.algorithm.gene.pre = quote(zScore * (1 - log10(Rfast::Pmax(matrix(1e-10, ncol = ncol(pv), nrow = nrow(pv)), pv)))))
    options(app.algorithm.gene.post = quote(zScore %>% abs %>% `*`(MFX_WEIGHT) %>% t %>% as.data.table))
  } else if(ID == 2) {
    options(app.algorithm.gene.pre = quote(pv <= 0.05))
    options(app.algorithm.gene.post = quote(zScore %>% `*`(MFX_WEIGHT) %>% t %>% as.data.table))
  } else if(ID == 3) {
    options(app.algorithm.gene.pre = quote(zScore %>% abs %>% `*`(1 - log10(Rfast::Pmax(matrix(1e-10, ncol = ncol(pv), nrow = nrow(pv)), pv)))))
    options(app.algorithm.gene.post = quote(zScore %>% t %>% as.data.table))
  } else if(ID == 4) {
    options(app.algorithm.gene.pre = quote(zScore %>% abs %>% pv <= 0.05))
    options(app.algorithm.gene.post = quote(zScore %>% t %>% as.data.table))
  } else if(ID == 5) {
    options(app.algorithm.gene.pre = quote(pv <= 0.05))
    options(app.algorithm.gene.post = quote(zScore %>% t %>% as.data.table))
  } else if(ID == 6) {
    #options(app.algorithm.gene.pre = quote(pv <= 0.05))
    #options(app.algorithm.gene.post = quote(zScore %>% t %>% as.data.table))
  }
  
  if(EID == 1)
    options(app.algorithm.experiment = quote(ee.q * (1 + f.IN) / (1 + 10^f.OUT)))
  else if(EID == 2)
    options(app.algorithm.experiment = quote((1 + f.IN) / (1 + 10^f.OUT)))
  else if(EID == 3)
    options(app.algorithm.experiment = quote(ee.q))
  else if(EID == 4)
    options(app.algorithm.experiment = quote(1))
}

for(g in 1:6) {
  for(e in 1:4) {
    updateAlgorithm(g, e)
    
    # TODO Need to do bootstrap first
  }
}
