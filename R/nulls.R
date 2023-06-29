# use a matrix with contrasts as columns to calculate the null distribution
generate_nulls = function(scores_contrast = "score_contrast",na_omit = TRUE ,max_cores = 8){
  tax = DATA.HOLDER %>% names
  
  parallel::mclapply(tax, function(t){
    dt = DATA.HOLDER[[t]]@data[[scores_contrast]]
    seq_len(ncol(dt)) %>% lapply(function(i){
      print(i)
      vals = dt[,i]
      if(na_omit){
        vals %<>% na.omit()
      } else{
        vals[is.na(vals)]=0
      }
      data.table(rn = dimnames(dt)[[2]][i],
                 score.mean = mean(vals),
                 score.sd = sd(vals))
    }) %>% do.call(rbind,.)  -> nulls
    
  },mc.cores = min(length(tax),max_cores)
  )-> out
  
  names(out) = tax
  return(out)
}
