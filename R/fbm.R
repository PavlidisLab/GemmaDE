# save an fbm of a data to a given directory path.
# it saves two copies of the data, one transposed along
# with their dimnames
create_fbm = function(data, path){
  file.remove(list.files(path,full.names = TRUE))
  dir.create(path,showWarnings = FALSE, recursive = TRUE)
  
  dimnames(data) %>% saveRDS(file.path(path,'dimnames.rds'))
  dimnames(data) %>% rev %>% saveRDS(file.path(path,'t.dimnames.rds'))
  
  bigstatsr::as_FBM(data,
                    backingfile = file.path(path, 'dat'),
                    is_read_only = T)$save()
  
  bigstatsr::as_FBM(data %>% t,
                    backingfile = file.path(path, 't.dat'),
                    is_read_only = T)$save()
  
}



# returns a list with all fbms within an fbm_path loaded
# use the parent path of inout of create_fbm
# use a suffix to differentiate between the regular and transposed
# forms of the data
load_fbms <-  function(path,suffix = '_contrast',t_suffix = ''){
  objects <- list.files(path)
  
  objects %>% lapply(function(x){
    out = list()
    
    out[[paste0(x,suffix)]] =  bigstatsr::big_attach(file.path(path,x,'dat.rds'))
    attr(out[[paste0(x,suffix)]], ".dimnames") <- readRDS(file.path(path,x,'dimnames.rds'))
     
    out[[paste0(x,t_suffix)]] =  bigstatsr::big_attach(file.path(path,x,'t.dat.rds'))
    attr(out[[paste0(x,t_suffix)]], ".dimnames") <- readRDS(file.path(path,x,'t.dimnames.rds'))
    
    return(out)
  }) %>% do.call(c,.) -> out
  
  return(out)
}



dimnames.FBM <- function(object, ...) {
  attr(object, ".dimnames")
}
