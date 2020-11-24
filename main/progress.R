advanceProgress <- function(detail) {
  n <- 1
  while((!exists('session', envir = parent.frame(n)) || is.null(parent.frame(n)$session)) && n < 10) { n <- n + 1 }
  
  if(exists('session', envir = parent.frame(n)))
    setProgress(parent.frame(n), parent.frame(n)$session$userData$progress + 1, detail)
  else
    warning('"session" didn\'t exist in calling chain.')
}

setProgress <- function(env, progress = NULL, detail = '', n.steps = NULL) {
  if(isTRUE(env)) {
    n <- 1
    while((!exists('session', envir = parent.frame(n)) || is.null(parent.frame(n)$session)) && n < 10) { n <- n + 1 }
    
    env <- parent.frame(n)
    if(!exists('session', envir = env))
      warning('"session" didn\'t exist in calling chain.')
  }
  
  if(!exists('session', envir = env)) {
    warning('"session" didn\'t exist in specified environment.')
    return(NULL)
  }
  
  if(is.null(progress)) {
    if(is.null(env$session$userData$progress.bar)) return(NULL)
    progress <- env$session$userData$progress.bar$getMax()
    n.steps <- progress
  }
  
  env$session$userData$progress <- progress
  
  if(progress == 0) {
    env$session$userData$progress.bar <- shiny::Progress$new(min = 0, max = n.steps)
    env$session$userData$progress.bar$set(message = 'Searching...', detail = detail)
  } else {
    env$session$userData$progress.bar$set(value = progress, detail = detail)
    
    if(progress >= env$session$userData$progress.bar$getMax())
      env$session$userData$progress.bar$close()
  }
}
