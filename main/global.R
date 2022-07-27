.START_TIME <- Sys.time()
source(here::here("main/requirements.R"))
source(here::here("main/dependencies.R"))
message(paste0("Finished loading in ", round(Sys.time() - .START_TIME, 2)))
