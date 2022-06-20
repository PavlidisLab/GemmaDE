# Visualizations
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders) # From jsicherman/shinycssloaders, NOT daattali devtools::install_github('jsicherman/shinycssloaders')
library(htmlwidgets)
library(DT)
library(heatmaply)
library(shinyHeatmaply)
library(shinypanels) # From jsicherman/shinypanels, NOT datasketch devtools::install_github('jsicherman/shinypanels')
library(circlepackeR) # devtools::install_github('jeromefroe/circlepackeR')
library(data.tree)
library(RColorBrewer)
library(sass)
library(shinydisconnect)

library(async) # devtools::install_github('gaborcsardi/async')
library(memoise)

# Data drivers
library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)
library(stringr)
library(bigstatsr)
# library(bit)
library(matrixTests)

# Parsing helpers
library(gemmaAPI, lib.loc = "/home/omancarci/R/x86_64-redhat-linux-gnu-library/3.6/")
library(ermineR) # devtools::install_github('PavlidisLab/ermineR')
library(mygene) # BiocManager::install("mygene")
library(homologene)
library(jsonlite)
library(XML)
library(sass)
library(stringdist)

# Concurrent users
library(promises)
library(future)

#' JSONify
#'
#' Assuming the input string is entirely unquoted, make it suitable for JSON.
jsonify <- function(str) {
  if (is.null(str) || length(grep("(\\[|\\])", str, value = F)) == 0) {
    return(str)
  }
  
  gsub("\\[", '["', gsub("\\]", '"]', gsub(",", '","', gsub(", ", ",", str)))) %>% parse_json(simplifyVector = T)
}

#' Parse List Entry
#'
#' Given an entry (that may or may not be NA) of semicolon delimited (; ) list (in which, any may be NA),
#' return a vector of those values, coercing NA strings into true NAs.
#'
#' @param entry NA or a character of semicolon delimited (; ) values.
#'
#' @return A character vector of values in the input character.
parseListEntry <- function(entry) {
  if (length(entry) == 0) {
    return(NA)
  }
  
  if (is.na(entry) || !grepl("; ", entry, fixed = T)) {
    return(entry)
  }
  
  cleaned <- strsplit(entry %>% as.character(), "; ", fixed = T) %>% unlist()
  ifelse("NA" == cleaned, NA, cleaned)
}
