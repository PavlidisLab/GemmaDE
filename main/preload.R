library(igraph)
library(data.table)
library(dplyr)

lapply(c('OBI', 'UBERON', 'DO', 'HP', 'MP', 'CL'), function(scope) {
  getTags(scope = scope)
})
