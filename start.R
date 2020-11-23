library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(htmlwidgets)
library(sparkline)
library(DT)
library(heatmaply)
library(shinyHeatmaply)
library(shinypanels)
library(RColorBrewer)

library(jsonlite)
library(async)

library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)

library(gemmaAPI)
library(ermineR)
library(mygene)

source('dependencies.R')

runApp('main', port = 18232, launch.browser = F)

# Check metadata for log transformation
# Internal weighting by baseline expression level

# TODO consider memoise::memoise -ing things

# ssh -L 12345:localhost:18232 jsicherman@nelson.msl.ubc.ca -p 22000
