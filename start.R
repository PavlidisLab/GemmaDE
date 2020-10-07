library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(htmlwidgets)
library(sparkline)
library(DT)
library(heatmaply)
library(shinyHeatmaply)

library(jsonlite)
library(async)

library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)

library(gemmaAPI)

source('dependencies.R')

runApp('main', port = 18232, launch.browser = F)

# ssh -L 12345:localhost:18232 jsicherman@nelson.msl.ubc.ca -p 22000
