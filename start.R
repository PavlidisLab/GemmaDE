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
library(stringr)
library(bit)
library(circlepackeR)
library(data.tree)

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
library(homologene)
library(XML)

source('dependencies.R')

runApp('main', port = 18232, launch.browser = F)

# TODO consider memoise::memoise -ing things

# Sz genes from Lilah's

# Technical details
# If we can set up lazy computation of Gemma linkages then we don't limit to limit rows returned
