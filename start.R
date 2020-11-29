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

# TODO consider memoise::memoise -ing things

# ssh -L 12345:localhost:18232 jsicherman@nelson.msl.ubc.ca -p 22000

# Sz genes from Lilah's
# Doping real data

# Technical details
# If we can set up lazy computation of Gemma linkages then we don't limit to limit rows returned

# Outstanding questions
# What to do with LINEAR, UNSCALED, OTHER, COUNT scales?
# Is ora only fair when the bg is only the detected genes?
