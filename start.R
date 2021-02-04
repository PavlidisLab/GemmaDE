# Visualizations
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders) # From jsicherman/shinycssloaders, NOT daattali
library(htmlwidgets)
library(DT)
library(heatmaply)
library(shinyHeatmaply)
library(shinypanels) # From jsicherman/shinypanels, NOT datasketch
library(circlepackeR)
library(d3wordcloud)
library(data.tree)
# library(sparkline)
library(RColorBrewer)
library(sass)

library(async)
library(memoise)

# Data drivers
library(matrixStats)
library(Rfast)
library(igraph)
library(dplyr)
library(data.table)
library(stringr)
library(bit)

# Parsing helpers
library(gemmaAPI, lib.loc = '/home/omancarci/R/x86_64-redhat-linux-gnu-library/3.6/')
library(ermineR)
library(mygene)
library(homologene)
library(jsonlite)
library(XML)

# Concurrent users
library(promises)
library(future)
plan(multicore, workers = 5)
options(future.globals.maxSize = 30 * 1000 * 1024^2)

source('dependencies.R')

runApp('main', port = 18232, launch.browser = F)

# ROADMAP
# [x] Add a biology-related loader

# [/] Support multisessions
# [ ]--- Needs to disable the search button
# [?]--- Cancel a process if the client disconnects

# [/] Get a good simulation framework set up for the new analyses
# [ ]--- Needs to be tested

# [-] Fix single gene queries

# [-] Output gene-wise contributions to scoring

# [-] Consider interpolating between null distributions

# [/] Consider result caching
# [ ]--- Make more decisions on what to cache

# [-] Release to lab
# [ ]--- Think of names
# [ ]--- Release name poll

# FOR WRITING UP
# Sex, cell type, tissue specific findings
# Anecdotes (drug-related, metabolic alterations in PD astrocytes)
# Demonstration through simulations
# Breakdown of gaps in knowledge incl. counts for other tools (per species)
