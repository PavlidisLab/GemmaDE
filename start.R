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

#Sys.setenv(PATH = paste('/home/jsicherman/miniconda2/bin:/home/jsicherman/miniconda2/condabin',
#                        Sys.getenv('PATH'), sep = ':'))
#Sys.setenv(LD_LIBRARY_PATH = paste('/home/jsicherman/centos/usr/lib64',
#                                   Sys.getenv('LD_LIBRARY_PATH'), sep = ':'))

source('dependencies.R')

runApp('main', port = 18232, launch.browser = F)

# ROADMAP
# [x] Add a biology-related loader

# [x] Support multisessions
# [x]--- Needs to disable the search button
# [x]--- Cancel a process if the client disconnects

# [/] Get a good simulation framework set up for the new analyses
# [ ]--- Needs to be tested

# [x] Fix single gene queries
# [x] Fix plot saving
# [x]--- Seems overly complicated and takes awhile
# [x]--- Changed from orca to plotly built-in svg
# [x] Fix table saving

# [/] Output gene-wise contributions to scoring
# [ ]--- These numbers may not be very meaningful as they only communicate amount of DE. If they could somehow
#        be modified to portray specificity, it would be ideal
# [ ]--- Need to have a condition selector to minimize legend overhead

# [-] Consider interpolating between null distributions

# [x] Consider result caching
# [ ]--- Make more decisions on what to cache

# [-] Release to lab
# [ ]--- Think of names
# [ ]--- Release name poll

# FOR WRITING UP
# Sex, cell type, tissue specific findings
# Anecdotes (drug-related, metabolic alterations in PD astrocytes)
# Demonstration through simulations
# Breakdown of gaps in knowledge incl. counts for other tools (per species)
