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

source('dependencies.R')

runApp('main', port = 18232, launch.browser = F)

# ROADMAP
# [/] Get a good simulation framework set up for the new analyses
# [-] Fix single gene queries
# [-] Output gene-wise contributions to scoring
# [-] Consider interpolating between null distributions
# [/] Consider result caching
# [-] Release to lab
# [-] Think of names
# [x] Add a biology-related loader
# [-] Support multisessions
# [-] Release name poll

# FOR WRITING UP
# Sex, cell type, tissue specific findings
# Anecdotes (drug-related, metabolic alterations in PD astrocytes)
# Demonstration through simulations
# Breakdown of gaps in knowledge incl. counts for other tools (per species)
