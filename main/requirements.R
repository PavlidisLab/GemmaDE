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
library(d3wordcloud) # devtools::install_github('jbkunst/d3wordcloud')
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
library(gemmaAPI, lib.loc = '/home/omancarci/R/x86_64-redhat-linux-gnu-library/3.6/')
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
