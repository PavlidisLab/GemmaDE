library(shiny)

source('main/load.R')
source('main/process.R')
source('main/server.R')
source('main/ui.R')

options('app.name' = "Jordan's Project")
options('app.description' = 'Pending')
options('app.tags' = 'genomics,bioinformatics,genetics,transcriptomes,rnaseq,microarrays,biotechnology,medicine,biomedical,meta-analysis,statistics,search,open source,database,software')
options('app.author' = 'Jordan Sicherman (jordan.sicherman@msl.ubc.ca)')

runApp('main', port = 18232, launch.browser = F)

# ssh -L 12345:localhost:18232 nelson -p 22000 -l jsicherman
