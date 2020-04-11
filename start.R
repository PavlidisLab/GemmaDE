library(shiny)

source('main/process.R')
source('main/load.R')
source('main/server.R')
source('main/ui.R')

options(app.name = 'EnriCh')
options(app.description = 'Pending')
options(app.tags = 'genomics,bioinformatics,genetics,transcriptomes,rnaseq,microarrays,biotechnology,medicine,biomedical,meta-analysis,statistics,search,open source,database,software')
options(app.author = 'Jordan Sicherman (jordan.sicherman@msl.ubc.ca)')
options(spinner.color = '#002145')
options(spinner.type = 6)

runApp('main', port = 18232, launch.browser = F)

# ssh -L 12345:localhost:18232 jsicherman@nelson -p 22000

# sbatch -C thrd64 "/home/jsicherman/Thesis Work/test/analyze_background.sh"
