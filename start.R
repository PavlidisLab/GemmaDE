library(shiny)

source('main/process.R')
source('main/server.R')
source('main/load.R')
source('main/ui.R')

options(app.name = 'EnriCh')
options(app.description = 'Pending')
options(app.tags = 'genomics,bioinformatics,genetics,transcriptomes,rnaseq,microarrays,biotechnology,medicine,biomedical,meta-analysis,statistics,search,open source,database,software')
options(app.author = 'Jordan Sicherman (jordan.sicherman@msl.ubc.ca)')
options(spinner.color = '#002145')
options(spinner.type = 6)

options(max.progress.steps = 8)
options(max.rows = 500)
options(app.pv = 0.05)
options(app.fc_lower = 0, app.fc_upper = 10)
options(app.distance_cutoff = 2.25)
options(app.mfx = T, app.geeq = T)
options(app.taxa = 'human', app.ontology = 'DO')
options(app.all_taxa = list(`H. sapiens` = 'human', `M. musculus` = 'mouse', `R. norvegicus` = 'rat', artificial = 'artificial'))

options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                               fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                               geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff')))

runApp('main', port = 18232, launch.browser = F)

# ssh -L 12345:localhost:18232 jsicherman@nelson.msl.ubc.ca -p 22000

# sbatch -C thrd64 "/home/jsicherman/Thesis Work/test/analyze_background.sh"
