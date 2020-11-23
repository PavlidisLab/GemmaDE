options(app.name = 'EnriCh')
options(app.description = 'Pending')
options(app.tags = 'genomics,bioinformatics,genetics,transcriptomes,rnaseq,microarrays,biotechnology,medicine,biomedical,meta-analysis,statistics,search,open source,database,software')
options(app.author = 'Jordan Sicherman (jordan.sicherman@msl.ubc.ca)')
options(spinner.color = '#002145')
options(spinner.type = 6)

options(max.progress.steps = 7)
options(max.gemma = 500)
options(app.max.rows = 100)
options(app.min.tags = 2)
options(app.pv = 0.05)
options(app.fc_lower = 0, app.fc_upper = 10)
options(app.distance_cutoff = 2.25)
options(app.mfx = T, app.geeq = T)
options(app.search_method = 'mvsm')
options(app.taxa = 'human',
        app.ontology = c('GO', 'CLO', 'CL', 'CHEBI', 'DO', 'EFO', 'TGEMO', 'HP', 'MP', 'OBI', 'UBERON'))
options(app.all_taxa = list(`H. sapiens` = 'human', `M. musculus` = 'mouse', `R. norvegicus` = 'rat', artificial = 'artificial'))
options(app.all_search_methods = list(`M-VSM` = 'mvsm', `Weighted Average` = 'zscore'))

options(app.all_options = list(pv = getOption('app.pv'), fc.lower = getOption('app.fc_lower'),
                               fc.upper = getOption('app.fc_upper'), mfx = getOption('app.mfx'),
                               geeq = getOption('app.geeq'), distance = getOption('app.distance_cutoff'),
                               min.tags = getOption('app.min.tags'),
                               max.rows = getOption('app.max.rows'), method = getOption('app.search_method')))

source('main/process.R')
source('main/progress.R')
source('main/renderTools.R')
source('main/gemmaAPI.R')
source('main/load.R')
