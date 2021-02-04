options(app.name = 'TODO')
options(app.description = 'TODO')
options(app.tags = 'genomics,bioinformatics,genetics,transcriptomes,rnaseq,microarrays,biotechnology,medicine,biomedical,meta-analysis,statistics,search,open source,database,software')
options(app.author = 'Jordan Sicherman (jordan.sicherman@msl.ubc.ca)')

options(max.progress.steps = 3)
options(max.gemma = 1000)

addConfig <- function(description, category, extras = NULL, ...) {
  mList <- getOption('app.registered')
  mNew <- list(...)
  if(names(mNew) %in% unlist(mList)) return(NULL)
  
  mEntry <- list(name = names(mNew), value = unname(unlist(mNew)), category = category, description = description)
  if(!is.null(extras)) {
    for(mName in names(extras)) {
      mEntry[[mName]] <- extras[[mName]]
    }
  }
  
  if(is.null(mList))
    mList <- list(mEntry)
  else
    mList[[length(mList) + 1]] <- mEntry
  
  names(mList)[length(mList)] <- names(mNew)
  
  options(app.registered = mList)
}

getConfig <- function(key = NULL, category = NULL, ...) {
  if(is.null(key) && is.null(category))
    ret <- getOption('app.registered')
  else if(!is.null(key))
    ret <- getOption('app.registered')[[key]]
  else if(!is.null(category))
    ret <- Filter(function(x) x$category == category, getOption('app.registered'))
  
  if(length(list(...)) > 0) {
    for(v in names(list(...))) {
      ret[[v]]$value <- unname(list(...)[[v]])
    }
  }
  
  ret
}

as.input <- function(entry) {
  if(is.numeric(entry$value)) {
    if(length(entry$value) == 1)
      numericInput(entry$name, entry$description, value = entry$value,
                   min = entry$min, max = entry$max, step = entry$step)
    else
      sliderInput(entry$name, entry$description, value = entry$value,
                  min = entry$min, max = entry$max, step = entry$step, ticks = entry$ticks)
  } else if(is.character(entry$value))
    selectInput(entry$name, entry$description, entry$choices, entry$value)
  else if(is.logical(entry$value))
    materialSwitch(entry$name, entry$description, entry$value, right = T)
}

do.update <- function(sess, entry, val) {
  if(is.numeric(entry$value)) {
    if(length(entry$value) == 1)
      updateNumericInput(sess, entry$name, value = val)
    else
      updateSliderInput(sess, entry$name, value = val)
  } else if(is.character(entry$value))
    updateSelectizeInput(sess, entry$name, selected = val)
  else if(is.logical(entry$value))
    updateCheckboxInput(sess, entry$name, value = val)
}

options(app.registered = NULL)
addConfig(pv = 0.05, description = 'Significance threshold', category = 'Scoring', extras = list(min = 0, max = 1, step = 0.01))
addConfig(req = 1, description = 'Required DEG count', category = 'Scoring', extras = list(min = 1, max = NA, step = 1))
addConfig(fc = c(0, 100), description = 'Fold-change threshold', category = 'Filtering', extras = list(min = 0, max = 100, step = 1, ticks = F))
addConfig(mfx = T, description = 'Score multifunctionality', category = 'Scoring')
addConfig(geeq = T, description = 'Score experiment quality (GEEQ)', category = 'Scoring')
addConfig(method = 'diff', description = 'Scoring function', category = 'Scoring', extras = list(choices = list(`M-VSM` = 'mvsm', `Difference` = 'diff', `Correlation` = 'cor')))
addConfig(gemmaLink = F, description = 'Add links to Gemma', category = 'Filtering')

addConfig(taxa = 'human', description = NA, category = NA,
          extras = list(choices = list(`H. sapiens` = 'human',
                                       `M. musculus` = 'mouse',
                                       `R. norvegicus` = 'rat',
                                       `Any` = 'any',
                                       `artificial` = 'artificial'), # TODO remove for prod
                        core = c('human', 'mouse', 'rat'),
                        mapping = c(human = 9606, mouse = 10090, rat = 10116)))
addConfig(sig = '', description = NA, category = NA)

source('/home/jsicherman/Thesis Work/main/process.R')
source('/home/jsicherman/Thesis Work/main/renderTools.R')
source('/home/jsicherman/Thesis Work/main/gemmaAPI.R')
source('/home/jsicherman/Thesis Work/main/load.R')
