options(app.name = "Gemma DE")
options(app.description = "TODO")
options(app.tags = "genomics,bioinformatics,genetics,transcriptomes,rnaseq,microarrays,biotechnology,medicine,biomedical,meta-analysis,statistics,search,open source,database,software")
options(app.author = "Jordan Sicherman (jordan.sicherman@msl.ubc.ca)")

options(max.progress.steps = 3) # How many different calls will be made to advanceProgress
options(max.gemma = 1000) # The maximum number of samples to pull from Gemma for expression visualization
options(chunk.size = 200) # For artificial data generation; not used in live app

options(app.algorithm.gene.pre = quote(zScore * (1 - log10(Rfast::Pmax(matrix(1e-10, ncol = ncol(pv), nrow = nrow(pv)), pv)))))
options(app.algorithm.gene.post = quote(zScore %>% abs() %>% `*`(MFX_WEIGHT) %>% t() %>% data.table::as.data.table()))
options(app.algorithm.experiment = quote(ee.q * (1 + f.IN) / (1 + 10^f.OUT)))

future::plan(future::multicore, workers = 32) # Maximum clients that Gemma DE can serve

# The signature field can potentially have duplicate entries, which Shiny's underlying
# library (selectize.js) doesn't support. Kinda hacky way to replace the selectize dependency
# with tom select. Didn't bother perfecting this, seems to work well enough
customSelectizeIt <- function(inputId, select, options, nonempty = FALSE) {
  res <- checkAsIs(options)
  selectizeDep <- htmlDependency("tomselect", "1.7.4", c(href = "libs/tomselect"),
    stylesheet = "css/tom-select.bootstrap3.min.css", head = format(tagList(tags$script(src = "libs/tomselect/js/tom-select.complete.min.js")))
  )

  select$children[[2]] <- tagAppendChild(
    select$children[[2]],
    tags$script(
      type = "application/json", `data-for` = inputId,
      `data-nonempty` = if (nonempty) {
        ""
      }, `data-eval` = if (length(res$eval)) {
        HTML(toJSON(res$eval))
      }, if (length(res$options)) {
        HTML(toJSON(res$options))
      } else {
        "{}"
      }
    )
  )
  attachDependencies(select, selectizeDep)
}
environment(customSelectizeIt) <- asNamespace("shiny")
assignInNamespace("selectizeIt", customSelectizeIt, ns = "shiny")

# Creates UI options to configure algorithm
addConfig <- function(description, tooltip = "", category, extras = NULL, ...) {
  mList <- getOption("app.registered")
  mNew <- list(...)
  if (names(mNew) %in% unlist(mList)) {
    return(NULL)
  }

  mEntry <- list(name = names(mNew), value = unname(unlist(mNew)), category = category, description = description)
  if (!is.null(extras)) {
    for (mName in names(extras)) {
      mEntry[[mName]] <- extras[[mName]]
    }
  }

  mEntry[["tooltip"]] <- tooltip

  if (is.null(mList)) {
    mList <- list(mEntry)
  } else {
    mList[[length(mList) + 1]] <- mEntry
  }

  names(mList)[length(mList)] <- names(mNew)

  options(app.registered = mList)
}

getConfig <- function(key = NULL, category = NULL, ...) {
  if (is.null(key) && is.null(category)) {
    ret <- getOption("app.registered")
  } else if (!is.null(key)) {
    ret <- getOption("app.registered")[[key]]
  } else if (!is.null(category)) {
    ret <- Filter(function(x) x$category == category, getOption("app.registered"))
  }

  if (length(list(...)) > 0) {
    for (v in names(list(...))) {
      ret[[v]]$value <- unname(list(...)[[v]])
    }
  }

  ret
}

as.input <- function(entry) {
  as_input <- function() {
    if (is.numeric(entry$value)) {
      if (length(entry$value) == 1) {
        numericInput(entry$name, entry$description,
          value = entry$value,
          min = entry$min, max = entry$max, step = entry$step
        )
      } else {
        sliderInput(entry$name, entry$description,
          value = entry$value,
          min = entry$min, max = entry$max, step = entry$step, ticks = entry$ticks
        )
      }
    } else if (is.character(entry$value)) {
      if (isTRUE(entry$selectize)) {
        selectInput(entry$name, entry$description, entry$choices, entry$value, switch(is.null(entry$multiple) + 1,
          entry$multiple,
          F
        ))
      } else {
        pickerInput(entry$name, entry$description, entry$choices, entry$value, switch(is.null(entry$multiple) + 1,
          entry$multiple,
          F
        ))
      }
    } else if (is.logical(entry$value)) {
      checkboxInput(entry$name, entry$description, entry$value)
    }
  }

  span(as_input(), `data-toggle` = "tooltip", title = entry$tooltip)
}

do.update <- function(sess, entry, val) {
  if (is.numeric(entry$value)) {
    if (length(entry$value) == 1) {
      updateNumericInput(sess, entry$name, value = val)
    } else {
      updateSliderInput(sess, entry$name, value = val)
    }
  } else if (is.character(entry$value)) {
    updateSelectizeInput(sess, entry$name, selected = val)
  } else if (is.logical(entry$value)) {
    updateCheckboxInput(sess, entry$name, value = val)
  }
}

options(app.registered = NULL)
addConfig(pv = 0.05, description = "Significance threshold", tooltip = "The maximum FDR-corrected p-value to consider a gene differentially expressed", category = "Scoring", extras = list(min = 0, max = 1, step = 0.01))
addConfig(dist = 1.5, description = "Maximum term distance", tooltip = "The maximum distance each condition comparison can be from the annotated one. Larger values will take longer to search, but will allow for greater grouping", category = "Filtering", extras = list(min = 0, max = 10, step = 0.25))
addConfig(confounds = F, description = "Include possibly confounded", tooltip = "Whether or not to include factors that exhibit a possible batch confound", category = "Scoring")
# wasaddConfig(mfx = FALSE, description = "Score multifunctionality", tooltip = "Whether or not to weight each gene's contribution by its multifunctionality rank", category = "Scoring")
addConfig(geeq = F, description = "Score experiment quality (GEEQ)", tooltip = "Whether or not to weight each experiment's contribution by its GEEQ score", category = "Scoring")
# addConfig(method = "diff", description = "Scoring function", tooltip = "Which scoring algorithm to use. You should use the default unless you want to search for conditions that resemble a specific DE signature", category = "Scoring", extras = list(choices = list(`M-VSM` = "mvsm", `Default` = "diff", `Correlation` = "cor")))

addConfig(sig = "", description = NA, category = NA)

source(paste(PROJDIR, 'main/renderTools.R', sep='/'))
source(paste(PROJDIR, 'main/gemmaAPI.R', sep='/'))
source(paste(PROJDIR, 'main/load.R', sep='/'))

CORPUS_STATS <- lapply(DATA.HOLDER, function(x) x@experiment.meta) %>%
  data.table::rbindlist() %>%
  {
    mComparisons <- nrow(.)
    .[!duplicated(ee.ID), .(
      assays = sum(ee.NumSample),
      studies = max(.N),
      comparisons = mComparisons
    )] %>%
      .[, lapply(.SD, format, big.mark = ",")]
  }

# Make sure we can serve clients asynchronously
options(future.globals.maxSize = (object.size(DATA.HOLDER) + object.size(CACHE.BACKGROUND) * 1.5) %>% as.double())

mChoices <- list(
  `H. sapiens` = "human",
  `M. musculus` = "mouse",
  `R. norvegicus` = "rat",
  `artificial` = "artificial"
)
mChoices <- Filter(function(x) x %in% names(DATA.HOLDER), mChoices)

categories <- lapply(DATA.HOLDER, function(x) x@experiment.meta[, unique(cf.Cat)]) %>%
  unlist() %>%
  unique() %>%
  as.character() %>%
  sort() %>%
  {
    setNames(as.list(.), .)
  }

subsets <- lapply(DATA.HOLDER, function(x) x@experiment.meta[, unique(sf.Val)]) %>%
  unlist() %>%
  unique() %>%
  as.character() %>%
  {
    Filter(function(x) !(x %in% c("DE_Include", "DE_Exclude")), .)
  } %>%
  sort() %>%
  {
    setNames(as.list(.), .)
  }

addConfig(
  categories = Filter(function(x) {
    !(x %in% c(
      "block", "cell line", "collection of material",
      "developmental stage", "disease staging", "dose",
      "generation", "growth condition", "individual",
      "population", "protocol", "strain", "study design", "timepoint"
    ))
  }, names(categories)),
  description = "Categories to display", tooltip = "Which topics to include. Limiting these will speed up the computation at the cost of potentially missing significant results", category = "Filtering",
  extras = list(choices = categories, multiple = T)
)

# addConfig(
#   subset = NA_character_,
#   description = "Subsets to include/exclude", tooltip = "Limit which biomaterial subsets to include or exclude", category = "Filtering",
#   extras = list(choices = subsets, multiple = T, selectize = T)
# )
# addConfig(
#   subsetExcludes = T,
#   description = "Blacklist", "Toggle whether the selected subsets act as a whitelist or blacklist", category = "Filtering"
# )

addConfig(
  taxa = "human", description = NA, category = NA,
  extras = list(
    choices = mChoices,
    core = c("human", "mouse", "rat"),
    multiple = T,
    mapping = TAX.DATA$id %>% {names(.) = c("human", "mouse", "rat");.}
  )
)

rm(mChoices, categories, subsets)

letterWrap <- function(n, depth = 1) {
  x <- do.call(
    paste0,
    do.call(expand.grid, args = list(lapply(1:depth, function(x) {
      return(LETTERS)
    }), stringsAsFactors = F)) %>%
      .[, rev(names(.[])), drop = F]
  )

  if (n <= length(x)) {
    return(x[1:n])
  }

  return(c(x, letterWrap(n - length(x), depth = depth + 1)))
}

