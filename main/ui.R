library(shiny)
library(shinyjs)
library(data.table)
library(DT)

generateResultsHeader <- function(title) {
  fluidRow(class = 'info-text', column(12, h2(title)))
}

generateResults <- function(taxa = 'human', scope = 'DO', experiments, conditions, options = DEFAULT_OPTIONS) {
  filterResults <- function(conditions) {
    do.call(rbind, lapply(1:(nrow(conditions) - 1), function(indx) {
      if(conditions$N.ranked[indx] != conditions$N.ranked[indx + 1])
        conditions[indx, ]
    })) %>% as.data.table
  }
  
  if(options$filter && nrow(conditions) > 1)
    conditions <- filterResults(conditions)
  
  mCache <- CACHE.BACKGROUND[[scope]][rsc.ID %in% rownames(experiments) & Definition %in% conditions$Definition]
  eDF <- experiments %>% as.data.frame
  rownames(eDF) <- rownames(experiments)
  colnames(eDF) <- colnames(experiments)
  
  geneScores <- data.table(def = unique(conditions$Definition))
  
  for(indx in 1:(ncol(experiments) - 1)) {
    geneScores[, colnames(experiments)[indx] := mean(eDF[mCache[Definition == def, rsc.ID], indx], na.rm = T),
               def]
  }
  
  conditions <- merge(conditions, geneScores, by.x = 'Definition', by.y = 'def') %>%
    setnames('pv', 'P-value')

  outputColumns <- c('P-value', colnames(experiments)[1:(ncol(experiments) - 1)])
  renderDataTable(datatable(conditions[, outputColumns, with = F] %>% as.data.frame,
                            extensions = c('ColReorder', 'FixedHeader', 'Buttons'),
                            rownames = conditions$Definition,
                            filter = 'top',
                            options = list(colReorder = T,
                                           pageLength = 10,
                                           order = list(
                                             list(1, 'asc')),
                                           language = list(lengthMenu = 'Show _MENU_ conditions per page'),
                                           fixedHeader = T,
                                           rowCallback = JS('asScientificNotation'),
                                           dom = 'lBfrtip',
                                           searchCols = list(NULL, list(search = '0 ... 0.05')),
                                           buttons = list(
                                             list(extend = 'csvHtml5',
                                                  text = 'Download',
                                                  title = 'data')))))
}

ui <- fluidPage(style = 'height: 100%;',
  useShinyjs(),
  
  # Add to the HTML head
  withTags(head(
    script(src = 'js/initialize.js'),
    link(rel = 'stylesheet', type = 'text/css', href = 'css/style.css'),
    
    title(options('app.name')),
    meta(name = 'description', content = options('app.description')),
    meta(name = 'keywords', content = options('app.tags')),
    meta(name = 'author', content = options('app.author')),
    
    link(rel = 'icon', href = 'https://gemma.msl.ubc.ca/images/favicon.ico')
  )),
  
  # A wrapper for the page content.
  div(class = 'content',
      
    # A pretty header
    fluidRow(id = 'header'),
    
    # The gene entry bar
    fluidRow(wellPanel(
      # Search parameters
      fluidRow(
      # Gene entry
      column(5,
             fluidRow(selectInput('genes', 'Input Gene(s) of Interest', list(`Enter query...` = ''), multiple = T)),
             fluidRow(helpText(HTML('Examples: <a genes=\'["RPS4Y1","XIST","KDM5D"]\'>RPS4Y1, XIST, KDM5D</a>, <a genes=\'"ENSG00000121410"\'>ENSG00000121410</a>')))),
      column(2, uiOutput('genes.csv.ui')),
      
      # Taxa entry
      column(2, selectInput('taxa', 'Taxon', TAXA)),
      
      # Ontology entry (with more options, as it's on the right)
      column(3,
             selectizeInput('scope', 'Ontologies', ONTOLOGIES[, unique(OntologyScope)], selected = 'DO', multiple = T,
                            options = list(selectOnTab = T)),
             helpText(style = 'float: right;', HTML('<a data-toggle="collapse" data-target="#options">More options...</a>')))
      ),
      
      # More options
      wellPanel(class = 'collapse', id = 'options',
                fluidRow(style = 'display: flex; flex-direction: row;',
                         column(6, wellPanel(`well-name` = 'Filtering',
                           checkboxInput('filter', 'Filter equal-scoring children', value = DEFAULT_OPTIONS$filterSame))),
                         column(6, wellPanel(`well-name` = 'Scoring',
                           checkboxInput('mfx', 'Multifunctionality', value = DEFAULT_OPTIONS$mfx),
                           numericInput('pv', 'P-value threshold', value = DEFAULT_OPTIONS$pv, step = 0.01, min = 0, max = 1),
                           sliderInput('fc', 'FC threshold', value = c(DEFAULT_OPTIONS$fc.lower, DEFAULT_OPTIONS$fc.upper), step = 0.1, min = 0, max = 10, ticks = F)))
                )),
      
      # Buttons
      fluidRow(
        # Core action buttons
        column(12,
               tags$form(style = 'float: right;',
                         actionButton('reset', 'Reset', class = 'btn-secondary'),
                         actionButton('search', 'Search', class = 'btn-primary', disabled = NA)
               )
        )
      )
    )),
    
    fluidRow(
      column(12, htmlOutput('results_header')),
      column(12, dataTableOutput('results'))
    )
  ),
  
  HTML('
<footer>
  <div class="text-center">© 2020 Copyright: <a href="mailto:jordan.sicherman@msl.ubc.ca">Jordan Sicherman</a></div>
</footer>'
      )
)
