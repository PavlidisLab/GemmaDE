library(shiny)
library(shinyjs)
library(data.table)
library(DT)

generateResultsHeader <- function(title) {
  fluidRow(class = 'info-text', h2(title))
}

generateResults <- function(experiments, conditions, options = DEFAULT_OPTIONS) {
  conditions <- conditions[`scaled N` > options$score.lower]
  if(options$score.lower != options$score.upper)
    conditions <- conditions[`scaled N` < options$score.upper]
  
  renderUI({
    lapply(0:min(nrow(conditions), options$n.display), function(row) {
      if(row == 0)
        return(generateResultsHeader(paste0('Found ', length(experiments), ' related experiment', ifelse(length(experiments) == 1, '', 's'), '.')))
      
      name.control <- gsub('( |-|\\+|,|\\.|:|;|\\/|\'|"|&|\\(|\\))', '', tolower(conditions$Definition[row]))
      div(class = 'results-row',
        fluidRow(
          HTML(paste0('<a data-toggle="collapse" aria-expanded="false" aria-controls="', name.control, '"',
                      '" class = "result-line" href="#', name.control, '">',
                      paste0(conditions$Definition[row], '  (', round(conditions$`scaled N`[row], 2), ')'), '</a>'))
      ),
      
      fluidRow(id = name.control, class = 'collapse result-info',
        p(paste('Enrichment Score (ES):', conditions$count[row] / conditions$background[row])),
        p(paste('Rank-adjusted Enrichment Score (RES):', conditions$`ranked N`[row])),
        p(paste('Scaled Rank-adjusted Enrichment Score (SRES):', conditions$`scaled N`[row]))
      )
    )})
  })
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
      column(2, selectInput('taxa', 'Taxa', TAXA)),
      
      # Ontology entry (with more options, as it's on the right)
      column(3,
             selectInput('scope', 'Ontologies', ONTOLOGIES[, unique(OntologyScope)], selected = 'DO', multiple = T),
             helpText(style = 'float: right;', HTML('<a data-toggle="collapse" data-target="#options">More options...</a>')))
      ),
      
      # More options
      wellPanel(class = 'collapse', id = 'options',
                fluidRow(style = 'display: flex; flex-direction: row;',
                         column(4, wellPanel(`well-name` = 'Filtering',
                           numericInput('top-n', 'Max to display', value = DEFAULT_OPTIONS$n.display, min = 0, max = 100),
                           sliderInput('score', 'Score threshold', value = c(DEFAULT_OPTIONS$score.lower, DEFAULT_OPTIONS$score.upper), step = 1, min = 0, max = 500, ticks = F))),
                         column(4, wellPanel(`well-name` = 'Scoring',
                           checkboxInput('mfx', 'Multifunctionality', value = DEFAULT_OPTIONS$mfx),
                           numericInput('pv', 'P-value threshold', value = DEFAULT_OPTIONS$pv, step = 0.01, min = 0, max = 1),
                           sliderInput('fc', 'FC threshold', value = c(DEFAULT_OPTIONS$fc.lower, DEFAULT_OPTIONS$fc.upper), step = 0.1, min = 0, max = 10, ticks = F))),
                         column(4, wellPanel(`well-name` = 'TODO',
                           textInput('some_text', 'Something', placeholder = 'Temporary'),
                           textInput('more_text', 'More Stuff', placeholder = 'Also Temporary')))
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
    
    fluidRow(htmlOutput('results', class = 'container-fluid'))
  ),
  
  HTML('
<footer>
  <div class="text-center">© 2020 Copyright: <a href="mailto:jordan.sicherman@msl.ubc.ca">Jordan Sicherman</a></div>
</footer>'
      )
)
