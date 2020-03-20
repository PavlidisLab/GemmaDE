library(shiny)
library(data.table)

ui <- fluidPage(
  
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
  
  # A pretty header
  fluidRow(id = 'header'),
  
  # The gene entry bar
  fluidRow(wellPanel(
    
    # Search parameters
    fluidRow(
      
    # Gene entry
    column(5,
           fluidRow(selectInput('genes', 'Input Gene(s) of Interest', list(`Enter query...` = ''), multiple = T)),
           fluidRow(helpText(HTML('Examples: <a genes="AKAP1">AKAP1</a>, <a genes="ENSG00000121410">ENSG00000121410</a>')))),
    column(2, uiOutput('genes.csv.ui')),
    
    # Taxa entry
    column(2, selectInput('taxa', 'Taxa', TAXA)),
    
    # Ontology entry
    column(3, selectInput('scope', 'Ontologies', ONTOLOGIES[, unique(OntologyScope)], selected = 'DO', multiple = T))
    ),
    
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
  ))
)
