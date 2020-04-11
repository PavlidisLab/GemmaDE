library(data.table)
library(shiny)
library(shinyjs)
library(DT)
library(heatmaply)
library(shinyHeatmaply)
library(shinyWidgets)
library(sparkline)
library(shinycssloaders)

generateResultsHeader <- function(title) {
  if(is.character(title))
    fluidRow(class = 'info-text', column(12, h2(title)))
  else
    fluidRow(class = 'info-text', column(12, title))
}

generateResultsPlot <- function(taxa = 'human', scope = 'DO', experiments, conditions, options = DEFAULT_OPTIONS, input) {
  if(length(unlist(input$plotData$selected)) > 0)
    conditions <- conditions[1 + unlist(input$plotData$selected), ]
  else if(input$plot_page)
    conditions <- conditions[1 + unlist(input$plotData$page), ]
  
  if(input$plot_type == 'Heatmap') {
    if(input$plot_data == 'Gene Score') {
      heatmaply_cor(conditions[, lapply(.SD, function(x) mean(as.numeric(unlist(strsplit(x, ','))))), by = Definition,
                               .SDcols = colnames(experiments)[1:ncol(experiments) - 1]
                               ][, colnames(experiments)[1:(ncol(experiments) - 1)], with = F] %>% as.data.frame %>% t,
                    limits = NULL, dendrogram = 'column', labCol = conditions$Definition)
    } else if(input$plot_data == 'P-value')
      heatmaply_cor(conditions[, `P-value (χ2)`] %>% as.data.frame %>% t, limits = NULL,
                dendrogram = 'column', labCol = conditions$Definition, labRow = 'P-value (χ2)', node_type = 'scatter',
                point_size_mat = -log10(conditions[, `P-value (χ2)`]) %>% as.data.frame %>% t, point_size_name = '-log10(P)')
  }
}

generateResults <- function(taxa = 'human', scope = 'DO', experiments, conditions, options = DEFAULT_OPTIONS) {
  outputColumns <- c('Evidence', 'E', 'P-value (χ2)', 'P-value (Fisher)', colnames(experiments)[1:(ncol(experiments) - 1)])
  
  conditions[, Evidence := paste0('<span data-toggle="popover" title="Experiments" data-html="true" data-content="',
                                  lapply(unlist(strsplit(Evidence, ',')), function(experiment) {
                                    paste0('<a target=_blank href=https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=',
                                           experiment, '>',
                                           DATA.HOLDER[[taxa]]@experiment.meta[ee.ID == experiment, unique(ee.Name)], '</a>')
                                  }) %>% paste0(collapse = ', '), '">', paste(length(unlist(strsplit(Evidence, ','))), 'Experiments', '<i class="fas fa-question-circle"></i>'), '</span>'), Definition]
  
  mTable <- datatable(conditions[, outputColumns, with = F] %>% as.data.frame,
                      extensions = c('FixedHeader', 'Buttons'),
                      rownames = conditions$Definition,
                      escape = -2,
                      filter = 'top',
                      options = list(pageLength = 10,
                                     order = list(
                                       list(3, 'asc')),
                                     language = list(lengthMenu = 'Show _MENU_ conditions per page',
                                                     processing = '',
                                                     emptyTable = 'No matching conditions found.',
                                                     infoEmpty = 'Showing 0 to 0 of 0 conditions',
                                                     info = 'Showing _START_ to _END_ of _TOTAL_ conditions',
                                                     infoFiltered = '(filtered from _MAX_ total conditions)'),
                                     fixedHeader = T,
                                     rowCallback = JS('asScientificNotation'),
                                     initComplete = JS('onTableCreated'),
                                     drawCallback = JS('onTableDraw'),
                                     dom = 'lBfrtip',
                                     searchCols = list(
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       list(search = '0 ... 0.05'),
                                                       list(search = '0 ... 0.05')),
                                     autoWidth = T,
                                     columnDefs = list(
                                      # list(targets = 0,
                                      #      width = '20%'),
                                       list(targets = 1,
                                            className = 'dt-right'),
                                       list(targets = 5:length(outputColumns),
                                            render = JS('asSparkline'), width = '1px', className = 'dt-center', searchable = F, orderable = F)
                                     ),
                                     search = list(
                                       list(regex = T)
                                     ),
                                     buttons = list(
                                       list(extend = 'collection',
                                            text = 'Visualize',
                                            action = JS('plotData')),
                                       list(extend = 'csvHtml5',
                                            text = 'Download',
                                            title = 'data',
                                            exportOptions = list(
                                              format = list(
                                                #format.header = JS('function(html, col, node) { return html; }'),
                                                #format.footer = JS('function(html, col, node) { return html; }'),
                                                body = JS('unformatSpark')
                                                #customizeData = JS('function(data) { console.log(data); }')
                                              )
                                            )
                                       )
                                     )
                                    )
                      )
  mTable$dependencies <- append(mTable$dependencies, htmlwidgets:::getDependency('sparkline'))
  
  renderDT(mTable)
}

ui <- fluidPage(style = 'height: 100%;',
  useShinyjs(),
  
  # Add to the HTML head
  withTags(head(
    script(src = 'js/initialize.js'),
    script(src = 'https://kit.fontawesome.com/33dcd9d8f9.js', crossorigin = 'anonymous'),
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
                           materialSwitch('filter', 'Filter equal-scoring children', value = DEFAULT_OPTIONS$filterSame, right = T))),
                         column(6, wellPanel(`well-name` = 'Scoring',
                           materialSwitch('mfx', 'Multifunctionality', value = DEFAULT_OPTIONS$mfx, right = T),
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
      column(12, dataTableOutput('results') %>% withSpinner)
    )
  ),
  
  HTML('
<footer>
  <div class="text-center">© 2020 Copyright: <a href="mailto:jordan.sicherman@msl.ubc.ca">Jordan Sicherman</a></div>
</footer>'
      )
)
