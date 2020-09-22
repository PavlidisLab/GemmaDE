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

generateResultsPlot <- function(taxa = getOption('app.taxa'), scope = getOption('app.ontology'), experiments, conditions, options = DEFAULT_OPTIONS, input) {
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

generateResults <- function(data, taxa = getOption('app.taxa'), scope = getOption('app.ontology'), experiments, conditions, options = DEFAULT_OPTIONS) {
  # outputColumns <- c('Evidence', 'E', 'P-value (χ2)', 'P-value (Fisher)', colnames(experiments)[1:(ncol(experiments) - 1)])
  outputColumns <- c('Evidence', 'P-value (χ2)', 'P-value (Fisher)', colnames(experiments)[1:(ncol(experiments) - 1)])
  
  conditions[, Evidence := paste0('<span data-toggle="popover" title="Experiments" data-html="true" data-content="',
                                  lapply(unlist(strsplit(Evidence, ',')), function(experiment) {
                                    paste0('<a target=_blank href=https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=',
                                           experiment, '>',
                                           data@experiment.meta[ee.ID == experiment, unique(ee.Name)], '</a>')
                                  }) %>% paste0(collapse = ', '), '">', paste(N, paste0('Experiment', ifelse(N > 1, 's', '')), '<i class="fas fa-question-circle" style="cursor: pointer;"></i>'), '</span>'),
             .(cf.BaseLongUri, cf.ValLongUri)]
  
  mTable <- datatable(conditions[, outputColumns, with = F] %>% as.data.frame,
                      extensions = c('FixedHeader', 'Buttons'),
                      rownames = paste0('<b>', conditions$cf.BaseLongUri, '</b> vs. <b>', conditions$cf.ValLongUri, '</b>'),
                      escape = -c(1, 2),
                      filter = 'top',
                      options = list(pageLength = 10,
                                     order = list(
                                       list(2, 'asc')),
                                     language = list(lengthMenu = 'Show _MENU_ conditions per page',
                                                     processing = '',
                                                     emptyTable = 'No matching conditions found.',
                                                     infoEmpty = 'Showing 0 to 0 of 0 condition-comparisons',
                                                     info = 'Showing _START_ to _END_ of _TOTAL_ condition-comparisons',
                                                     infoFiltered = '(filtered from over _MAX_ total)'),
                                     fixedHeader = T,
                                     rowCallback = JS('asScientificNotation'),
                                     initComplete = JS('onTableCreated'),
                                     drawCallback = JS('onTableDraw'),
                                     dom = 'lBfrtip',
                                     searchCols = list(
                                                       NULL,
                                                       NULL,
                                                       list(search = '0 ... 0.05'),
                                                       list(search = '0 ... 0.05')),
                                     autoWidth = T,
                                     columnDefs = list(
                                       list(targets = 0,
                                            width = '25%'),
                                       list(targets = 1,
                                            className = 'dt-right',
                                            searchable = F, orderable = F),
                                       list(targets = 2:3,
                                            width = '8%'),
                                       list(targets = 4:length(outputColumns),
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
    
    title(getOption('app.name')),
    meta(name = 'description', content = getOption('app.description')),
    meta(name = 'keywords', content = getOption('app.tags')),
    meta(name = 'author', content = getOption('app.author')),
    
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
      column(2, selectInput('taxa', 'Taxon', getOption('app.all_taxa'))),
      
      # Ontology entry (with more options, as it's on the right)
      column(3,
             selectizeInput('scope', 'Ontologies', ONTOLOGIES[, unique(OntologyScope)], selected = getOption('app.ontology'), multiple = T,
                            options = list(selectOnTab = T)),
             helpText(style = 'float: right;', HTML('<a data-toggle="collapse" data-target="#options">More options...</a>')))
      ),
      
      # More options
      wellPanel(class = 'collapse', id = 'options',
                fluidRow(style = 'display: flex; flex-direction: row;',
                         column(6, wellPanel(`well-name` = 'Filtering',
                           numericInput('distance', 'Ontology expansion limit', value = getOption('app.distance_cutoff'), step = 0.25, min = 0, max = 10))),
                         column(6, wellPanel(`well-name` = 'Scoring',
                           materialSwitch('mfx', 'Multifunctionality', value = getOption('app.mfx'), right = T),
                           materialSwitch('geeq', 'Experiment quality', value = getOption('app.geeq'), right = T),
                           numericInput('pv', 'Significance threshold', value = getOption('app.pv'), step = 0.01, min = 0, max = 1),
                           sliderInput('fc', 'FC threshold', value = c(getOption('app.fc_lower'), getOption('app.fc_upper')), step = 0.1, min = 0, max = 10, ticks = F)))
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
