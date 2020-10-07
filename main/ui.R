generateResultsHeader <- function(title) {
  if('html' %in% class(title))
    fluidRow(class = 'info-text', column(12, title))
  else
    fluidRow(class = 'info-text', column(12, h2(title)))
}

generateResultsPlot <- function(genes, conditions, expr, options = getOption('app.all_options'), input, session) {
  if(is.null(input$plot_genes)) {
    updatePickerInput(session, 'plot_genes', choices = as.list(intersect(genes, rownames(expr$expr))),
                      selected = intersect(genes, rownames(expr$expr)))
    updatePickerInput(session, 'plot_conditions', choices = as.list(conditions), selected = conditions)
  }
  
  expr$expr <- expr$expr[input$plot_genes, ]
  
  if(length(input$plot_genes) == 1)
    expr$expr <- t(expr$expr)
  
  thresh <- 1.5 * apply(expr$expr, 1, iqr, na.rm = T)
  quants <- rowQuantiles(expr$expr, probs = c(0.25, 0.75), na.rm = T)
  
  inliers <- suppressWarnings((expr$expr < (quants[, 1] - thresh) | expr$expr > (quants[, 2] + thresh)) %>%
    apply(2, max, na.rm = T)) %>% {
      names(.[. == 0])
    }
  
  expr$metadata <- expr$metadata[name %in% inliers] %>% setorder(baseline)
  expr$expr <- expr$expr[, expr$metadata$name]
  
  if(length(input$plot_genes) == 1)
    expr$expr <- t(expr$expr)
  
  if(any(dim(expr$expr) == 0))
    return(NULL)
  
  if(input$plot_type == 'Heatmap') {
    if(input$plot_data == 'Gene Expression')
      heatmaply(expr$expr, labRow = NULL, showticklabels = c(F, T), scale = 'row', Rowv = NULL, dendrogram = 'none', col_side_colors = expr$metadata[, .(Contrast = baseline)])
  } else if(input$plot_type == 'Scatterplot') {
    if(input$plot_data == 'Gene Expression')
      ggplot(expr$expr %>% reshape2::melt(value.name = 'Expression', varnames = c('Gene', 'Sample')) %>%
               merge(expr$metadata, by.x = 'Sample', by.y = 'name') %>%
               mutate(Gene = as.factor(Gene), Contrast = baseline),
             aes(Sample, Expression, color = Gene)) +
      geom_point(aes(shape = Contrast), size = 2) + geom_line(aes(group = interaction(Gene, Contrast))) +
      theme_classic() + theme(axis.text.x = element_blank())
  }
}

generateResults <- function(data, experiments, conditions, options = getOption('app.all_options')) {
  outputColumns <- c('Contrast', 'Direction', 'Evidence', 'P-value', colnames(experiments)[1:(ncol(experiments) - 2)])
  
  conditions[, Evidence := paste0('<span data-toggle="popover" title="Experiments" data-html="true" data-content="',
                                  lapply(unlist(strsplit(Evidence, ',')), function(experiment) {
                                    paste0('<a target=_blank href=https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=',
                                           experiment, '>',
                                           data@experiment.meta[ee.ID == experiment, unique(ee.Name)], '</a>')
                                  }) %>% paste0(collapse = ', '), '">', paste(N, paste0('Experiment', ifelse(N > 1, 's', '')), '<i class="fas fa-question-circle" style="cursor: pointer;"></i>'), '</span>'),
             .(cf.BaseLongUri, cf.ValLongUri)]
  
  conditions[, Contrast := paste0('<b>', cf.BaseLongUri, '</b> vs. <b>', cf.ValLongUri, '</b>')]
  
  mTable <- datatable(conditions[, outputColumns, with = F] %>% as.data.frame,
                      extensions = c('FixedHeader', 'Buttons'),
                      rownames = conditions[, cf.Cat],
                      escape = -(c(which(outputColumns == 'Contrast'), which(outputColumns == 'Evidence')) + 1),
                      filter = 'top',
                      options = list(pageLength = 10,
                                     order = list(
                                       list(which(outputColumns == 'P-value'), 'asc')),
                                     language = list(lengthMenu = 'Show _MENU_ conditions per page',
                                                     processing = '',
                                                     emptyTable = 'No matching conditions found.',
                                                     infoEmpty = 'Showing 0 to 0 of 0 over condition-comparisons',
                                                     info = 'Showing _START_ to _END_ of _TOTAL_ condition-comparisons',
                                                     infoFiltered = '(filtered from over _MAX_)'),
                                     fixedHeader = T,
                                     initComplete = JS('onTableCreated'),
                                     drawCallback = JS('onTableDraw'),
                                     dom = 'lBfrtip',
                                     autoWidth = T,
                                     columnDefs = list(
                                       list(targets = 0,
                                            width = '10%', className = 'cf-cat'),
                                       list(targets = which(outputColumns == 'Contrast'),
                                            width = '25%',
                                            searchable = F, orderable = F),
                                       list(targets = which(outputColumns == 'Evidence'),
                                            className = 'dt-right',
                                            searchable = F, orderable = F),
                                       list(targets = which(outputColumns == 'P-value'),
                                            render = JS('asPval'), width = '12.5%'),
                                       list(targets = which(outputColumns == 'Direction'),
                                            render = JS('asSparkline2'), width = '1px', className = 'dt-center', searchable = F, orderable = F),
                                       list(targets = (which(outputColumns == 'P-value') + 1):length(outputColumns),
                                            render = JS('asSparkline'), width = '1px', className = 'dt-center', searchable = F, orderable = F)
                                     ),
                                     search = list(
                                       list(regex = T),
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
                  script(src = 'https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'),
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
                               pickerInput('scope', 'Ontologies', ONTOLOGIES[, unique(as.character(OntologyScope))], selected = getOption('app.ontology'), multiple = T,
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
                      mainPanel(
                        tabsetPanel(
                          tabPanel('Conditions', column(12, style = 'margin-top: 16px', dataTableOutput('results') %>% withSpinner)),
                          tabPanel('Genes', column(12, style = 'margin-top: 16px', verbatimTextOutput('test') %>% withSpinner))
                        ),
                        width = 12
                      )
                    )
                ),
                
                HTML('
<footer>
  <div class="text-center">© 2020 Copyright: <a href="mailto:jordan.sicherman@msl.ubc.ca">Jordan Sicherman</a></div>
</footer>'
                )
)
