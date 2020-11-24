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
                  
                  link(rel = 'icon', href = 'https://gemma.msl.ubc.ca/images/favicon.ico'),
                  useShinyPanels()
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
                               fluidRow(helpText(HTML('Examples: <a genes=\'["RPS4Y1","XIST","KDM5D"]\'>RPS4Y1, XIST, KDM5D</a>, <a genes=\'"ENSG00000121410"\'>ENSG00000121410</a>, <a style="color: #002145" genes=\'"random()"\'>I\'m feeling lucky</a>')))),
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
                                                             numericInput('distance', 'Ontology expansion limit', value = getOption('app.distance_cutoff'), step = 0.25, min = 0, max = 10),
                                                             numericInput('min.tags', 'Minimum augmented tag count', value = getOption('app.min.tags'), step = 1, min = 1, max = 100000),
                                                             numericInput('max.rows', 'Maximum contrasts to compute', value = getOption('app.max.rows'), step = 10, min = 10, max = 1000))),
                                         column(6, wellPanel(`well-name` = 'Scoring',
                                                             selectInput('method', 'Scoring function', getOption('app.all_search_methods')),
                                                             materialSwitch('mfx', 'Include multifunctionality', value = getOption('app.mfx'), right = T),
                                                             materialSwitch('geeq', 'Include experiment quality (GEEQ)', value = getOption('app.geeq'), right = T),
                                                             numericInput('pv', 'Significance threshold', value = getOption('app.pv'), step = 0.01, min = 0, max = 1),
                                                             sliderInput('fc', 'FC threshold', value = c(getOption('app.fc_lower'), getOption('app.fc_upper')), step = 0.1, min = 0, max = 100, ticks = F)))
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
                        tabsetPanel(id = 'tabs',
                          tabPanel('Conditions', column(12, style = 'margin-top: 16px', dataTableOutput('results') %>% withSpinner)),
                          tabPanel('Gene Info', column(12, style = 'margin-top: 16px', htmlOutput('results_genes') %>% withSpinner)),
                          tabPanel('GO Enrichment', column(12, style = 'margin-top: 16px', dataTableOutput('results_go') %>% withSpinner))
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
