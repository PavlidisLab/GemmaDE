ui <- fluidPage(style = 'height: 100%;',
                useShinyjs(),
                
                # Add to the HTML head
                withTags(head(
                  script(src = 'js/initialize.js'),
                  script(src = 'https://kit.fontawesome.com/33dcd9d8f9.js', crossorigin = 'anonymous'),
                  script(src = 'https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'),
                  link(rel = 'stylesheet', type = 'text/css', href = 'css/style.css'),
                  style(sass(sass_file('www/css/dna_loader.scss'), options = sass_options(output_style = "compressed"))),
                  
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
                    fluidRow(wellPanel(style = 'margin-bottom: 0;',
                      # Search parameters
                      fluidRow(
                        # Gene entry
                        column(6,
                               fluidRow(selectInput('genes', 'Input gene(s) of interest', list(`Enter query...` = ''), multiple = T)),
                               fluidRow(helpText(HTML('Examples: <a genes=\'["RPS4Y1","EIF1AY","DDX3Y","KDM5D","XIST"]\'>RPS4Y1, EIF1AY, DDX3Y, KDM5D, XIST</a>, <a genes=\'["ENSG00000131095","ENSG00000110436"]\'>ENSG00000131095, ENSG00000110436</a>, <a style="color: #002145" genes=\'"random()"\'>I\'m feeling lucky</a>')))),
                        column(2, uiOutput('genes.csv.ui')),
                        
                        # Signature entry
                        column(2, textInput('sig', 'DE signature', placeholder = '(optional)')),
                        
                        # Taxa entry
                        column(2,
                               pickerInput('taxa', 'Taxon', getConfig('taxa')$choices, getConfig('taxa')$value, multiple = T),
                               helpText(style = 'float: right;', HTML('<a data-toggle="collapse" data-target="#options">More options...</a>')))
                      ),
                      
                      # More options
                      wellPanel(class = 'collapse', id = 'options',
                                fluidRow(style = 'display: flex; flex-direction: row;',
                                         column(6, wellPanel(`well-name` = 'Filtering',
                                                             lapply(getConfig(category = 'Filtering'), as.input))),
                                         column(6, wellPanel(`well-name` = 'Scoring',
                                                             lapply(getConfig(category = 'Scoring'), as.input)))
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
                    
                    # Hidden, but hijacked for downloading table data
                    downloadButton('dataDownload', '', style = 'visibility: hidden; height: 0'),
                    
                    fluidRow(
                      column(12, htmlOutput('results_header')),
                      mainPanel(
                        tabsetPanel(id = 'tabs',
                                    tabPanel('Conditions', column(12, style = 'margin-top: 16px', dataTableOutput('results') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase')))))),
                                    tabPanel('Gene Contributions', column(12, style = 'margin-top: 16px', plotlyOutput('results_contribs', height = '50vw'))),
                                    #tabPanel('Hierarchical View (beta)', column(12, style = 'margin-top: 16px', circlepackeROutput('results_tree', height = '90vh') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase')))))),
                                    tabPanel('Word Cloud', column(12, style = 'margin-top: 16px', d3wordcloudOutput('results_cloud', height = '90vh') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase')))))),
                                    tabPanel('Gene Info', column(12, style = 'margin-top: 16px', htmlOutput('results_genes') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase')))))),
                                    tabPanel('GO Enrichment', column(12, style = 'margin-top: 16px', dataTableOutput('results_go') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase'))))))
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
