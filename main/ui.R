ui <- fluidPage(style = 'height: 100%;',
                useShinyjs(),
                disconnectMessage(text = 'Your session has been disconnected.',
                                  refresh = '',
                                  background = '#FA1919E6',
                                  colour = '#FFFFFF',
                                  overlayColour = '#999999',
                                  overlayOpacity = 0.7,
                                  width = 'full',
                                  top = 'center',
                                  size = 60,
                                  css = 'padding: 15px !important; box-shadow: none !important'),
                
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
                    fluidRow(id = 'header', HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 268.12 59.45" xmlns:v="https://vecta.io/nano"><path d="M79.46 37.27v-5.91h15.26v13.97c-1.48 1.44-3.63 2.7-6.45 3.79s-5.67 1.64-8.55 1.64c-3.67 0-6.87-.77-9.59-2.31-2.73-1.54-4.78-3.74-6.15-6.6s-2.06-5.98-2.06-9.34c0-3.65.77-6.9 2.3-9.74s3.77-5.02 6.72-6.53c2.25-1.16 5.05-1.75 8.4-1.75 4.35 0 7.76.91 10.2 2.74 2.45 1.83 4.02 4.35 4.73 7.57l-7.03 1.32c-.49-1.72-1.42-3.08-2.79-4.08-1.36-1-3.07-1.5-5.11-1.5-3.09 0-5.55.98-7.38 2.94s-2.74 4.87-2.74 8.73c0 4.16.92 7.29 2.78 9.37 1.85 2.08 4.27 3.12 7.27 3.12 1.48 0 2.97-.29 4.46-.87s2.77-1.29 3.84-2.12v-4.45h-8.11zm34.22 4.81l6.7 1.12c-.86 2.46-2.22 4.33-4.08 5.61s-4.18 1.93-6.97 1.93c-4.42 0-7.69-1.44-9.81-4.33-1.68-2.32-2.52-5.24-2.52-8.76 0-4.21 1.1-7.51 3.3-9.89s4.98-3.58 8.35-3.58c3.78 0 6.76 1.25 8.95 3.74 2.19 2.5 3.23 6.32 3.13 11.47h-16.84c.05 1.99.59 3.54 1.63 4.65s2.33 1.66 3.88 1.66c1.05 0 1.94-.29 2.66-.86s1.25-1.49 1.62-2.76zm.38-6.8c-.05-1.95-.55-3.42-1.51-4.44-.96-1.01-2.12-1.52-3.49-1.52a4.7 4.7 0 0 0-3.64 1.6c-.96 1.06-1.43 2.52-1.41 4.35h10.05zm9.09-10.53h6.2v3.47c2.22-2.7 4.86-4.04 7.92-4.04 1.63 0 3.04.33 4.23 1 1.2.67 2.18 1.68 2.94 3.04 1.12-1.36 2.32-2.37 3.61-3.04s2.67-1 4.14-1c1.87 0 3.45.38 4.74 1.14s2.26 1.87 2.9 3.34c.46 1.08.69 2.84.69 5.26v16.25h-6.72V35.64c0-2.52-.23-4.15-.69-4.88-.62-.96-1.58-1.44-2.87-1.44-.94 0-1.83.29-2.66.86s-1.43 1.42-1.79 2.52c-.37 1.11-.55 2.86-.55 5.25v12.2h-6.72V36.24c0-2.47-.12-4.07-.36-4.79s-.61-1.25-1.11-1.6-1.18-.53-2.05-.53c-1.04 0-1.97.28-2.8.84s-1.42 1.36-1.78 2.42c-.36 1.05-.54 2.8-.54 5.24v12.35h-6.72V24.75zm41.01 0h6.2v3.47c2.22-2.7 4.86-4.04 7.92-4.04 1.63 0 3.04.33 4.23 1 1.2.67 2.18 1.68 2.94 3.04 1.12-1.36 2.32-2.37 3.61-3.04s2.67-1 4.14-1c1.87 0 3.45.38 4.74 1.14s2.26 1.87 2.89 3.34c.46 1.08.69 2.84.69 5.26v16.25h-6.72V35.64c0-2.52-.23-4.15-.69-4.88-.62-.96-1.58-1.44-2.87-1.44-.94 0-1.83.29-2.66.86s-1.43 1.42-1.79 2.52c-.37 1.11-.55 2.86-.55 5.25v12.2h-6.72V36.24c0-2.47-.12-4.07-.36-4.79s-.61-1.25-1.11-1.6-1.19-.53-2.05-.53c-1.04 0-1.97.28-2.8.84s-1.42 1.36-1.78 2.42c-.36 1.05-.54 2.8-.54 5.24v12.35h-6.72V24.75zm46.41 7.76l-6.1-1.1c.69-2.46 1.87-4.27 3.54-5.46s4.16-1.77 7.46-1.77c3 0 5.23.36 6.7 1.06 1.47.71 2.5 1.61 3.1 2.7s.9 3.1.9 6.02l-.07 7.85c0 2.23.11 3.88.32 4.94s.62 2.2 1.21 3.41h-6.65c-.18-.45-.39-1.11-.65-1.99l-.24-.79c-1.15 1.12-2.38 1.95-3.69 2.51s-2.7.84-4.19.84c-2.62 0-4.68-.71-6.18-2.13-1.51-1.42-2.26-3.21-2.26-5.38 0-1.44.34-2.72 1.03-3.84a6.74 6.74 0 0 1 2.88-2.58c1.24-.6 3.02-1.12 5.35-1.57 3.14-.59 5.32-1.14 6.53-1.65v-.67c0-1.29-.32-2.21-.96-2.76s-1.84-.83-3.61-.83c-1.2 0-2.13.24-2.8.71-.66.47-1.2 1.29-1.62 2.48zm9 5.45c-.86.29-2.22.63-4.09 1.03s-3.09.79-3.66 1.17c-.88.62-1.32 1.41-1.32 2.37 0 .94.35 1.75 1.05 2.44s1.59 1.03 2.68 1.03c1.21 0 2.37-.4 3.47-1.2.81-.61 1.35-1.35 1.6-2.22.18-.57.26-1.67.26-3.28v-1.34z" fill="#fff"/><path d="M29.72 6.87L9.94 18.3v22.85l19.78 11.42 19.79-11.42V18.3z" fill="#00718f"/><path d="M26.6 13.32l-9.52 5.5c-1.93 1.12-3.13 3.18-3.13 5.41v10.99c0 2.23 1.19 4.3 3.13 5.41l9.52 5.5c1.93 1.12 4.32 1.12 6.25 0l9.52-5.5c1.93-1.12 3.13-3.18 3.13-5.41V24.23c0-2.23-1.19-4.3-3.13-5.41l-9.52-5.5c-1.94-1.12-4.32-1.12-6.25 0z" fill="#f8f6db"/><circle cx="29.72" cy="29.72" r="5.16" fill="#fcb534"/><g fill="#fddb5c"><circle cx="29.72" cy="17.71" r="5.16"/><circle cx="29.72" cy="41.74" r="5.16"/></g><circle cx="40.13" cy="23.72" r="5.16" fill="#b63926"/><circle cx="19.32" cy="35.73" r="5.16" fill="#f58020"/><circle cx="19.32" cy="23.72" r="5.16" fill="#fddb5c"/><circle cx="40.13" cy="35.73" r="5.16" fill="#fcb534"/><path d="M29.72 0l-2.69 10.61h5.39zM14.86 3.98l2.97 10.54 4.67-2.7zM3.98 14.86l7.84 7.64 2.7-4.67zM0 29.72l10.61 2.7v-5.39zm3.98 14.86l10.54-2.96-2.7-4.67zm10.88 10.88l7.64-7.83-4.67-2.7zm14.86 3.99l2.7-10.61h-5.39zm14.86-3.99l-2.96-10.53-4.67 2.7zm10.88-10.88l-7.83-7.63-2.7 4.67zm3.99-14.86l-10.61-2.69v5.39zm-3.99-14.86l-10.53 2.97 2.7 4.67zM44.58 3.98l-7.63 7.84 4.67 2.7z" fill="#00718f"/><path d="M234.19 28.35h-2.07c-1.1-1.65-1.93-3.37-2.5-5.16s-.86-3.51-.86-5.18c0-2.07.35-4.03 1.06-5.88a21.33 21.33 0 0 1 2.34-4.44h2.06c-.98 2.17-1.66 4.02-2.03 5.54s-.55 3.14-.55 4.84a20.54 20.54 0 0 0 .33 3.61c.22 1.23.52 2.4.9 3.51.25.74.69 1.79 1.32 3.16zm1.03-20.38h5.81c1.31 0 2.31.1 3 .3a5.3 5.3 0 0 1 2.37 1.45c.66.69 1.16 1.55 1.5 2.55.34 1.01.52 2.25.52 3.72 0 1.3-.16 2.41-.48 3.35-.39 1.15-.96 2.07-1.69 2.78-.55.54-1.3.96-2.23 1.26-.7.22-1.64.33-2.81.33h-5.98V7.97zm3.18 2.66v10.43h2.37c.89 0 1.53-.05 1.92-.15.52-.13.94-.35 1.28-.66s.62-.81.83-1.52.32-1.67.32-2.88c0-1.22-.11-2.15-.32-2.8s-.52-1.16-.9-1.53-.88-.61-1.47-.74c-.44-.1-1.31-.15-2.61-.15h-1.42zm11.68 13.09V7.97h11.68v2.66h-8.5v3.49h7.91v2.65h-7.91v4.29h8.8v2.65h-11.98zm12.6 4.63l1.26-2.93c.25-.68.47-1.46.68-2.35a23.92 23.92 0 0 0 .46-2.53c.1-.8.15-1.62.15-2.45 0-1.7-.18-3.32-.55-4.84s-1.04-3.37-2.02-5.54h2.05c1.08 1.54 1.92 3.17 2.52 4.9s.9 3.48.9 5.25c0 1.5-.24 3.1-.71 4.81-.54 1.92-1.42 3.81-2.65 5.68h-2.09z" fill="#fff"/></svg>')),
                    
                    # The gene entry bar
                    fluidRow(wellPanel(style = 'margin-bottom: 0;',
                      # Search parameters
                      fluidRow(
                        # Gene entry
                        column(6,
                               fluidRow(span(selectInput('genes', 'Input gene(s) of interest', list(`Enter query...` = ''), multiple = T), `data-toggle` = 'tooltip', title = 'May be NCBI IDs, Ensembl IDs, official symbols or GO groups')),
                               fluidRow(helpText(HTML('Examples: <a genes=\'["RPS4Y1","EIF1AY","DDX3Y","KDM5D","XIST"]\'>RPS4Y1, EIF1AY, DDX3Y, KDM5D, XIST</a>&nbsp;&nbsp;|&nbsp;&nbsp;<a genes=\'["ENSG00000131095","ENSG00000110436"]\'>ENSG00000131095, ENSG00000110436</a>&nbsp;&nbsp;|&nbsp;&nbsp;<a style="color: #002145" genes=\'"random()"\'>I\'m feeling lucky</a>')))),
                        column(2, uiOutput('genes.csv.ui')),
                        
                        # Signature entry
                        column(2, span(selectInput('sig', 'DE signature', list(`(optional)` = ''), multiple = T), `data-toggle` = 'tooltip', title = 'If providing a signature, you should use M-VSM or Correlation as the search algorithm')),
                        
                        # Taxa entry
                        column(2,
                               span(pickerInput('taxa', 'Taxon', getConfig('taxa')$choices, getConfig('taxa')$value, multiple = T), `data-toggle` = 'tooltip', title = 'Taxa to include. Homologs are automatically detected'),
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
                      column(12, htmlOutput('results_suggestions')),
                      column(12, htmlOutput('results_header')),
                      mainPanel(
                        tabsetPanel(id = 'tabs',
                                    tabPanel('Conditions',
                                             column(12, style = 'margin-top: 16px', dataTableOutput('results') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase')))))),
                                    tabPanel('Gene Contributions', column(12, style = 'margin-top: 16px', plotlyOutput('results_contribs', height = '30vw') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase')))))),
                                    tabPanel('Gene Info', column(12, style = 'margin-top: 16px', htmlOutput('results_genes') %>% withSpinner(custom.class = 'DNA_cont', custom.html = div(lapply(1:10, function(x) div(class = 'nucleobase')))))),
                        ),
                        width = 12
                      )
                    )
                ),
                
                HTML('
<footer>
  <div class="text-center">© 2020-2021 <a href="mailto:jordan.sicherman@msl.ubc.ca">Jordan Sicherman</a></div>
</footer>'
                )
)
