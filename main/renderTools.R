#' Generate Results Header
#'
#' @param title The text or HTML to render
generateResultsHeader <- function(title) {
  if('html' %in% class(title))
    fluidRow(class = 'info-text', column(12, title))
  else
    fluidRow(class = 'info-text', column(12, h2(title)))
}

# TODO Implement condition picker to not overwhelm
generateGeneContribs <- function(data, options, plot_conditions = NULL) {
  mData <- data[, !c('Test Statistic', 'Ontology Steps')] %>%
    head(20) %>% # TODO Stopgap for picker
    melt(id.vars = 1:3) %>%
    .[, contrast := paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)]
  
  fig <- plot_ly(type = 'scatterpolar', fill = 'toself', mode = 'markers')
  for(group in unique(mData[, contrast])) {
    fig <- fig %>%
      add_trace(r = mData[contrast == group, value],
                theta = mData[contrast == group, variable],
                visible = 'legendonly',
                name = group)
  }
  
  renderPlotly(fig %>% layout(polar = list(radialaxis = list(visible = T, range = c(0, 1)))))
}

#' Generate Results Plot
#'
#' @param genes Genes that can be visualized
#' @param conditions Conditions that can be visualized
#' @param expr Expression information to plot
#' @param options Any additional options 
#' @param plot_genes Genes to visualize
#' @param plot_conditions Conditions to visualize
#' @param plot_type The plot type (heatmap, boxplot, etc.)
#' @param plot_data The data to plot (gene expression data, etc.)
generateResultsPlot <- function(genes, conditions, expr, options = getConfig(),
                                plot_genes, plot_conditions, plot_type, plot_data) {
  if(is.null(expr) || is.null(plot_conditions))
    return(NULL)
  
  plot_conditions <- expr$metadata[grepl(paste0(gsub('([.|()\\^{}+$*?]|\\[|\\])', '\\\\\\1', unlist(strsplit(plot_conditions, ' vs. ', T))), collapse = '|'), baseline), name]
  
  expr$expr <- expr$expr[plot_genes, plot_conditions]
  
  if(length(plot_genes) == 1)
    expr$expr <- t(expr$expr)
  
  if(plot_type != 'Boxplot') {
    thresh <- 1.5 * apply(expr$expr, 1, iqr, na.rm = T)
    quants <- rowQuantiles(expr$expr, probs = c(0.25, 0.75), na.rm = T)
    
    inliers <- suppressWarnings((expr$expr < (quants[, 1] - thresh) |
                                   expr$expr > (quants[, 2] + thresh)) %>%
                                  apply(2, max, na.rm = T)) %>% {
                                    names(.[. == 0])
                                  }
    
    expr$metadata <- expr$metadata[name %in% inliers] %>% setorder(baseline)
    expr$expr <- expr$expr[, expr$metadata$name]
    
    if(length(plot_genes) == 1)
      expr$expr <- t(expr$expr)
  }
  
  if(any(dim(expr$expr) == 0))
    return(NULL)
  
  if(plot_type == 'Heatmap') {
    if(plot_data == 'Gene Expression')
      heatmaply(expr$expr, labRow = NULL, showticklabels = c(F, T), scale = 'row',
                Rowv = NULL, dendrogram = 'none', colors = cool_warm,
                col_side_colors = expr$metadata[, .(Contrast = baseline)])
  } else {
    if(plot_data == 'Gene Expression') {
      data <- expr$expr %>% reshape2::melt(value.name = 'Expression', varnames = c('Gene', 'Sample')) %>%
        merge(expr$metadata, by.x = 'Sample', by.y = 'name') %>%
        mutate(Gene = as.factor(Gene), Contrast = baseline)
    }
    
    if(plot_type == 'Scatterplot') {
      (data %>% ggplot(aes(color = Gene, x = Sample, y = Expression)) +
         geom_point(aes(text = paste('Accession:', accession), shape = Contrast), size = 2) +
         scale_color_brewer(palette = 'Dark2') +
         geom_line(aes(group = interaction(Gene, Contrast))) +
         theme_classic() + theme(axis.text.x = element_blank())) %>%
        ggplotly %>%
        layout(yaxis = list(title = 'Expression (log<sub>2</sub> CPM)'))
    } else if(plot_type == 'Boxplot') {
      suppressWarnings(suppressMessages((data %>% ggplot(aes(fill = Contrast,
                                                             x = Gene, y = Expression)) +
                                           geom_boxplot() +
                                           scale_fill_brewer(palette = 'Dark2') +
                                           theme_classic()) %>%
                                          ggplotly(dynamicTicks = T) %>%
                                          layout(boxmode = 'group',
                                                 yaxis = list(title = 'Expression (log<sub>2</sub> CPM)'))))
    } else if(plot_type == 'Jitterplot') {
      data <- data %>% mutate(GeneID = Gene, Gene = as.numeric(Gene))
      
      suppressWarnings(suppressMessages((data %>% ggplot(aes(text = paste0('Accession: ', accession, '<br>Gene: ', GeneID),
                                                             fill = Contrast, group = interaction(Gene, Contrast),
                                                             x = Gene, y = Expression)) +
                                           geom_jitter(position = position_jitterdodge(), alpha = 0.8) +
                                           scale_fill_brewer(palette = 'Dark2') +
                                           theme_classic()) %>%
                                          ggplotly(dynamicTicks = T, tooltip = c('text', 'fill', 'y')) %>%
                                          layout(boxmode = 'group',
                                                 xaxis = list(title = 'Gene', tickmode = 'array', autotick = F, tickvals = 1:(length(unique(data$Gene))), ticktext = unique(data$GeneID)),
                                                 yaxis = list(title = 'Expression (log<sub>2</sub> CPM)'))))
    } else if(plot_type == 'Violin plot') {
      suppressWarnings(suppressMessages((data %>% ggplot(aes(fill = Contrast,
                                                             x = Gene, y = Expression)) +
                                           geom_violin(alpha = 0.9) +
                                           scale_fill_brewer(palette = 'Dark2') +
                                           theme_classic()) %>%
                                          ggplotly %>%
                                          layout(yaxis = list(title = 'Expression (log<sub>2</sub> CPM)'))))
    }
  }
}

makeNetwork <- memoise::memoise(function(data, options) {
  mTag <- getTags(options$taxa$value) %>%
    .[cf.Cat %in% data[, unique(cf.Cat)] &
        cf.BaseLongUri %in% data[, unique(cf.BaseLongUri)] &
        cf.ValLongUri %in% data[, unique(cf.ValLongUri)]] %>%
    .[, Contrast := paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri)] %>%
    .[, .(rsc.ID = as.integer(as.factor(rsc.ID)), cf.Cat, Contrast, distance)] %>%
    .[distance == 0, N := .N, Contrast] %>% .[N > 1, rsc.ID := head(.SD$rsc.ID, 1), .(cf.Cat, Contrast)] %>%
    .[, distance := mean(distance), .(rsc.ID, Contrast)] %>%
    .[, unique(.SD[, !'N']), rsc.ID] %>%
    .[, setorder(.SD, distance), rsc.ID] %>%
    .[, .(level = as.integer(distance / max(distance) * length(unique(distance)) - min(distance)), cf.Cat, Contrast), rsc.ID] %>%
    .[, parent := level - 1] %>%
    .[, Parent := factor(parent, levels = unique(parent),
                         labels = c(NA, .SD[level < max(level), .(Contrast = head(Contrast, 1)), level][, Contrast])) %>%
        as.character, rsc.ID] %>% .[, .(cf.Cat, Parent, Contrast)] %>%
    .[is.na(Parent), Parent := cf.Cat] %>% unique %>%
    merge(data[`Test Statistic` > 0, .(cf.Cat, Contrast = paste0(cf.BaseLongUri, ' vs. ', cf.ValLongUri), `Test Statistic`)],
          by = c('cf.Cat', 'Contrast'), all = T) %>% {
            rbind(data.table(Parent = NA, Contrast = 'Data', `Test Statistic` = NA, children = NA),
                  rbind(.[, .(Parent = 'Data', Contrast = unique(cf.Cat), `Test Statistic` = NA, children = NA)],
                        .[is.na(Parent), Parent := cf.Cat] %>%
                          .[, children := .N, cf.Cat] %>%
                          .[, .(Parent, Contrast, `Test Statistic`, children)]))
          }
  
  mGraph <- suppressWarnings(
    igraph::graph_from_data_frame(mTag[, .(Parent, Contrast, weight = 1 / (`Test Statistic` + abs(min(`Test Statistic`) + 1)))] %>%
                                    .[is.na(weight), weight := 1]) %>%
      simplify(edge.attr.comb = 'sum')
  )
  
  if(!igraph::is_dag(mGraph))
    mGraph <- suppressWarnings(igraph::induced_subgraph(mGraph, igraph::topo_sort(mGraph)))
  
  vertices <- mTag[, .(`Test Statistic`, Contrast)] %>% .[, head(.SD, 1), `Test Statistic`] %>%
    setorder(-`Test Statistic`, na.last = T) %>% .[, Contrast] %>% intersect(names(igraph::V(mGraph))) %>%
    head(500)
  
  paths <- shortest_paths(mGraph, from = 'NA',
                          to = vertices,
                          mode = 'out')$vpath %>% unlist %>% unique
  
  mSimpleGraph <- igraph::induced_subgraph(mGraph, paths)
  dGraph <- igraph::as_data_frame(mSimpleGraph)
  
  mTag[is.na(Parent) | (Parent %in% dGraph$from & Contrast %in% dGraph$to)] %>%
    .[Contrast %in% Parent & !(Contrast %in% vertices), `Test Statistic` := NA] %>%
    .[, `Test Statistic` := `Test Statistic` / children] %>%
    .[, .(Parent, Contrast, `Test Statistic`)]
})

generateResultsCloud <- function(data, options) {
  mCloudValues <- data[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, `Test Statistic`)] %>%
    melt(measure.vars = 1:3) %>%
    .[, .(`Test Statistic` = mean(`Test Statistic`), variable), value] %>%
    setorder(-`Test Statistic`) %>% .[`Test Statistic` > 0] %>%
    head(300)
  
  mVals <- mCloudValues[, value]
  mCounts <- mCloudValues[, `Test Statistic`]
  
  mCloud <- d3wordcloud(mVals, mCounts)
  renderD3wordcloud(mCloud)
}

generateResultsTree <- function(data, options) {
  mCircle <- makeNetwork(data, options) %>%
    as.data.frame %>% .[!is.na(.[, 1]), ] %>%
    data.tree::FromDataFrameNetwork() %>%
    circlepackeR(size = 'Test Statistic',
                 color_max = 'hsl(211, 100%, 14%)',
                 color_min = 'hsl(210, 53%, 74%)')
   
   renderCirclepackeR(mCircle)
}

#' Generate Results
#' 
#' Make a pretty results table and render it
#'
#' @param session The Shiny session storing our data
generateResults <- function(data) {
  outputColumns <- c('Contrast', 'Evidence', 'Observations', 'Ontology Steps', 'Test Statistic')
  
  mTable <- datatable(data[, outputColumns, with = F] %>% as.data.frame,
                      extensions = 'Buttons',
                      rownames = data[, cf.Cat],
                      escape = -(c(which(outputColumns == 'Contrast'), which(outputColumns == 'Evidence')) + 1),
                      filter = 'top',
                      options = list(pageLength = 10,
                                     order = list(
                                       list(which(outputColumns == 'Test Statistic'), 'desc'),
                                       list(which(outputColumns == 'Ontology Steps'), 'asc'),
                                       list(which(outputColumns == 'Observations'), 'desc')),
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
                                     deferRender = T,
                                     serverSide = T,
                                     columnDefs = list(
                                       list(targets = 0,
                                            width = '10%',
                                            searchable = T,
                                            className = 'cf-cat'),
                                       list(targets = which(outputColumns == 'Contrast'),
                                            width = '50%',
                                            searchable = T, orderable = F),
                                       list(targets = which(outputColumns == 'Evidence'),
                                            width = '10%',
                                            className = 'dt-right',
                                            searchable = F, orderable = F),
                                       list(targets = which(outputColumns == 'Test Statistic'),
                                            render = JS('asPval'),
                                            width = '10%',
                                            searchable = F),
                                       list(targets = which(outputColumns == 'Observations'),
                                            width = '5%',
                                            searchable = F),
                                       list(targets = which(outputColumns == 'Ontology Steps'),
                                            width = '5%')#,
                                       #list(targets = which(outputColumns == 'Direction'),
                                       #     width = '5%',
                                       #     render = JS('asSparkline2'), width = '1px', className = 'dt-center', searchable = F, orderable = F)#,
                                       #list(targets = (which(outputColumns == 'P-value') + 1):length(outputColumns),
                                       #     render = JS('asSparkline'), width = '1px', className = 'dt-center', searchable = F, orderable = F)
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
                                            title = 'data'#,
                                            #exportOptions = list(
                                            #  format = list(
                                            #body = JS('unformatSpark')
                                            # TODO Will eventually need to fix this so html data is fixed
                                            #customizeData = JS('function(data) { console.log(data); }')
                                            #  )
                                            #)
                                       )
                                     )
                      )
  )
  mTable$dependencies <- append(mTable$dependencies, htmlwidgets:::getDependency('sparkline'))
  
  renderDT(mTable)
}

#' Generate GO Page
#'
#' @param ora The result of an ermineR::ora
generateGOPage <- function(ora) {
  ora$results <- ora$results %>% rename(c('Pval' = 'P-value', 'CorrectedPvalue' = 'FDR'))
  mTable <- datatable(ora$results[, c('Name', 'P-value', 'FDR')],
                      rownames = ora$results$ID,
                      filter = 'top',
                      options = list(pageLength = 10,
                                     order = list(
                                       list(2, 'asc')
                                     ),
                                     language = list(lengthMenu = 'Show _MENU_ enrichments per page',
                                                     processing = '',
                                                     emptyTable = 'No matching enrichments found.',
                                                     infoEmpty = 'Showing 0 to 0 of 0 over enrichments',
                                                     info = 'Showing _START_ to _END_ of _TOTAL_ enrichments',
                                                     infoFiltered = '(filtered from over _MAX_)'),
                                     fixedHeader = T,
                                     autoWidth = T,
                                     drawCallback = JS('addMathJax'),
                                     search = list(
                                       list(regex = T)
                                     ),
                                     columnDefs = list(
                                       list(targets = 2:3,
                                            render = JS('asPval'), width = '10%')
                                     )))
  renderDT(mTable)
} # TODO May memoise

#' Generate Gene Page
#'
#' @param evidence The evidence that was fetched from Gemma (@seealso geneEvidence)
generateGenePage <- function(evidence) {
  evidence <- evidence[sapply(evidence, Negate(is.null))]
  
  evidenceCodes <- function(code) {
    switch(code,
           EXP = c('Inferred from Experiment', 'http://wiki.geneontology.org/index.php/Inferred_from_Experiment_(EXP)'),
           IDA = c('Inferred from Direct Assay', 'http://wiki.geneontology.org/index.php/Inferred_from_Direct_Assay_(IDA)'),
           IPI = c('Inferred from Physical Interaction', 'http://wiki.geneontology.org/index.php/Inferred_from_Physical_Interaction_(IPI)'),
           IMP = c('Inferred from Mutant Phenotype', 'http://wiki.geneontology.org/index.php/Inferred_from_Mutant_Phenotype_(IMP)'),
           IGI = c('Inferred from Genetic Interaction', 'http://wiki.geneontology.org/index.php/Inferred_from_Genetic_Interaction_(IGI)'),
           IEP = c('Inferred from Expression Pattern', 'http://wiki.geneontology.org/index.php/Inferred_from_Expression_Pattern_(IEP)'),
           
           HTP = c('Inferred from High Throughput Experiment', 'http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Experiment_(HTP)'),
           HDA = c('Inferred from High Throughput Direct Assay', 'http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Direct_Assay_(HDA)'),
           HMP = c('Inferred from High Throughput Mutant Phenotype', 'http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Mutant_Phenotype_(HMP)'),
           HGI = c('Inferred from High Throughput Genetic Interaction', 'http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Genetic_Interaction_(HGI)'),
           HEP = c('Inferred from High Throughput Expression Pattern', 'http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Expression_Pattern_(HEP)'),
           
           IBA = c('Inferred from Biological aspect of Ancestor', 'http://wiki.geneontology.org/index.php/Inferred_from_Biological_aspect_of_Ancestor_(IBA)'),
           IBD = c('Inferred from Biological aspect of Descendant', 'http://wiki.geneontology.org/index.php/Inferred_from_Biological_aspect_of_Descendant_(IBD)'),
           IKR = c('Inferred from Key Residues', 'http://wiki.geneontology.org/index.php/Inferred_from_Key_Residues_(IKR)'),
           IRD = c('Inferred from Rapid Divergence', 'http://wiki.geneontology.org/index.php/Inferred_from_Rapid_Divergence(IRD)'),
           
           ISS = c('Inferred from Sequence or structural Similarity', 'http://wiki.geneontology.org/index.php/Inferred_from_Sequence_or_structural_Similarity_(ISS)'),
           ISO = c('Inferred from Sequence Orthology', 'http://wiki.geneontology.org/index.php/Inferred_from_Sequence_Orthology_(ISO)'),
           ISA = c('Inferred from Sequence Alignment', 'http://wiki.geneontology.org/index.php/Inferred_from_Sequence_Alignment_(ISA)'),
           ISM = c('Inferred from Sequence Model', 'http://wiki.geneontology.org/index.php/Inferred_from_Sequence_Model_(ISM)'),
           IGC = c('Inferred from Genomic Context', 'http://wiki.geneontology.org/index.php/Inferred_from_Genomic_Context_(IGC)'),
           RCA = c('Inferred from Reviewed Computational Analysis', 'http://wiki.geneontology.org/index.php/Inferred_from_Reviewed_Computational_Analysis_(RCA)'),
           
           TAS = c('Traceable Author Statement', 'http://wiki.geneontology.org/index.php/Traceable_Author_Statement_(TAS)'),
           NAS = c('Non-traceable Author Statement', 'http://wiki.geneontology.org/index.php/Non-traceable_Author_Statement_(NAS)'),
           
           IC = c('Inferred by Curator', 'http://wiki.geneontology.org/index.php/Inferred_by_Curator_(IC)'),
           ND = c('No biological Data available', 'http://wiki.geneontology.org/index.php/No_biological_Data_available_(ND)_evidence_code'),
           
           IEA = c('Inferred from Electronic Annotation', 'http://wiki.geneontology.org/index.php/Inferred_from_Electronic_Annotation_(IEA)'),
           
           IED = c('Inferred from Experimental Data', 'Inferred from Experimental Data'),
           IAGP = c('Inferred from Association of Genotype and Phenotype', 'Inferred from Association of Genotype and Phenotype'),
           IPM = c('Inferred from Phenotype Manipulation', 'Inferred from Phenotype Manipulation'),
           QTM = c('Quantitative Trait Measurement', 'Quantitative Trait Measurement'),
           IIA = c('Inferred from Imported Annotation', 'Inferred from Imported Annotation')
    )
  }
  
  if(length(evidence) <= 8)
    colors <- brewer.pal(max(3, length(evidence)), 'Dark2')
  else
    colors <- colorRampPalette(brewer.pal(8, 'Dark2'))(length(evidence))
  
  content <- lapply(1:length(evidence), function(indx) {
    gene <- evidence[[indx]]
    
    if(is.null(gene)) return(NULL)
    else
      shinypanels::panel(title = gene[[1]]$symbol, collapsed = T, width = 500, color = colors[indx],
                         lapply(gene, function(evidence) {
                           cat <- ONTOLOGIES.DEFS[Node_Long == evidence$cf.CatLongUri, data.table::first(Definition)]
                           if(!is.character(cat))
                             cat <- evidence$cf.CatLongUri
                           
                           val <- ONTOLOGIES.DEFS[Node_Long == evidence$cf.ValLongUri, data.table::first(Definition)]
                           if(!is.character(val))
                             val <- evidence$cf.ValLongUri
                           
                           div(p(tags$b('Category: '), cat),
                               p(tags$b('Value: '), val),
                               p(tags$b('Relationship: '), evidence$relationship),
                               p(tags$b('Evidence Code: '), paste0(evidence$evidence.Code, ' '),
                                 a(target = '_blank', href = evidenceCodes(evidence$evidence.Code)[2], tags$i(class = 'fas fa-question-circle'))),
                               p(tags$b(paste0(ifelse(evidence$score.Name == '', 'Score', evidence$score.Name), ': ')),
                                 ifelse(is.null(evidence$score.Value), evidence$score.Strength, evidence$score.Value)),
                               p(tags$b('Citation: '), a(target = '_blank', href = evidence$citation.Url, evidence$citation.Name))
                           )
                         }))
  })
  
  renderUI(div(class = 'layout-container',
               div(class = 'layout-panels',
                   div(class = 'app-container', content)
               )
  ))
}
