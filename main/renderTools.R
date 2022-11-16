generatePlotContainer <- function() {
  fluidRow(
    style = "height: 85vh;",
    column(10, plotlyOutput("plot") %>% withSpinner(custom.class = "DNA_cont", custom.html = div(lapply(1:10, function(x) div(class = "nucleobase"))))),
    column(2,
      style = "padding-top: 15px; padding-right: 30px;",
      selectInput("plot_type", "Type", list("Boxplot", "Violin plot", "Heatmap", "Scatterplot", "Jitterplot")),
      selectInput("plot_data", "Data", list("Gene Expression")),
      pickerInput("plot_genes", "Genes", list("Loading..."), multiple = T),
      pickerInput("plot_conditions", "Contrasts", list("Loading..."), multiple = T),
      pickerInput("plot_taxa", "Taxa", list("Loading..."), multiple = T)
    )
  )
}

#' Generate Results Header
#'
#' @param title The text or HTML to render
generateResultsHeader <- function(title) {
  if ("html" %in% class(title)) {
    fluidRow(class = "info-text", column(12, title))
  } else {
    fluidRow(class = "info-text", column(12, h3(title)))
  }
}

# TODO Implement condition picker to not overwhelm
generateGeneContribs <- function(data, options, plot_conditions = NULL) {
  mData <- data %>%
    setorder(-`Test Statistic`) %>%
    .[, !c(
      "Test Statistic", "Effect Size", "Ontology Steps", "N", "Evidence",
      "cf.Cat", "cf.BaseLongUri", "cf.ValLongUri"
    )] %>%
    head(20) %>% # TODO This is a stopgap for picker
    melt(id.vars = "Condition Comparison") %>%
    .[!is.finite(value), value := 0] # %>%
  # .[, value := log10(1 + value / sum(value)), `Condition Comparison`] %>%
  # .[, value := (value - min(value)) / (max(value) - min(value)), `Condition Comparison`]
  
  # Change plot height depending of # genes
  if (length(unique(mData$variable)) < 15)
    height <- NULL
  else
    height <- 40 * length(unique(mData$variable))
    
  fig <- plot_ly(mData, height = height)
  
  for(condition in unique(mData$`Condition Comparison`)) {
    fig <- fig %>%
      add_trace(
        type = "bar",
        x = mData[`Condition Comparison` == condition, value],
        y = mData[`Condition Comparison` == condition, variable],
        name = condition,
        hoverinfo = "skip",
        visible = "legendonly"
      )
  }
  
  fig %>%
    style(visible = TRUE, traces = c(1)) %>% # Show first condition comparison
    layout(
      yaxis = list(
        type = "category",
        categoryorder = "total ascending"),
      legend = list(
        title=list(text='<b> Condition Comparison </b>'))
      ) %>%
    config(displaylogo = FALSE,
           toImageButtonOptions = list(format = "svg"),
           displayModeBar = FALSE,
           modeBarButtonsToRemove = c("toggleSpikelines", "hoverCompareCartesian")) %>%
    renderPlotly() 
}

#' Generate Results Plot
#'
#' @param genes Genes that can be visualized
#' @param conditions Conditions that can be visualized
#' @param expr Expression information to plot
#' @param options Any additional options
#' @param plot_taxa Taxa to visualize
#' @param plot_genes Genes to visualize
#' @param plot_conditions Conditions to visualize
#' @param plot_type The plot type (heatmap, boxplot, etc.)
#' @param plot_data The data to plot (gene expression data, etc.)
generateResultsPlot <- function(genes, conditions, expr, options = getConfig(),
                                plot_taxa, plot_genes, plot_conditions, plot_type, plot_data) {
  # TODO Update this for plot_taxa and new genes (if taxa is any)
  if (is.null(expr) || is.null(plot_conditions)) {
    return(NULL)
  }

  # I hate this a lot. It turns plot_conditions from something like "A vs. B" to the names of the samples with
  # annotations either or B
  plot_conditions <- expr$metadata[grepl(paste0(gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", unlist(strsplit(plot_conditions, " vs. ", T))), collapse = "|"), baseline), name]

  expr$expr <- expr$expr[plot_genes, plot_conditions, drop = F]

  if (plot_type != "Boxplot") {
    thresh <- 1.5 * apply(expr$expr, 1, iqr, na.rm = T)
    quants <- rowQuantiles(expr$expr, probs = c(0.25, 0.75), na.rm = T)

    inliers <- suppressWarnings((expr$expr < (quants[, 1] - thresh) |
      expr$expr > (quants[, 2] + thresh)) %>%
      apply(2, max, na.rm = T)) %>%
      {
        names(.[. == 0])
      }

    expr$metadata <- expr$metadata[name %in% inliers] %>% setorder(baseline)
    expr$expr <- expr$expr[, expr$metadata$name, drop = F]
  }

  if (any(dim(expr$expr) == 0)) {
    return(NULL)
  }

  if (plot_type == "Heatmap") {
    if (plot_data == "Gene Expression") {
      if (length(unique(expr$metadata$ee.Name)) > 1) {
        for (i in unique(expr$metadata$ee.Name)) {
          expr$expr[, expr$metadata[ee.Name == i, name]] <- t(scale(t(expr$expr[, expr$metadata[ee.Name == i, name]])))
        }

        ret <- heatmaply(expr$expr,
          labRow = NULL, showticklabels = c(F, T), scale = "none",
          colors = PuOr, fontsize_row = 15,
          Rowv = NULL, dendrogram = "none",
          col_side_colors = expr$metadata[, .(
            `Condition Comparison` = baseline,
            Experiment = ee.Name
          )]
        )
      } else {
        ret <- heatmaply(expr$expr,
          labRow = NULL, showticklabels = c(F, T), scale = "row",
          colors = PuOr, fontsize_row = 15,
          Rowv = NULL, dendrogram = "none",
          col_side_colors = expr$metadata[, .(`Condition Comparison` = baseline)]
        )
      }
    }
  } else {
    if (plot_data == "Gene Expression") {
      data <- expr$expr %>%
        reshape2::melt(value.name = "Expression", varnames = c("Gene", "Sample")) %>%
        merge(expr$metadata, by.x = "Sample", by.y = "name") %>%
        mutate(Gene = as.factor(Gene)) %>%
        dplyr::rename(`Condition Comparison` = baseline)
    }

    if (plot_type == "Scatterplot") {
      ret <- (data %>% ggplot(aes(color = Gene, x = Sample, y = Expression)) +
        geom_point(aes(text = paste("Accession:", accession), shape = `Condition Comparison`), size = 2) +
        scale_color_brewer(palette = "Dark2") +
        geom_line(aes(group = interaction(Gene, `Condition Comparison`))) +
        theme_classic() +
        theme(axis.text.x = element_blank())) %>%
        ggplotly() %>%
        layout(yaxis = list(title = "Expression (log<sub>2</sub> CPM)"))
    } else if (plot_type == "Boxplot") {
      ret <- suppressWarnings(suppressMessages((data %>% ggplot(aes(
        fill = `Condition Comparison`,
        x = Gene, y = Expression
      )) +
        geom_boxplot() +
        scale_fill_brewer(palette = "Dark2") +
        theme_classic()) %>%
        ggplotly(dynamicTicks = T) %>%
        layout(
          boxmode = "group",
          xaxis = list(title = "Gene", tickmode = "array", autotick = F, tickangle = -45, tickvals = 1:(length(unique(data$Gene))), ticktext = unique(data$GeneID)),
          yaxis = list(title = "Expression (log<sub>2</sub> CPM)")
        )))
    } else if (plot_type == "Jitterplot") {
      data <- data %>% mutate(GeneID = Gene, Gene = as.numeric(Gene))

      ret <- suppressWarnings(suppressMessages((data %>% ggplot(aes(
        text = paste0("Accession: ", accession, "<br>Gene: ", GeneID),
        fill = `Condition Comparison`, group = interaction(Gene, `Condition Comparison`),
        x = Gene, y = Expression
      )) +
        geom_jitter(position = position_jitterdodge(), alpha = 0.8) +
        scale_fill_brewer(palette = "Dark2") +
        theme_classic()) %>%
        ggplotly(dynamicTicks = T, tooltip = c("text", "fill", "y")) %>%
        layout(
          boxmode = "group",
          xaxis = list(title = "Gene", tickmode = "array", autotick = F, tickangle = -45, tickvals = 1:(length(unique(data$Gene))), ticktext = unique(data$GeneID)),
          yaxis = list(title = "Expression (log<sub>2</sub> CPM)")
        )))
    } else if (plot_type == "Violin plot") {
      ret <- suppressWarnings(suppressMessages((data %>% ggplot(aes(
        fill = `Condition Comparison`,
        x = Gene, y = Expression
      )) +
        geom_violin(alpha = 0.9) +
        scale_fill_brewer(palette = "Dark2") +
        theme_classic()) %>%
        ggplotly() %>%
        layout(
          xaxis = list(title = "Gene", tickmode = "array", autotick = F, tickangle = -45, tickvals = 1:(length(unique(data$Gene))), ticktext = unique(data$GeneID)),
          yaxis = list(title = "Expression (log<sub>2</sub> CPM)")
        )))
    }
  }

  ret %>% config(
    displaylogo = F,
    toImageButtonOptions = list(format = "svg"),
    modeBarButtonsToRemove = c("toggleSpikelines", "hoverCompareCartesian")
  )
}

#' Generate Results
#'
#' Make a pretty results table and render it
#'
#' @param data Our data
generateResults <- function(data) {
  outputColumns <- c("Condition Comparison", "Evidence", "Ontology Steps", "Effect Size", "Test Statistic")

  mTable <- datatable(
    data[, outputColumns, with = F] %>% 
      .[, c("Effect Size", "Test Statistic") := list(round(`Effect Size`, 2), round(`Test Statistic`, 2))] %>% 
      as.data.frame(),
    extensions = "Buttons",
    selection = "none",
    rownames = data[, as.character(cf.Cat)],
    colnames = c(
      "Condition Comparison" = "Condition Comparison",
      "<span title ='Shows list of implicated experiments.' data-toggle='tooltip'>Evidence</span>" = "Evidence",
      "<span title ='Indicates a mean of distances between the baseline and contrasting factors in ontology' data-toggle='tooltip'>Ontology Steps</span>" = "Ontology Steps",
      "<span title ='Indicates how strongly differentally-expressed the query gene list is' data-toggle='tooltip'>Effect Size</span>" = "Effect Size",
      "<span title ='Indicates how much the query genes contribute in a statistically significant way. Ranges from zero to one.' data-toggle='tooltip'>Test Statistic</span>" = "Test Statistic"
    ),
    callback = JS(
      "var a = document.createElement('a');",
      "$(a).addClass('dt-button');",
      "$(a).click(function() { window.open($('#dataDownload').attr('href')); });",
      "$(a).text('Download CSV');",
      "$('div.dwnld').append(a);"
    ),
    escape = F,
    filter = "top",
    options = list(
      bFilter = 0,
      pageLength = 10,
      order = list(
        list(which(outputColumns == "Test Statistic"), "desc"),
        list(which(outputColumns == "Ontology Steps"), "asc")
      ),
      language = list(
        lengthMenu = "Show _MENU_ condition comparisons per page",
        processing = "",
        emptyTable = "No matching condition comparisons found.",
        infoEmpty = "Showing 0 to 0 of 0 over condition comparisons",
        info = "Showing _START_ to _END_ of _TOTAL_ condition comparisons",
        infoFiltered = "(filtered from over _MAX_)"
      ),
      fixedHeader = T,
      initComplete = JS("onTableCreated"),
      drawCallback = JS("onTableDraw"),
      dom = 'lB<"dwnld">frtip',
      autoWidth = T,
      deferRender = T,
      serverSide = T,
      columnDefs = list(
        list(
          targets = 0,
          width = "10%",
          className = "cf-cat"
        ),
        list(
          targets = which(outputColumns == "Condition Comparison"),
          width = "46%",
          searchable = T, orderable = F
        ),
        list(
          targets = which(outputColumns == "Evidence"),
          width = "10%",
          className = "dt-right",
          searchable = F, orderable = F
        ),
        list(
          targets = which(outputColumns == "Test Statistic"),
          # render = JS('asPval'),
          width = "12%"
        ),
        list(
          targets = which(outputColumns == "Effect Size"),
          # render = JS('asPval'),
          width = "12%"
        ),
        list(
          targets = which(outputColumns == "Ontology Steps"),
          width = "10%"
        )
      ),
      search = list(
        list(regex = T)
      ),
      buttons = list(
       #  list(
       #    extend = "collection",
       #    text = "Visualize Expression",
       #    action = JS("plotData")
       # )
      )
    )
  )

  renderDT(mTable)
}

#' Generate Gene Page
#'
#' @param evidence The evidence that was fetched from Gemma (@seealso geneEvidence)
generateGenePage <- function(evidence) {
  evidence <- evidence[sapply(evidence, Negate(is.null))]

  evidenceCodes <- function(code) {
    switch(code,
      EXP = c("Inferred from Experiment", "http://wiki.geneontology.org/index.php/Inferred_from_Experiment_(EXP)"),
      IDA = c("Inferred from Direct Assay", "http://wiki.geneontology.org/index.php/Inferred_from_Direct_Assay_(IDA)"),
      IPI = c("Inferred from Physical Interaction", "http://wiki.geneontology.org/index.php/Inferred_from_Physical_Interaction_(IPI)"),
      IMP = c("Inferred from Mutant Phenotype", "http://wiki.geneontology.org/index.php/Inferred_from_Mutant_Phenotype_(IMP)"),
      IGI = c("Inferred from Genetic Interaction", "http://wiki.geneontology.org/index.php/Inferred_from_Genetic_Interaction_(IGI)"),
      IEP = c("Inferred from Expression Pattern", "http://wiki.geneontology.org/index.php/Inferred_from_Expression_Pattern_(IEP)"),
      HTP = c("Inferred from High Throughput Experiment", "http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Experiment_(HTP)"),
      HDA = c("Inferred from High Throughput Direct Assay", "http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Direct_Assay_(HDA)"),
      HMP = c("Inferred from High Throughput Mutant Phenotype", "http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Mutant_Phenotype_(HMP)"),
      HGI = c("Inferred from High Throughput Genetic Interaction", "http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Genetic_Interaction_(HGI)"),
      HEP = c("Inferred from High Throughput Expression Pattern", "http://wiki.geneontology.org/index.php/Inferred_from_High_Throughput_Expression_Pattern_(HEP)"),
      IBA = c("Inferred from Biological aspect of Ancestor", "http://wiki.geneontology.org/index.php/Inferred_from_Biological_aspect_of_Ancestor_(IBA)"),
      IBD = c("Inferred from Biological aspect of Descendant", "http://wiki.geneontology.org/index.php/Inferred_from_Biological_aspect_of_Descendant_(IBD)"),
      IKR = c("Inferred from Key Residues", "http://wiki.geneontology.org/index.php/Inferred_from_Key_Residues_(IKR)"),
      IRD = c("Inferred from Rapid Divergence", "http://wiki.geneontology.org/index.php/Inferred_from_Rapid_Divergence(IRD)"),
      ISS = c("Inferred from Sequence or structural Similarity", "http://wiki.geneontology.org/index.php/Inferred_from_Sequence_or_structural_Similarity_(ISS)"),
      ISO = c("Inferred from Sequence Orthology", "http://wiki.geneontology.org/index.php/Inferred_from_Sequence_Orthology_(ISO)"),
      ISA = c("Inferred from Sequence Alignment", "http://wiki.geneontology.org/index.php/Inferred_from_Sequence_Alignment_(ISA)"),
      ISM = c("Inferred from Sequence Model", "http://wiki.geneontology.org/index.php/Inferred_from_Sequence_Model_(ISM)"),
      IGC = c("Inferred from Genomic Context", "http://wiki.geneontology.org/index.php/Inferred_from_Genomic_Context_(IGC)"),
      RCA = c("Inferred from Reviewed Computational Analysis", "http://wiki.geneontology.org/index.php/Inferred_from_Reviewed_Computational_Analysis_(RCA)"),
      TAS = c("Traceable Author Statement", "http://wiki.geneontology.org/index.php/Traceable_Author_Statement_(TAS)"),
      NAS = c("Non-traceable Author Statement", "http://wiki.geneontology.org/index.php/Non-traceable_Author_Statement_(NAS)"),
      IC = c("Inferred by Curator", "http://wiki.geneontology.org/index.php/Inferred_by_Curator_(IC)"),
      ND = c("No biological Data available", "http://wiki.geneontology.org/index.php/No_biological_Data_available_(ND)_evidence_code"),
      IEA = c("Inferred from Electronic Annotation", "http://wiki.geneontology.org/index.php/Inferred_from_Electronic_Annotation_(IEA)"),
      IED = c("Inferred from Experimental Data", "Inferred from Experimental Data"),
      IAGP = c("Inferred from Association of Genotype and Phenotype", "Inferred from Association of Genotype and Phenotype"),
      IPM = c("Inferred from Phenotype Manipulation", "Inferred from Phenotype Manipulation"),
      QTM = c("Quantitative Trait Measurement", "Quantitative Trait Measurement"),
      IIA = c("Inferred from Imported Annotation", "Inferred from Imported Annotation")
    )
  }

  if (length(evidence) <= 8) {
    colors <- brewer.pal(max(3, length(evidence)), "Dark2")
  } else {
    colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(evidence))
  }

  content <- lapply(1:length(evidence), function(indx) {
    gene <- evidence[[indx]]

    if (is.null(gene)) {
      return(NULL)
    } else {
      shinypanels::panel(
        title = gene[[1]]$symbol, collapsed = T, width = 500, color = colors[indx],
        lapply(gene, function(evidence) {
          cat <- ONTOLOGIES.DEFS[Node_Long == evidence$cf.CatLongUri, data.table::first(Definition)]
          if (!is.character(cat)) {
            cat <- evidence$cf.CatLongUri
          }

          val <- ONTOLOGIES.DEFS[Node_Long == evidence$cf.ValLongUri, data.table::first(Definition)]
          if (!is.character(val)) {
            val <- evidence$cf.ValLongUri
          }

          div(
            p(tags$b("Category: "), cat),
            p(tags$b("Value: "), val),
            p(tags$b("Relationship: "), evidence$relationship),
            p(
              tags$b("Evidence Code: "), paste0(evidence$evidence.Code, " "),
              a(target = "_blank", href = evidenceCodes(evidence$evidence.Code)[2], tags$i(class = "fas fa-question-circle"))
            ),
            p(
              tags$b(paste0(ifelse(evidence$score.Name == "", "Score", evidence$score.Name), ": ")),
              ifelse(is.null(evidence$score.Value), evidence$score.Strength, evidence$score.Value)
            ),
            p(tags$b("Citation: "), a(target = "_blank", href = evidence$citation.Url, evidence$citation.Name))
          )
        })
      )
    }
  })

  renderUI(div(
    class = "layout-container",
    div(
      class = "layout-panels",
      div(class = "app-container", content)
    )
  ))
}
