#' Server
#'
#' Query strings are supported as follows:
#' ?genes=[list of genes]&taxa=[human|mouse|rat]&...
server <- function(input, output, session) {
  session$onSessionEnded(function() {
    session$userData$INTERRUPT <- TRUE
  })

  #' Advance the progress bar one step for the client
  #'
  #' @param detail Any text to show in the progress bar
  advanceProgress <- function(detail) {
    setProgress(session$userData$progress + 1, detail)
  }

  #' Set the progress bar to a specified step
  #'
  #' @param progress The step to move to or `NULL` to complete the progress bar
  #' @param detail Any text to show in the progress bar
  #' @param n.steps The maximum number of steps this progress bar can take
  setProgress <- function(progress = NULL, detail = "", n.steps = NULL) {
    if (is.null(progress)) {
      if (is.null(session$userData$progress.bar)) {
        return(NULL)
      }
      progress <- session$userData$progress.bar$getMax()
      n.steps <- progress
    }

    session$userData$progress <- progress

    if (progress == 0) {
      session$userData$progress.bar <- shiny::Progress$new(min = 0, max = n.steps)
      session$userData$progress.bar$set(message = "Working...", detail = detail)
    } else {
      session$userData$progress.bar$set(value = progress, detail = detail)

      # Enable the search and reset buttons after reaching completion
      if (progress >= session$userData$progress.bar$getMax()) {
        session$userData$progress.bar$close()
        shinyjs::enable("search")
        shinyjs::enable("reset")
      }
    }
  }

  # On connect, observe the query string. Search if any genes are specified
  observeEvent(input$LOAD,
    {
      output$results_header <- renderUI(generateResultsHeader({
        HTML(paste0(
          '<div style="margin-bottom: 10px"><h3 style="display: inline">',
          paste0(
            "Search across <span>", length(DATA.HOLDER), " species</span>, <span>",
            CORPUS_STATS$studies, " studies</span>, <span>",
            CORPUS_STATS$comparisons, " condition comparisons</span>, <span>",
            CORPUS_STATS$assays, " assays</span>"
          ),
          "</h3></div>"
        ))
      }))

      query <- getQueryString(session)
      options <- getConfig()

      genes <- query$genes %>% jsonify()
      taxa <- query$taxa %>% jsonify()

      sig <- switch(is.null(query$sig) + 1,
        jsonify(query$sig) %>% as.numeric(),
        options$sig$value
      )
      fc <- switch(is.null(query$fc) + 1,
        jsonify(query$fc) %>% as.numeric(),
        options$fc$value
      )
      categories <- switch(is.null(query$categories) + 1,
        jsonify(query$categories),
        options$categories$value
      )

      # Update UI elements from the query string
      updateSelectizeInput(session, "genes")
      updatePickerInput(session, "taxa", selected = taxa)
      session$sendCustomMessage("querySet", genes)

      updateSelectizeInput(session, "sig", selected = sig)
      updateSliderInput(session, "fc", value = fc)
      updateSelectizeInput(session, "categories", selected = categories)

      for (mName in Filter(function(x) !(x %in% c("sig", "fc", "categories", "taxa")), names(options))) {
        do.update(session, options[[mName]], ifelse(is.null(query[[mName]]), options[[mName]]$value, query[[mName]]))
      }
    },
    once = T,
    ignoreNULL = T
  )

  # Send a message to the client if a file is uploaded
  observeEvent(input$genes.csv, {
    if (!is.null(input$genes.csv)) {
      session$sendCustomMessage("fileUpload", T)
    }
  })

  # Not exactly sure why I implemented it this way. Maybe to make sure it resets?
  output$genes.csv.ui <- renderUI({
    input$reset
    fileInput("genes.csv", "Or Upload CSV", accept = "text/csv", placeholder = "N/A")
  })

  # Reset the search bar
  observeEvent(input$reset, {
    session$sendCustomMessage("queryReset", list())
    session$sendCustomMessage("fileUpload", F)
    updatePickerInput(session, "taxa", selected = getConfig("taxa")$value)
    updateSelectizeInput(session, "method", selected = NULL)
    updateSelectizeInput(session, "sig", selected = NULL)

    options <- getConfig()
    for (mName in Filter(function(x) !(x %in% c("taxa", "method", "sig")), names(options))) {
      do.update(session, options[[mName]], options[[mName]]$value)
    }
  })
  
  # Update gene choices depending on taxon(s). Commented out because only works with selectize.js
  # observeEvent(input$taxa, {
  #   updateSelectizeInput(session, "genes", choices = ALL.GENES[input$taxa], server = TRUE)
  # })
  
  # Search
  observeEvent(input$search, {
    genes <- input$genes
    if (!is.null(input$genes.csv)) {
      genes <- read.csv(input$genes.csv$datapath, header = F)$V1 %>% as.character()
    }

    searchGenes(genes)
  })

  #' Open other (non-table) tabs
  #'
  #' @param tab The tab to open
  #' @param delay A delay (in ms) before opening
  openTab <- function(tab, delay = 0) {
    if (delay > 0) {
      shinyjs::delay(delay, {
        openTab(tab)
      })
    } else if (!is.null(session$userData$plotData)) {
      if (tab == "Gene Info" && is.null(session$userData$genesRendered)) {
        session$userData$genesRendered <- T

        synchronise({
          geneEvidence(
            session$userData$plotData$genes,
            session$userData$plotData$options$taxa$value
          )$
            then(function(evidence) {
            output$results_genes <- generateGenePage(evidence)
          })
        })
      } else if (tab == "Gene Contributions" && is.null(session$userData$contribsRendered)) {
        session$userData$contribsRendered <- T

        output$results_contribs <- generateGeneContribs(session$userData$plotData$conditions, session$userData$plotData$options)
      }
    }
  }

  # Open other tabs
  observeEvent(input$tabs, {
    openTab(input$tabs)
  })

  # Force an updated of the gene view when we get new search results
  observeEvent(input$UPDATED, {
    session$userData$genesRendered <- NULL
    session$userData$contribsRendered <- NULL

    # Delay if we're opening the gene info tab (I forget why)
    openTab(input$tabs, 100 * (input$tabs == "Gene Info"))
  })

  # Populate with random genes
  observeEvent(input$RANDOM_GENES, {
    updateSelectizeInput(session, "genes", options = list(persist = F, create = T, createOnBlur = T))
    while (length(genes <- DATA.HOLDER[[ifelse(length(input$taxa) > 1, "human", input$taxa)]]@gene.meta[sample(1:.N, sample(1:10, 1))] %>%
      as.data.frame() %>% .[, sample(c(1, 3, 4, 6), 1)] %>%
      .[. != "None"]) == 0) {
    }

    genes <- genes %>%
      {
        paste0("[", paste0(., collapse = ","), "]")
      } %>%
      jsonify()

    session$sendCustomMessage("queryReset", genes)
  })

  # Open plot
  observeEvent(input$plotData, {
    shiny::showModal(modalDialog(
      easyClose = T,
      size = "l",
      footer = NULL,
      generatePlotContainer()
    ), session)

    sendPlot <- function() {
      output$plot <- renderPlotly({
        generateResultsPlot(
          session$userData$plotData$genes,
          session$userData$plotData$conditions[`Ontology Steps` == 0] %>%
            .[1:5, paste0(cf.BaseLongUri, " vs. ", cf.ValLongUri)], # TODO selection strategy
          session$userData$plotData$exprInfo,
          session$userData$plotData$options,
          input$plot_taxa, input$plot_genes, input$plot_conditions, input$plot_type, input$plot_data
        )
      })
    }

    # We don't have plot data; we need to make it
    if (is.null(session$userData$plotData$exprInfo)) {
      if (!is.null(session$userData$plotData$conditions)) {
        # Give the UI a tick to update and then fetch the data.
        shinyjs::delay(1, {
          synchronise({
            mExpData <- rbindlist(lapply(
              session$userData$plotData$options$taxa$value,
              function(i) DATA.HOLDER[[i]]@experiment.meta[, .(rsc.ID, ee.ID, ee.NumSample)]
            ))

            # Choose experiments to display in the plot
            getTags(session$userData$plotData$options$taxa$value) %>% # TODO Selection strategy
              .[cf.Cat %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.Cat)]) &
                cf.BaseLongUri %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.BaseLongUri)]) &
                cf.ValLongUri %in% (session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, unique(cf.ValLongUri)])] %>%
              unique() %>%
              merge(mExpData, by = c("rsc.ID", "ee.ID"), sort = F, allow.cartesian = T) %>%
              {
                if (nrow(.) == 0) {
                  return(NULL)
                }

                # Select some samples to queue for gene expression visualization
                while ((tmp <- .[sample(1:nrow(.))])[1, ee.NumSample] > getOption("max.gemma")) {
                }

                tmp[cumsum(ee.NumSample) <= getOption("max.gemma")]
              } %>%
              {
                geneExpression(
                  .[, unique(ee.ID)],
                  .[, unique(rsc.ID)],
                  session$userData$plotData$options$taxa$value,
                  session$userData$plotData$genes
                )$then(function(value) {
                  # Delegate the plot renderer
                  session$userData$plotData$exprInfo <- value

                  # Send options to the UI
                  updatePickerInput(session, "plot_taxa",
                    choices = as.list(session$userData$plotData$options$taxa$value),
                    selected = session$userData$plotData$options$taxa$value[1]
                  )

                  if (length(session$userData$plotData$options$taxa$value) > 1) {
                    mGenes <- unique(session$userData$plotData$genes[, gene.realName])
                  } else {
                    mGenes <- DATA.HOLDER[[session$userData$plotData$options$taxa$value]]@gene.meta[entrez.ID %in% unique(session$userData$plotData$genes), gene.Name]
                  }

                  updatePickerInput(session, "plot_genes",
                    choices = as.list(intersect(mGenes, rownames(value$expr))),
                    selected = intersect(mGenes, rownames(value$expr))
                  )
                  updatePickerInput(session, "plot_conditions", # TODO Selection strategy
                    choices = as.list(session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, " vs. ", cf.ValLongUri)]),
                    selected = data.table::first(session$userData$plotData$conditions[`Ontology Steps` == 0] %>% .[1:5, paste0(cf.BaseLongUri, " vs. ", cf.ValLongUri)])
                  )

                  # Send the plot to the UI
                  sendPlot()
                })
              }
          })
        })
      }
    } else { # We already have plot data; send it
      sendPlot()
    }
  })

  #' Display a failure message
  #'
  #' Ends the search protocol with a failure message
  endFailure <- function() {
    setProgress()
    output$results_header <- renderUI({
      generateResultsHeader(HTML('<h3 style="margin-top: 0;">No results</h2>'))
    })
    output$results <- NULL
  }

  #' Make suggestions on taxa/genes for users
  #'
  #' @param original The original gene query
  #' @param taxon The original taxa query
  #' @param suggestions Suggestions for closest gene matches or otherwise mapped gene IDs
  makeSuggestions <- function(original, taxon, suggestions) {
    if (!is.null(suggestions)) {
      # If we searched only one taxon, make sure it looks appropriate
      if (length(taxon) == 1) {
        suggested <- lapply(getConfig("taxa")$core, function(t) {
          i <- 0
          oGenes <- NA
          tmp <- tidyGenes(original, t)
          if (!is.null(tmp)) {
            i <- length(tmp$genes)
            oGenes <- tmp$genes
          }

          list(
            n = i,
            genes = oGenes,
            taxon = t
          )
        }) %>%
          {
            # If a different taxon matches our gene set better, suggest that taxon
            suggestedTaxon <- lapply(., "[[", "n") %>% which.max()
            list(taxon = .[[suggestedTaxon]]$taxon, genes = .[[suggestedTaxon]]$genes)
          }

        if (suggested$taxon != taxon) {
          # Update the UI
          output$results_suggestions <- renderUI({
            fluidRow(
              class = "info-text",
              column(12, HTML(paste0(
                "Did you mean to search ", paste0("<a search taxon='", suggested$taxon, "'>", suggested$taxon, "</a>"),
                " experiments? Otherwise, try using <a search genes='[", paste0(paste0('"', tidyGenes(suggested$genes, c(taxon, suggested$taxon))[taxon != suggested$taxon, gene.realName], '"'), collapse = ","), "]'>", taxon, " orthologs", "</a> of these ", suggested$taxon, " genes."
              )))
            )
          })
          shinyjs::delay(100, {
            shinyjs::runjs("loadExamples();")
          })

          return(NULL)
        }
      }

      # Substitute the corrected genes in place for a new search
      revisedSearch <- original
      for (i in names(suggestions)) {
        revisedSearch <- gsub(i, suggestions[i], revisedSearch, fixed = T)
      }

      # Update the UI
      output$results_suggestions <- renderUI({
        fluidRow(class = "info-text", column(12, HTML(paste0("Did you mean: ", paste0("<a search genes='[", paste0(paste0('"', revisedSearch, '"'), collapse = ","), "]'>", paste0(suggestions, collapse = ", "), "</a>"), "?"))))
      })
      shinyjs::delay(100, {
        shinyjs::runjs("loadExamples();")
      })
    }
  }

  #' Display an empty message
  #'
  #' Ends the search protocol with no results
  endEmpty <- function() {
    setProgress()
    output$results_header <- renderUI({
      generateResultsHeader("No conditions found in scope.")
    })
    output$results <- NULL
  }

  #' Display a success message
  #'
  #' Ends the search protocol successfully
  #'
  #' @param genes The genes that were searched
  #' @param experiments The experiment rankings that were obtained
  #' @param conditions The condition rankings that were obtained
  #' @param options Any additional options that were passed
  endSuccess <- function(genes, experiments, conditions, taxa) {
    if (is.data.table(experiments)) {
      exps <- experiments$rn
    } else {
      exps <- lapply(experiments, "[[", "rn") %>% unlist()
    }

    tmp <- rbindlist(lapply(taxa, function(i) {
      DATA.HOLDER[[i]]@experiment.meta[rsc.ID %in% exps, .(rsc.ID, ee.ID, ee.Name, ee.NumSample, ef.IsBatchConfounded)]
    }))

    # Generate the results header
    if (exists("output")) {
      studies <- length(unique(tmp$ee.ID))
      assays <- tmp[!duplicated(ee.ID), sum(ee.NumSample)]

      n_exp <- length(exps)

      output$results_header <- renderUI({
        generateResultsHeader(HTML(paste0(
          '<div style="margin-bottom: 10px"><h2 style="display: inline">Examined ',
          format(n_exp, big.mark = ","), " condition comparison", ifelse(n_exp > 1, "s", ""),
          " (", format(studies, big.mark = ","), " studies, ", format(assays, big.mark = ","), " assays) for ",
          ifelse(nrow(genes) == 1, genes[, gene.Name],
            paste0(
              '<span data-toggle="tooltip" data-placement="top" title="',
              paste0(unique(genes[, gene.Name]), collapse = ", "), '">',
              nrow(genes), " gene", ifelse(nrow(genes) > 1, "s", ""), "</span>"
            )
          ),
          ifelse(length(taxa) > 1,
            paste0(
              ' across <span data-toggle="tooltip" data-placement="top" title="', paste0(taxa, collapse = ", "), '">',
              length(taxa), " taxa</span>"
            ), ""
          ),
          '</h2><span class="timestamp">in ',
          format(difftime(session$userData$endTime, session$userData$startTime), digits = 3), ".</span></span>"
        )))
      })
    }

    # Some cleaning
    conditions[, `Condition Comparison` := stringi::stri_c("<b>", cf.BaseLongUri, "</b> vs. <b>", cf.ValLongUri, "</b>")]

    advanceProgress("Cross-linking")

    # Associating experiments with tags
    tmp <- tmp %>%
      merge(getTags(taxa, exps), sort = F) %>%
      .[, N := length(unique(ee.ID)), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)] %>%
      .[N < 0.03 * nrow(tmp)] # Get rid of contrasts that overlap in more than 3% experiments

    # Add links to Gemma (short form here to not create long strings which can be memory hungry)
    # these are dealt with in JS when the table pages show
    tmp[
      , Evidence := ifelse(N == 1,
        stringi::stri_c('<span data-id="', ee.ID[1], '"', ifelse(any(ef.IsBatchConfounded), stringi::stri_c(' data-conf="', stringi::stri_c(which(ef.IsBatchConfounded)), '"'), ""), ">1 Experiment</span>"),
        {
          dedup <- !duplicated(ee.ID)
          stringi::stri_c(
            '<span data-id="', stringi::stri_c(ee.ID[dedup], collapse = ","),
            '" data-ee="', stringi::stri_c(ee.Name[dedup], collapse = ","), '"',
            ifelse(any(ef.IsBatchConfounded), stringi::stri_c(' data-conf="', stringi::stri_c(which(ef.IsBatchConfounded)), '"'), ""),
            ">", N, " Experiments</span>"
          )
        }
      ),
      .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)
    ]
    
    tmp[, EvidencePlain := {
      dedup <- !duplicated(ee.ID)
      stringi::stri_c(ee.Name[dedup], collapse = ",")
      },.(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]

    tmp[, .(cf.Cat, cf.BaseLongUri, cf.ValLongUri, N, Evidence, EvidencePlain)] %>%
      unique() %>%
      merge(conditions, by = c("cf.Cat", "cf.BaseLongUri", "cf.ValLongUri"), sort = F) %>%
      setnames(c("stat", "score", "distance"), c("Effect Size", "Test Statistic", "Ontology Steps"))
  }

  #' Handle Search
  #'
  #' Call the processing function to handle searching and update the display.
  #'
  #' @param genes A character vector of entrez IDs
  #' @param options The search options
  handleSearch <- function(genes, rawGenes, options) {
    advanceProgress("Ranking experiments")

    # For some reason genes becomes a Promise and this resolves it?
    # TODO this is dirty and needs to be re-examined
    genes
    
    future(
      {
        # First run a search asynchronously
        # print(options$mfx$value)
        if (!is.null(genes) && length(options$taxa$value) > 1 && "data.table" %in% class(genes)) {
          lapply(options$taxa$value, function(t) {
            mOp <- options
            mOp$taxa$value <- t
            vsmSearch(genes[taxon == t, entrez.ID], taxa = mOp$taxa$value, confounds = mOp$confounds$value, filter = mOp$filter$value, mfx = FALSE, geeq = mOp$geeq$value, p_threshold = mOp$pv$value)
          })
        } else if (length(options$taxa$value) == 1) {
          vsmSearch(genes$genes, taxa = options$taxa$value, confounds = options$confounds$value, filter = options$filter$value,  mfx = FALSE, geeq = options$geeq$value, p_threshold = options$pv$value)
        }
      },
      globals = "DATA.HOLDER",
      seed = T
    ) %...>% (function(experiments) {
      if (!is.null(genes) && !("data.table" %in% class(genes))) {
        makeSuggestions(rawGenes, options$taxa$value, genes$suggestions)
      }

      if (!is.null(session$userData$INTERRUPT)) {
        return(NULL)
      }

      if (all(sapply(experiments, is.null))) {
        endFailure()
      } else {
        advanceProgress("Enriching")

        future(
          {
            # Do enrichments asynchronously too.
            if (length(options$taxa$value) == 1) {
              copy(enrich(experiments, taxa = options$taxa$value, dist = options$dist$value,categories = options$categories$value))
            } else {
              # TODO it's maybe reasonable to mclapply this as long as the
              # multicore workers (@seealso dependencies.R) is sufficiently low
              lapply(1:length(options$taxa$value), function(t) {
                mOp <- options
                mOp$taxa$value <- options$taxa$value[t]
                mSearchable <- experiments[[t]]

                enrich(mSearchable, taxa = mOp$taxa$value, dist = mOp$dist$value,categories = mOp$categories$value) %>%
                  setnames(genes[taxon == options$taxa$value[t], gene.realName],
                    genes[taxon == options$taxa$value[t], identifier],
                    skip_absent = T
                  ) # TODO ?
              }) %>%
                rbindlist(fill = T) %>%
                reorderTags3() %>%
                .[, lapply(.SD, mean, na.rm = T), .(cf.Cat, cf.BaseLongUri, cf.ValLongUri)]
            }
          },
          globals = c("CACHE.BACKGROUND", "NULLS"),
          seed = TRUE
        ) %...>% (function(conditions) {
          if (!is.null(session$userData$INTERRUPT)) {
            return(NULL)
          }

          session$userData$endTime <- Sys.time()

          if (is.null(conditions)) {
            endEmpty()
          } else {
            if (length(options$taxa$value) > 1) {
              geneInfo <- genes %>%
                copy() %>%
                setnames("identifier", "gene.Name")
              mGenes <- genes %>% copy()
            } else {
              geneInfo <- DATA.HOLDER[[options$taxa$value]]@gene.meta %>%
                .[, .(.I, entrez.ID, gene.Name)] %>%
                .[entrez.ID %in% genes$genes, .(I, entrez.ID, gene.Name)]
              mGenes <- genes$genes
            }


            conditions <- endSuccess(geneInfo, experiments, conditions, options$taxa$value)

            
            if (!is.null(session$userData$INTERRUPT)) {
              return(NULL)
            }

            getPercentageStat <- function(x, n = 1){
              x / n
            }
            conditions[,'Test Statistic'] <- apply(conditions[,'Test Statistic'], 2, getPercentageStat, n = nrow(geneInfo))
            output$results <- generateResults(conditions)
            # Prepare some plotting information.
            session$userData$plotData <- list(
              options = options,
              genes = mGenes,
              conditions = conditions %>% copy()
            )

            # Prepare a data download handler
            output$dataDownload <- downloadHandler(function() {
              paste0("enrichment-", Sys.Date(), ".csv")
            }, function(file) {
              write.csv(session$userData$plotData$conditions %>%
                .[, !c("Condition Comparison", "Evidence")] %>%
                data.table::setorder(-`Test Statistic`, `Ontology Steps`) %>%
                data.table::setnames(c("cf.Cat", "cf.BaseLongUri", "cf.ValLongUri","EvidencePlain"), c("Category", "Baseline", "Value","Evidence")), file)
            }, "text/csv")
          }
        })
      }
    })
  }

  #' Search Genes
  #'
  #' Begin the gene search process. Cleans the parameters as necessary and calls a function to
  #' perform the actual search and update the display.
  #'
  #' @param genes A character vector of genes (entrez ID, ensembl ID, name or keyword)
  #' @param signature A DE signature to search
  #' @param taxa The taxon ([human|mouse|rat])
  #' @param update Whether or not to update the browser's query string.
  searchGenes <- function(genes = input$genes,
                          signature = input$signature,
                          taxa = input$taxa,
                          update = T) {
    output$results_suggestions <- NULL

    # These will be re-enabled when the progress bar completes (@seealso setProgress)
    shinyjs::disable("search")
    shinyjs::disable("reset")

    session$userData$startTime <- as.POSIXct(Sys.time())
    setProgress(0, "Validating input", n.steps = getOption("max.progress.steps"))

    # Populate options
    options <- lapply(names(getConfig()), function(x) {
      ret <- getConfig(x)
      if (x == "sig") {
        v <- signature
      } else {
        v <- input[[x]]
      }

      if (!is.null(v)) {
        ret$value <- v
      }
      ret
    }) %>% `names<-`(names(getConfig()))

    # Update the query string
    if (update) {
      query <- paste0("?genes=", switch(min(2, length(genes)),
        genes,
        paste0("[", paste0(genes, collapse = ","), "]")
      ))

      # TODO have to look into categories and subsets, might need to be treated specially
      extras <- lapply(names(getConfig()), function(x) {
        if (x == "sig") {
          v <- signature
        } else {
          v <- input[[x]]
        }

        if (!is.null(v) && v != "" && (!isTRUE(all.equal(v, getConfig(key = x)$value)) || length(v) != length(getConfig(key = x)$value))) {
          if (length(v) > 1) {
            paste0(x, "=", paste0("[", paste0(v, collapse = ","), "]"))
          } else {
            paste0(x, "=", v)
          }
        }
      }) %>%
        unlist() %>%
        paste0(collapse = "&")

      if (nchar(extras) > 0) {
        query <- paste0(query, "&", extras)
      }

      updateQueryString(query, "push")
    }

    if (is.null(taxa)) taxa <- options$taxa$value

    # Done processing, handle the search.
    # TODO should prevent signature from working for length(taxa) > 1
    # TODO This { . } is a weird way to prevent lazy evaluation. Probably a better way
    tidyGenes(genes, taxa) %>%
      {
        .
      } %>%
      handleSearch(genes, options)
  }
}
