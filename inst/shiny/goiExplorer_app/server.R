`%||%` <- function(x, y) if (is.null(x)) y else x

app_server <- function(input, output, session) {

    # set up roots for shinyFiles
    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
    shinyFiles::shinyDirChoose(input, "salmon_folder", roots = volumes, session = session)
    shinyFiles::shinyDirChoose(input, "save_dir", roots = volumes, session = session)

    # where the current run wrote its output; also drives the download buttons
    output_dir_rv <- shiny::reactiveVal(NULL)

    pipeline_res <- shiny::eventReactive(input$run, {
        message("DEBUG: run observer firing")

        # Safely parse the Salmon folder selection
        salmon_dir_sel <- shinyFiles::parseDirPath(volumes, input$salmon_folder)
        salmon_folder  <- if (length(salmon_dir_sel) >= 1 && nzchar(salmon_dir_sel[1])) {
            salmon_dir_sel[[1]]
        } else {
            NULL
        }

        # check that the required inputs are provided
        shiny::req(input$goi)  # require gene of interest

        # 1) Assemble a clean temp subfolder
        data_dir <- base::tempfile("goi_data_")
        base::dir.create(data_dir)

        # 2) Branch on dataType
        if (input$dataType == "counts") {
            shiny::req(input$counts_file, input$gr_file)
            # copy uploaded files into data_dir
            base::file.copy(input$counts_file$datapath, base::file.path(data_dir, "counts.txt"))
            base::file.copy(input$gr_file$datapath,     base::file.path(data_dir, "gr.txt"))
        } else {
            # salmon branch: enforce folder selection
            shiny::req(input$salmon_folder)
            data_dir <- salmon_folder
            # copy ID.txt into place
            shiny::req(input$id_file)
            base::file.copy(input$id_file$datapath, file.path(data_dir, "ID.txt"))
        }

        # Parse the mandatory output directory
        save_dir_input <- input$save_dir
        if (is.null(save_dir_input)) stop("No directory selected.")

        # shinyFiles returns a list with 'root' and 'path'
        if (!is.null(save_dir_input$root) && !is.null(save_dir_input$path)) {
            root_path <- volumes[[save_dir_input$root]]
            sub_path <- do.call(file.path, as.list(save_dir_input$path))
            save_folder <- file.path(root_path, sub_path)
        } else {
            save_folder <- NULL
        }
        shiny::req(save_folder)
        output_dir <- normalizePath(save_folder, mustWork = FALSE)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        # 3) Run the pipeline
        shiny::withProgress(message = "Running analysis…", value = 0, {
            shiny::incProgress(0.2, detail = "Initializing data…")
            shiny::showNotification("Running pipeline…", type = "message")
            # Parse palette as a character vector
            palette_vals <- trimws(unlist(strsplit(input$palette, ",")))
            result <- tryCatch(
                goiExplorer::run_pipeline(
                  input         = data_dir,
                  dataType      = input$dataType,
                  goi           = input$goi,
                  parent_outdir = output_dir,
                  lfcCutoff     = input$lfcCutoff,
                  pCutoff       = input$pCutoff,
                  pAdjustMethod = input$pAdjustMethod,
                  palette       = palette_vals,
                  rowNamesOfCounts = input$rowNamesOfCounts,
                  extra_plots   = isTRUE(input$extra_plots),
                  report        = isTRUE(input$report)
                ),
                error = function(e) {
                    shiny::showNotification(paste0(
                      "Failed to query Ensembl: ", e$message,
                      ".\nTry again later or set a different mirror, e.g.:\n",
                      "biomaRt::useMart(host='www.ensembl.org', ",
                      "biomart='ENSEMBL_MART_ENSEMBL', dataset=input$species)"
                    ), type="error", duration = NULL)
                    return(NULL)
                }
            )
            shiny::incProgress(0.8, detail = "Finalizing results…")
        })
        if (!is.null(result)) {
            showNotification("Pipeline complete!", type = "message")
            output_dir_rv(output_dir)
        }

        result
    })

    output$status <- renderText({
      if (!isTruthy(input$run)) {
        "Waiting for you to click Run…"
      } else if (is.null(pipeline_res())) {
        "Pipeline failed or still running…"
      } else {
        "Pipeline complete!"
      }
    })

    output$output_path <- shiny::renderText({
      d <- output_dir_rv()
      if (is.null(d)) "" else d
    })


    # ---- headline numbers --------------------------------------------------
    output$summary_cards <- shiny::renderUI({
      res <- pipeline_res()
      if (is.null(res)) return(NULL)

      card <- function(value, label) {
        shiny::tags$div(class = "goi-stat",
          shiny::tags$div(class = "v", value),
          shiny::tags$div(class = "l", label)
        )
      }
      degs <- res$degs
      n_deg <- if (is.null(degs)) 0 else nrow(degs)
      n_up  <- if (n_deg > 0) sum(degs$log2FoldChange > 0, na.rm = TRUE) else 0
      n_dn  <- if (n_deg > 0) sum(degs$log2FoldChange < 0, na.rm = TRUE) else 0

      cards <- list(
        card(if (is.null(res$res.df)) "-" else format(nrow(res$res.df), big.mark = ","), "genes tested"),
        card(format(n_deg, big.mark = ","), "significant"),
        card(n_up, "up"),
        card(n_dn, "down")
      )
      if (!is.null(res$res.df) && !is.null(res$goi) && res$goi %in% rownames(res$res.df)) {
        row <- res$res.df[res$goi, ]
        cards <- c(cards, list(
          card(round(row$log2FoldChange[1], 2), paste0(res$goi, " log2FC")),
          card(signif(row$padj[1], 3), paste0(res$goi, " padj"))
        ))
      }
      shiny::tags$div(class = "goi-stats", cards)
    })


    # ---- downloads ---------------------------------------------------------
    output$download_ui <- shiny::renderUI({
      shiny::req(pipeline_res())
      shiny::tagList(
        shiny::downloadButton("download_zip", "Download all results (.zip)"),
        if (!is.null(pipeline_res()$reportPath)) {
          shiny::downloadButton("download_report", "Download HTML report")
        }
      )
    })

    output$download_zip <- shiny::downloadHandler(
      filename = function() {
        paste0("goiExplorer_", pipeline_res()$goi, "_",
               format(Sys.Date(), "%Y%m%d"), ".zip")
      },
      content = function(file) {
        d <- output_dir_rv()
        shiny::req(d)
        # zip relative to the output directory so the archive has no
        # machine-specific leading path
        owd <- setwd(d)
        on.exit(setwd(owd), add = TRUE)
        tryCatch(
          utils::zip(zipfile = file, files = list.files(".", recursive = TRUE)),
          error = function(e) {
            shiny::showNotification(
              paste0("Could not build the archive (is a `zip` command available?): ", e$message),
              type = "error", duration = NULL
            )
          }
        )
      },
      contentType = "application/zip"
    )

    output$download_report <- shiny::downloadHandler(
      filename = function() paste0("goiExplorer_report_", pipeline_res()$goi, ".html"),
      content = function(file) {
        p <- pipeline_res()$reportPath
        shiny::req(p, file.exists(p))
        file.copy(p, file, overwrite = TRUE)
      },
      contentType = "text/html"
    )


    # ---- plots -------------------------------------------------------------
    # each tab pulls a ggplot straight out of the pipeline result
    plot_slots <- list(
      barplot          = "Barplot",
      boxplot          = "Boxplot",
      countplot        = "Countplot",
      ma_plot          = "plotMA",
      volcano_plot     = "Volcanoplot",
      volcano_labelled = "VolcanoLabelled",
      deg_heatmap      = "DEGHeatmap",
      rank_plot        = "RankPlot",
      correlated_plot  = "CorrelatedGenesPlot",
      disease_plot     = "DiseasePlot",
      pca_plot         = "PCAplot",
      distance_plot    = "SampleDistancePlot",
      libsize_plot     = "LibrarySizePlot",
      pvalue_plot      = "PvalueHistogram",
      dispersion_plot  = "DispersionPlot"
    )
    for (id in names(plot_slots)) {
      local({
        out_id <- id
        slot <- plot_slots[[id]]
        output[[out_id]] <- shiny::renderPlot({
          res <- pipeline_res()
          shiny::req(res)
          g <- res[[slot]]
          shiny::validate(shiny::need(
            !is.null(g),
            "This plot was not produced for this run (it may have been switched off, or the data it needs was unavailable)."
          ))
          g
        })
      })
    }


    # ---- tables ------------------------------------------------------------
    output$goi_des <- shiny::renderTable({
      shiny::req(pipeline_res())
      des <- pipeline_res()$des
      if (is.null(des) || nrow(des) == 0) return(NULL)
      des
    })

    output$goi_test <- shiny::renderTable({
      shiny::req(pipeline_res())
      test <- pipeline_res()$test
      if (is.null(test) || nrow(test) == 0) return(NULL)
      test
    })

    output$run_summary <- shiny::renderText({
      shiny::req(pipeline_res())
      goiExplorer::summarise_run(pipeline_res())
    })

    output$degs_table <- DT::renderDataTable({
      shiny::req(pipeline_res())
      degs <- pipeline_res()$degs
      if (is.null(degs) || nrow(degs) == 0) return(NULL)

      sci_cols <- c("lfcSE", "pvalue", "padj")
      skip_cols <- c("entrezgene_id", "start_position", "end_position", "strand")

      for (col in names(degs)) {
        if (is.numeric(degs[[col]]) && !(col %in% skip_cols)) {
          if (col %in% sci_cols) {
            degs[[col]] <- formatC(degs[[col]], format = "e", digits = 3)
          } else {
            degs[[col]] <- round(degs[[col]], 1)
          }
        }
      }

      DT::datatable(
        degs,
        options = list(
          pageLength = 10,
          lengthMenu = c(5, 10, 25, 50),
          searchHighlight = TRUE,
          scrollX = TRUE   # Enable horizontal scrolling
        ),
        filter = "top",
        rownames = FALSE
      )
    })


    # ---- AI chat -----------------------------------------------------------
    chat <- shiny::reactiveVal(list())   # list(role=, content=)

    # the example links just prefill the question box
    examples <- list(
      ai_ex1 = "Interpret this run for me: what happened to the gene of interest, and how convincing is it?",
      ai_ex2 = "Is my gene of interest a strong hit here, or only just significant? Compare it to the rest of the results.",
      ai_ex3 = "What do the top differentially expressed genes have in common? Any recurring pathway or process?",
      ai_ex4 = "Write a figure legend for the boxplot of my gene of interest, suitable for a paper."
    )
    for (ex in names(examples)) {
      local({
        id <- ex
        txt <- examples[[ex]]
        shiny::observeEvent(input[[id]], {
          shiny::updateTextAreaInput(session, "ai_query", value = txt)
        })
      })
    }

    shiny::observeEvent(input$ai_clear, chat(list()))

    shiny::observeEvent(input$ai_ask, {
      q <- trimws(input$ai_query %||% "")
      shiny::req(nzchar(q))

      history <- chat()
      chat(c(history, list(list(role = "user", content = q))))
      shiny::updateTextAreaInput(session, "ai_query", value = "")

      answer <- shiny::withProgress(message = "Asking the model…", value = 0.5, {
        goiExplorer::ai_agent(
          prompt     = q,
          res_output = pipeline_res_safe(),
          history    = history,
          provider   = input$ai_provider,
          model      = if (nzchar(input$ai_model %||% "")) input$ai_model else NULL,
          api_key    = if (nzchar(input$ai_key %||% "")) input$ai_key else NULL
        )
      })

      chat(c(chat(), list(list(role = "assistant", content = answer))))
    })

    # the chat works with or without a finished run
    pipeline_res_safe <- function() {
      tryCatch(if (isTruthy(input$run)) pipeline_res() else NULL, error = function(e) NULL)
    }

    output$ai_chat <- shiny::renderUI({
      msgs <- chat()
      if (length(msgs) == 0) {
        return(shiny::tags$div(
          class = "hint",
          shiny::tags$p("Run the pipeline first, then ask anything about the result."),
          shiny::tags$p("Set ANTHROPIC_API_KEY or OPENAI_API_KEY before starting the app, or paste a key on the left.")
        ))
      }
      shiny::tags$div(
        class = "chat",
        lapply(msgs, function(m) {
          shiny::tags$div(
            class = paste("msg", if (m$role == "user") "user" else "bot"),
            shiny::tags$div(class = "who", if (m$role == "user") "you" else "assistant"),
            m$content
          )
        })
      )
    })
}

