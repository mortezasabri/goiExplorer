app_server <- function(input, output, session) {

    # set up roots for shinyFiles
    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
    shinyFiles::shinyDirChoose(input, "salmon_folder", roots = volumes, session = session)
    shinyFiles::shinyDirChoose(input, "save_dir", roots = volumes, session = session)

    pipeline_res <- shiny::eventReactive(input$run, {
        shiny::showNotification("Button clicked!", type = "message")
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
            output$status <- shiny::renderText("Running pipeline…")
            result <- tryCatch(
                goiExplorer::run_pipeline(
                  input         = data_dir,
                  dataType      = input$dataType,
                  goi           = input$goi,
                  parent_outdir = output_dir
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
        }

        # Expose the chosen output path
        output$output_path <- shiny::renderText({ output_dir })

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
  
    output$results_table <- shiny::renderTable({
        shiny::req(pipeline_res())
        pipeline_res()$terms
    })
  
    output$go_plot <- renderPlot({
        shiny::req(pipeline_res())
        pipeline_res()$plot
    })
  
    # Removed the observeEvent for opening the output folder
  
    # AI query reactive
    ai_response <- shiny::eventReactive(input$ai_ask, {
        req(nzchar(input$ai_query))
        goiExplorer::ai_agent(input$ai_query)
    })

    # Render AI answer
    output$ai_answer <- shiny::renderText({ ai_response() })
}