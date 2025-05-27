
app_server <- function(input, output, session) {

    # set up roots for shinyFiles
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyDirChoose(input, "salmon_folder", roots = volumes, session = session)

  pipeline_res <- shiny::eventReactive(input$run, {
    shiny::showNotification("Button clicked!", type = "message")
    # or in console
    message("DEBUG: run observer firing")
    
    # check that the required inputs are provided
    shiny::req(input$species)
    

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
      # salmon
      # parseDirPath will give a named vector; take the first element
      folder <- shinyFiles::parseDirPath(volumes, input$salmon_folder)[[1]]
      # copy ID.txt
      shiny::req(input$id_file)
      base::file.copy(input$id_file$datapath, base::file.path(folder, "ID.txt"))
      data_dir <- folder
    }

    # 3) Run the pipeline
    shiny::withProgress(message="Running analysis…", value=0, {
      shiny::showNotification("Running pipeline…", type = "message")
      output$status <- shiny::renderText("Running pipeline…")
      res <- tryCatch(
        goiExplorer::run_pipeline(data_dir,
                                  dataType = input$dataType,
                                  species = input$species,
                                  parent_outdir = tempfile("goi_out_")),
        error = function(e) {
          shiny::showNotification(paste("Error:", e$message), type="error")
          return(NULL)
        }
      )
    })
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
    shiny::req(pipeline_res()) # this won’t execute until pipeline_res() is non-NULL
    pipeline_res()$terms
  })
  
  output$go_plot <- renderPlot({
    shiny::req(pipeline_res()) # nothing will render until pipeline_res() is ready
    pipeline_res()$plot
  })
  
  shiny::observeEvent(input$ai_ask, {
    shiny::req(base::nzchar(base::Sys.getenv("OPENAI_API_KEY")),
        cancelOutput = TRUE)
    output$ai_answer <- shiny::renderText({
      goiExplorer::ai_agent(input$ai_query)
    })
  })
}