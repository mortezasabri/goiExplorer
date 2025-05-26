
app_server <- function(input, output, session) {
  pipeline_res <- shiny::eventReactive(input$run, {
    shiny::showNotification("Button clicked!", type = "message")
    # or in console
    message("DEBUG: run observer firing")
    
    # Block until we have the two required inputs
    shiny::req(input$counts_file, input$gr_file, input$species)
    
    # create a dedicated temporary directory
    data_dir <- base::tempfile("goi_data_")
    base::dir.create(data_dir)
    # copy the uploaded files in with the expected names
    file.copy(input$counts_file$datapath, base::file.path(data_dir, "counts.txt"))
    file.copy(input$gr_file$datapath,     base::file.path(data_dir, "gr.txt"))
    
    shiny::withProgress(message="Running analysis…", value=0, {
      output$status <- shiny::renderText("Running pipeline…")
      res <- tryCatch(
        #counts <- read.csv(input$infile$datapath)[[1]]
        goiExplorer::run_pipeline(data_dir,
                                  dataType      = input$dataType,
                                  species = input$species,
                                  parent_outdir = withr::local_tempdir()),
        error = function(e) {
          shiny::showNotification(paste("Error:", e$message), type="error")
          return(NULL)
        }
      )
      # output$status <- shiny::renderText("Done.")
      # res
    })
    
  })
  
  output$status <- shiny::renderText({
    shiny::req(pipeline_res())
    "Finished!"
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


