app_ui <- shiny::fluidPage(
  shiny::titlePanel("GOI Explorer"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      # 1) Let user pick counts vs. salmon
      radioButtons(
        "dataType",
        "Input type:",
        choices = c("Counts" = "counts", "Salmon" = "salmon"),
        selected = "counts"
      ),
      
      # 2a) When counts: upload two files
      shiny::conditionalPanel(
        condition = "input.dataType == 'counts'",
        shiny::fileInput("counts_file", "counts.txt", accept = ".txt"),
        shiny::fileInput("gr_file",     "gr.txt",     accept = ".txt")
      ),
      
      # 2b) When salmon: choose folder + upload ID.txt
      shiny::conditionalPanel(
        condition = "input.dataType == 'salmon'",
        shinyFiles::shinyDirButton(
          "salmon_folder", 
          "Select Salmon folder", 
          "Please select a directory"
        ),
        fileInput("id_file", "Upload ID.txt", accept = ".txt")
      ),

      # 3a) Gene of interest
      shiny::textInput("goi", "Gene of interest", value = ""),
      # 3b) Mandatory output directory chooser
      shinyFiles::shinyDirButton(
        "save_dir",
        "Choose output directory",
        "Select a directory"
      ),
      # 4) Run pipeline & show outputs
      shiny::actionButton("run", "Run Pipeline"),
      shiny::verbatimTextOutput("output_path")
    ),
    shiny::mainPanel(
      shiny::verbatimTextOutput("status"),
      shiny::tableOutput("results_table"),
      shiny::plotOutput("go_plot"),
      shiny::textInput("ai_query", "Ask the AI agent"),
      shiny::actionButton("ai_ask", "Send to AI"),
      verbatimTextOutput("ai_answer"),
      shiny::conditionalPanel(
        "input.run > 0",
        shiny::downloadButton("dl", "Download results")
      )
    )
  )
)

