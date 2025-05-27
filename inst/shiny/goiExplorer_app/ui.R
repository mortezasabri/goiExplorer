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
      shiny::textInput("goi", "Gene of interest (e.g. CYLD)", value = ""),
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
      tags$div(
        style = "font-size: 1.2em; font-weight: bold; margin-bottom: 10px; color: #2c3e50;",
        shiny::verbatimTextOutput("output_path")
      ),
      shiny::tabsetPanel(
        id = "plot_tabs",
        tabPanel("Barplot",   shiny::plotOutput("barplot")),
        tabPanel("Boxplot",   shiny::plotOutput("boxplot")),
        tabPanel("Countplot", shiny::plotOutput("countplot")),
        tabPanel("MA Plot",   shiny::plotOutput("ma_plot")),
        tabPanel("Volcano",   shiny::plotOutput("volcano_plot")),
        tabPanel(
          "Differentially Expressed Genes",
          DT::dataTableOutput("degs_table")
        )
      ),
      shiny::tableOutput("goi_des"),
      shiny::verbatimTextOutput("goi_entrez"),
      shiny::tableOutput("goi_test"),
      shiny::tableOutput("results_table")
    )
  )
)

