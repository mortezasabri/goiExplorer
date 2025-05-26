

app_ui <- shiny::fluidPage(
  shiny::titlePanel("GOI Explorer"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      radioButtons(
        "dataType",
        "Input type:",
        choices = c("counts", "salmon"),
        selected = "counts"
      ),
      #shiny::fileInput("infile", "counts table (txt)"),
      shiny::fileInput("counts_file", "Upload counts.txt", accept = ".txt"),
      shiny::fileInput("gr_file", "Upload gr.txt",  accept = ".txt"),
      shiny::textInput("ensemblSpecies", "ensembl name for the species (e.g. hsapiens_gene_ensembl)"),
      shiny::actionButton("run", "Run Pipeline")
      # 
      # 
      # shiny::verbatimTextOutput("ai_answer")
    ),
    shiny::mainPanel(
      shiny::verbatimTextOutput("status"),
      shiny::tableOutput("results_table"),
      shiny::plotOutput("plot"),
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

