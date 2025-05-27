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
      shiny::verbatimTextOutput("output_path"),
      shiny::actionButton("show_advanced", "Show advanced options"),
      shiny::conditionalPanel(
        condition = "input.show_advanced % 2 == 1",
        tags$div(style = "margin-top: 10px; margin-bottom: 10px; border: 1px solid #eee; padding: 10px; background: #fafafa;",
          shiny::numericInput("lfcCutoff", "logFC cutoff (default is 1)", value = 1, min = 0, step = 0.1),
          shiny::numericInput("pCutoff", "Adjusted Pvalue (default is 0.05)", value = 0.05, min = 0, max = 1, step = 0.01),
          shiny::selectInput("pAdjustMethod", "Method for padj (default is fdr)",
            choices = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"),
            selected = "fdr"
          ),
          shiny::textInput("palette", "Two Colorhexa codes (default is #4f8832 and #f79c18)", value = "#4f8832,#f79c18"),
          shiny::textInput("rowNamesOfCounts", "The name of the rows in counts.txt", value = "external_gene_name")
        )
      ),
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

