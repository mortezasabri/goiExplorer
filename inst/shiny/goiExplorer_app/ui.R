app_ui <- shiny::fluidPage(
  shiny::tags$head(shiny::tags$style(shiny::HTML("
    .goi-stats { display:flex; flex-wrap:wrap; gap:10px; margin:10px 0 18px; }
    .goi-stat { flex:1 1 120px; background:#f8f9fa; border:1px solid #e4e7eb;
      border-radius:8px; padding:10px 12px; }
    .goi-stat .v { font-size:1.4em; font-weight:650; color:#2c3e50; }
    .goi-stat .l { font-size:.75em; color:#7b8794; text-transform:uppercase;
      letter-spacing:.05em; }
    .chat { max-height:460px; overflow-y:auto; padding:4px; margin-bottom:12px; }
    .msg { padding:10px 13px; border-radius:10px; margin-bottom:9px;
      white-space:pre-wrap; line-height:1.5; }
    .msg.user { background:#eef4ea; border-left:3px solid #4f8832; }
    .msg.bot { background:#fdf4e6; border-left:3px solid #f79c18; }
    .msg .who { font-size:.72em; text-transform:uppercase; letter-spacing:.06em;
      color:#7b8794; margin-bottom:4px; }
    .hint { color:#7b8794; font-size:.9em; }
  "))),

  shiny::titlePanel("GOI Explorer"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      # 1) Let user pick counts vs. salmon
      shiny::radioButtons(
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
        shiny::fileInput("id_file", "Upload ID.txt", accept = ".txt")
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
      shiny::actionButton("run", "Run Pipeline", class = "btn-primary"),
      shiny::actionButton("show_advanced", "Show advanced options"),
      shiny::conditionalPanel(
        condition = "input.show_advanced % 2 == 1",
        shiny::tags$div(
          style = "margin-top: 10px; margin-bottom: 10px; border: 1px solid #eee; padding: 10px; background: #fafafa;",
          shiny::numericInput("lfcCutoff", "logFC cutoff (default is 1)", value = 1, min = 0, step = 0.1),
          shiny::numericInput("pCutoff", "Adjusted Pvalue (default is 0.05)", value = 0.05, min = 0, max = 1, step = 0.01),
          shiny::selectInput("pAdjustMethod", "Method for padj (default is fdr)",
            choices = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"),
            selected = "fdr"
          ),
          shiny::textInput("palette", "Two Colorhexa codes (default is #4f8832 and #f79c18)", value = "#4f8832,#f79c18"),
          shiny::textInput("rowNamesOfCounts", "The name of the rows in counts.txt", value = "external_gene_name"),
          shiny::checkboxInput("extra_plots", "Build QC / exploration plots", value = TRUE),
          shiny::checkboxInput("report", "Write an HTML report", value = TRUE)
        )
      ),
      shiny::tags$hr(),
      shiny::verbatimTextOutput("status"),
      shiny::verbatimTextOutput("output_path"),
      shiny::uiOutput("download_ui")
    ),

    shiny::mainPanel(
      shiny::uiOutput("summary_cards"),
      shiny::tabsetPanel(
        id = "plot_tabs",

        shiny::tabPanel(
          "Gene of interest",
          shiny::tabsetPanel(
            shiny::tabPanel("Boxplot",       shiny::plotOutput("boxplot", height = "520px")),
            shiny::tabPanel("Barplot",       shiny::plotOutput("barplot", height = "520px")),
            shiny::tabPanel("Countplot",     shiny::plotOutput("countplot", height = "520px")),
            shiny::tabPanel("Fold-change rank", shiny::plotOutput("rank_plot", height = "520px")),
            shiny::tabPanel("Co-expression", shiny::plotOutput("correlated_plot", height = "600px")),
            shiny::tabPanel("Diseases",      shiny::plotOutput("disease_plot", height = "520px"))
          ),
          shiny::tags$hr(),
          shiny::tags$h5("Annotation"),
          shiny::tableOutput("goi_des"),
          shiny::tags$h5("Test result"),
          shiny::tableOutput("goi_test")
        ),

        shiny::tabPanel(
          "Differential expression",
          shiny::tabsetPanel(
            shiny::tabPanel("Volcano (labelled)", shiny::plotOutput("volcano_labelled", height = "620px")),
            shiny::tabPanel("Volcano (GOI)",      shiny::plotOutput("volcano_plot", height = "560px")),
            shiny::tabPanel("MA plot",            shiny::plotOutput("ma_plot", height = "560px")),
            shiny::tabPanel("Heatmap",            shiny::plotOutput("deg_heatmap", height = "800px")),
            shiny::tabPanel(
              "DEG table",
              shiny::tags$br(),
              DT::dataTableOutput("degs_table")
            )
          )
        ),

        shiny::tabPanel(
          "Quality control",
          shiny::tabsetPanel(
            shiny::tabPanel("PCA",               shiny::plotOutput("pca_plot", height = "560px")),
            shiny::tabPanel("Sample distances",  shiny::plotOutput("distance_plot", height = "560px")),
            shiny::tabPanel("Library sizes",     shiny::plotOutput("libsize_plot", height = "480px")),
            shiny::tabPanel("p-value histogram", shiny::plotOutput("pvalue_plot", height = "480px")),
            shiny::tabPanel("Dispersion",        shiny::plotOutput("dispersion_plot", height = "560px"))
          )
        ),

        shiny::tabPanel(
          "Ask AI",
          shiny::tags$br(),
          shiny::fluidRow(
            shiny::column(
              4,
              shiny::selectInput(
                "ai_provider", "Provider",
                choices = c("Anthropic (Claude)" = "anthropic", "OpenAI" = "openai"),
                selected = "anthropic"
              ),
              shiny::passwordInput(
                "ai_key", "API key",
                placeholder = "leave empty to use the environment variable"
              ),
              shiny::textInput(
                "ai_model", "Model (optional)",
                placeholder = "provider default"
              ),
              shiny::tags$p(
                class = "hint",
                "Your question and a text summary of the run are sent to the",
                "provider you pick. The count matrix and the results table stay",
                "on this machine. Answers are a reading aid — check anything",
                "you plan to publish against the tables the pipeline wrote."
              )
            ),
            shiny::column(
              8,
              shiny::uiOutput("ai_chat"),
              shiny::textAreaInput(
                "ai_query", NULL,
                placeholder = "Ask something about this run…",
                width = "100%", height = "90px"
              ),
              shiny::div(
                shiny::actionButton("ai_ask", "Ask", class = "btn-primary"),
                shiny::actionButton("ai_clear", "Clear chat")
              ),
              shiny::tags$br(),
              shiny::tags$div(
                class = "hint",
                "Try: ",
                shiny::actionLink("ai_ex1", "interpret this run"), " · ",
                shiny::actionLink("ai_ex2", "is my gene a strong hit?"), " · ",
                shiny::actionLink("ai_ex3", "what do the top genes have in common?"), " · ",
                shiny::actionLink("ai_ex4", "write a figure legend")
              )
            )
          )
        ),

        shiny::tabPanel(
          "Run summary",
          shiny::tags$br(),
          shiny::verbatimTextOutput("run_summary")
        )
      )
    )
  )
)
