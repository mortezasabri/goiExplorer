# goiExplorer
Exploring a gene of interest

An R package for exploring Genes Of Interest (GOIs) via differential‐expression pipelines, pathway analysis, and a Shiny interface with an AI agent.

## Installation

```r
# from GitHub
# devtools::install_github("mortezasabri/goiExplorer")
```

## Quickstart

```r
library(goiExplorer)

res_output <- run_pipeline(
  input         = "path/to/counts_dir",
  dataType      = "counts",
  goi           = "TP53",
  parent_outdir = "path/to/output"
)
head(res_output$counts)
res_output$Volcanoplot  # volcano plot
```

## Shiny App

```r
goiExplorer::run_app()
```

## Contributing

Pull requests and issues welcome:
https://github.com/mortezasabri/goiExplorer