# goiExplorer
Exploring a gene of interest

An R package for exploring Genes Of Interest (GOIs) via differential‐expression pipelines, pathway analysis, and a Shiny interface. See the tutorial on [YouTube](https://youtu.be/pq5Wg64rjAU)

Input files:
- Count Table (counts.txt)
- Salmon (quant directory)

## Installation

```r
# from GitHub
devtools::install_github("mortezasabri/goiExplorer")
```

## Shiny App (User-friendly)

```r
goiExplorer::run_app()
```

## Running on R console 

With more freedom to change the arguments

```r
library(goiExplorer)

res_output <- run_pipeline(
  input         = "path/to/counts_dir",
  dataType      = "counts",
  goi           = "TP53",
  parent_outdir = "path/to/output"
)
head(res_output$counts)
p <- res_output$Volcanoplot  # volcano plot
p
```

Other plots:
- Barplot
- Boxplot
- Countplot
- plotMA



## Contributing

Pull requests and issues welcome:
https://github.com/mortezasabri/goiExplorer
