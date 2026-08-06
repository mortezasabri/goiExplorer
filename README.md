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

## Plots

Each plot is also a standalone function, so you can redraw any of them from a
finished run without re-running the analysis.

**About the gene of interest**

| Slot | Function | What it shows |
| --- | --- | --- |
| `Barplot` | | group means ± SD with the adjusted p-value |
| `Boxplot` | | per-group distribution with every sample drawn |
| `Countplot` | | normalised counts per sample, log scale |
| `RankPlot` | `plot_goi_rank()` | where the GOI sits among all genes by fold change, with its rank and percentile |
| `CorrelatedGenesPlot` | `plot_top_correlated()` | the genes whose expression tracks the GOI across samples |
| `DiseasePlot` | `plot_disease_associations()` | associated diseases and their scores |

**About the experiment**

| Slot | Function | What it shows |
| --- | --- | --- |
| `Volcanoplot` | | volcano plot highlighting the GOI |
| `VolcanoLabelled` | `plot_volcano_top()` | volcano plot naming the strongest genes, with up/down counts |
| `plotMA` | | MA plot |
| `DEGHeatmap` | `plot_deg_heatmap()` | top DE genes as row z-scores, clustered, GOI outlined |
| `KEGGpaths` | | KEGG pathway diagrams containing the GOI |

**Quality control** — worth a look before trusting any of the above

| Slot | Function | What it shows |
| --- | --- | --- |
| `PCAplot` | `plot_pca()` | sample PCA with the variance explained per axis |
| `SampleDistancePlot` | `plot_sample_distances()` | sample-to-sample distances, clustered |
| `LibrarySizePlot` | `plot_library_sizes()` | library sizes and DESeq2 size factors |
| `PvalueHistogram` | `plot_pvalue_histogram()` | p-value distribution — flat with a peak at zero is healthy |
| `DispersionPlot` | `plot_dispersion()` | dispersion estimates against the fitted trend |

```r
# redraw anything from a finished run
plot_pca(res_output, ntop = 2000)
plot_deg_heatmap(res_output, n = 80)
plot_top_correlated(res_output, n = 40)

# or rebuild the whole set with different settings
res_output <- build_extra_plots(res_output, n_heatmap = 60, pCutoff = 0.01)
```

## HTML report

Every run writes a single self-contained HTML file — numbers, figures and
tables in one document you can email or archive. Figures are embedded, so
there is nothing else to keep track of.

```r
write_report(res_output)                              # to the run's output dir
write_report(res_output, file = "~/report.html")      # somewhere specific
write_report(res_output, ai = TRUE)                   # with an LLM interpretation
```

Pass `report = FALSE` to `run_pipeline()` to skip it.

## Ask about your results

`ai_agent()` sends your question, together with a text summary of the run, to
Claude or to OpenAI and returns the answer. It can talk about *your* numbers
rather than guessing.

```r
Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")   # or OPENAI_API_KEY

ai_agent("Is TP53 among the strongest changes here, or just significant?", res_output)
ai_interpret(res_output)                       # a ready-made "explain this run" prompt

cat(summarise_run(res_output))                 # exactly what gets sent as context
```

The same thing is available in the app under the **Ask AI** tab. Only the
summary is transmitted — the count matrix and the results table stay on your
machine. Treat the answers as a reading aid and verify anything you plan to
publish against the tables the pipeline wrote to disk.

## Contributing

Pull requests and issues welcome:
https://github.com/mortezasabri/goiExplorer
