# goiExplorer 0.2.0

## New plots

- Ten new figures, each also available as a standalone function so they can be
  redrawn from a finished run: `plot_pca()`, `plot_sample_distances()`,
  `plot_deg_heatmap()`, `plot_pvalue_histogram()`, `plot_dispersion()`,
  `plot_goi_rank()`, `plot_top_correlated()`, `plot_library_sizes()`,
  `plot_disease_associations()` and `plot_volcano_top()`.
- `build_extra_plots()` builds the whole set at once; `pipeline()` calls it
  automatically unless `extra_plots = FALSE`. A plot that cannot be drawn warns
  and is skipped instead of failing the run.
- `theme_goi()` gives every figure the same look.
- The variance-stabilised counts the pipeline already computed are now actually
  used, by the PCA, the distance heatmap, the DEG heatmap and the co-expression
  plot.
- Plots are saved on a white background, so they stay readable when dropped
  into a dark-themed document.

## HTML report

- `write_report()` writes a single self-contained HTML file with the headline
  numbers, every figure (embedded, not linked) and the annotation tables.
  Produced automatically by `run_pipeline()`; pass `report = FALSE` to skip it.

## AI assistant

- `ai_agent()` rewritten. It now sends a summary of the actual run as context,
  supports Anthropic (Claude) as well as OpenAI, holds multi-turn
  conversations, and talks to the APIs over `httr` so no provider package is
  needed. Only the summary is transmitted; the count matrix and results table
  never leave the machine.
- `summarise_run()` produces that summary as plain text, which is useful to
  read on its own.
- `ai_interpret()` is a ready-made "explain this run" prompt.
- Fixed: `R/ai_agent.R` ran `stop()` at the top level, so installing the
  package failed unless the `openai` package happened to be present.

## Shiny app

- Added the **Ask AI** tab. The server already handled `ai_ask`, `ai_query` and
  `openai_api_key`, but the UI never defined them, so the assistant could not
  be reached at all. There is now a provider selector, an API key field, a chat
  transcript and example questions.
- All the new plots have tabs, grouped into *Gene of interest*, *Differential
  expression* and *Quality control*.
- Headline numbers (genes tested, significant, up, down, the GOI's fold change
  and adjusted p-value) are shown above the tabs.
- Download buttons for the whole output directory as a `.zip` and for the HTML
  report.
- A *Run summary* tab showing the text digest.
- Fixed a duplicate `output_path` output that meant the chosen output directory
  was only ever rendered in one of the two places it was declared.
- The API key is passed straight to `ai_agent()` instead of being written into
  the process environment.

## Pipeline

- The disease, KEGG, Open Targets and Manhattan steps are wrapped
  individually: one unreachable service now costs you that section, not the
  whole run and every result computed before it.
- The KEGG step restores the working directory even when it fails part way.
- Manhattan plots are drawn whenever colocalising studies exist. Previously the
  code only handled exactly one or exactly two studies and silently drew
  nothing for three or more.
- `run_pipeline()` reported where the output went *after* returning, so the
  message was never shown.
- The p-value histogram uses the p-values as DESeq2 reported them, before the
  `NA` values are replaced with 1.
- `res_output$plotMAPath` records where the MA plot was written, like the other
  plots do.

## Packaging

- `DESCRIPTION` had the placeholder title and description from `usethis`.
- Declared the packages the app actually uses (`DT`, `fs`) and the ones the AI
  client needs (`httr`, `jsonlite`); dropped `openai`, which is no longer used.
- Added `SummarizedExperiment`, `S4Vectors`, `grDevices`, `stats` and `utils`
  to `Imports`.
- Tests for the plots, the summary, the report and the app's chat logic, none
  of which need the network.

# goiExplorer 0.0.0.9000

- Fixed bar‐plot summary to import rstatix::mean_sd correctly.
- Ensured `baseMean` is always present by recomputing it from normalized counts.
