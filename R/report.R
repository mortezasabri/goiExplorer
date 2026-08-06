## ---------------------------------------------------------------------------
## A single-file HTML report for a finished run. No rmarkdown / pandoc needed:
## the figures are embedded as data URIs when knitr is available, otherwise
## they are linked relatively (which works because the report is written next
## to the PNGs).
## ---------------------------------------------------------------------------

.html_escape <- function(x) {
  x <- base::as.character(x)
  x[base::is.na(x)] <- ""
  x <- base::gsub("&", "&amp;", x, fixed = TRUE)
  x <- base::gsub("<", "&lt;", x, fixed = TRUE)
  x <- base::gsub(">", "&gt;", x, fixed = TRUE)
  x <- base::gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

# <img> tag for a PNG: embedded when possible, linked otherwise.
.html_img <- function(path, alt) {
  if (base::is.null(path) || !base::file.exists(path)) {
    return("")
  }
  src <- NULL
  if (base::requireNamespace("knitr", quietly = TRUE)) {
    src <- base::tryCatch(knitr::image_uri(path), error = function(e) NULL)
  }
  if (base::is.null(src)) src <- base::basename(path)
  base::paste0(
    "<img src=\"", src, "\" alt=\"", .html_escape(alt), "\" loading=\"lazy\">"
  )
}

# Render a data.frame as an HTML table (first `n` rows, selected columns).
.html_table <- function(df, n = 25, cols = NULL, digits = 3) {
  if (base::is.null(df) || base::nrow(df) == 0) {
    return("<p class=\"muted\">No rows.</p>")
  }
  if (!base::is.null(cols)) {
    cols <- base::intersect(cols, base::names(df))
    if (base::length(cols) > 0) df <- df[, cols, drop = FALSE]
  }
  keep_names <- base::rownames(df)
  df <- df[base::seq_len(base::min(n, base::nrow(df))), , drop = FALSE]
  for (j in base::seq_along(df)) {
    if (base::is.numeric(df[[j]])) df[[j]] <- base::signif(df[[j]], digits)
  }
  has_names <- !base::is.null(keep_names) && !base::identical(keep_names, base::as.character(base::seq_along(keep_names)))

  head_cells <- base::paste0("<th>", .html_escape(base::names(df)), "</th>", collapse = "")
  if (has_names) head_cells <- base::paste0("<th>gene</th>", head_cells)

  rows <- base::vapply(base::seq_len(base::nrow(df)), function(i) {
    cells <- base::paste0(
      "<td>",
      .html_escape(base::vapply(df, function(col) base::as.character(col[i]), base::character(1))),
      "</td>",
      collapse = ""
    )
    if (has_names) {
      cells <- base::paste0("<td class=\"gene\">", .html_escape(base::rownames(df)[i]), "</td>", cells)
    }
    base::paste0("<tr>", cells, "</tr>")
  }, base::character(1))

  base::paste0(
    "<div class=\"tablewrap\"><table><thead><tr>", head_cells, "</tr></thead><tbody>",
    base::paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

.html_stat <- function(label, value, note = "") {
  base::paste0(
    "<div class=\"stat\"><div class=\"stat-value\">", .html_escape(value),
    "</div><div class=\"stat-label\">", .html_escape(label), "</div>",
    if (base::nzchar(note)) base::paste0("<div class=\"stat-note\">", .html_escape(note), "</div>") else "",
    "</div>"
  )
}

.report_css <- function() {
  "
:root {
  --ink:#1f2933; --muted:#7b8794; --line:#e4e7eb; --bg:#ffffff;
  --accent:#4f8832; --accent2:#f79c18; --card:#f8f9fa;
}
* { box-sizing:border-box; }
body { margin:0; padding:0 1.25rem 4rem; background:var(--bg); color:var(--ink);
  font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
  line-height:1.55; }
.wrap { max-width:1000px; margin:0 auto; }
header { border-bottom:3px solid var(--accent); padding:2.5rem 0 1.25rem; margin-bottom:2rem; }
h1 { margin:0 0 .25rem; font-size:1.9rem; letter-spacing:-.02em; }
h2 { margin:2.5rem 0 .75rem; font-size:1.25rem; padding-bottom:.4rem;
  border-bottom:1px solid var(--line); }
h3 { margin:1.5rem 0 .5rem; font-size:1rem; color:var(--muted);
  text-transform:uppercase; letter-spacing:.06em; }
.sub { color:var(--muted); font-size:.95rem; }
.muted { color:var(--muted); }
.stats { display:flex; flex-wrap:wrap; gap:.75rem; margin:1.25rem 0; }
.stat { flex:1 1 150px; background:var(--card); border:1px solid var(--line);
  border-radius:10px; padding:.9rem 1rem; }
.stat-value { font-size:1.5rem; font-weight:650; letter-spacing:-.02em; }
.stat-label { font-size:.8rem; color:var(--muted); text-transform:uppercase;
  letter-spacing:.05em; margin-top:.15rem; }
.stat-note { font-size:.8rem; color:var(--muted); margin-top:.35rem; }
figure { margin:1.5rem 0; }
figure img { width:100%; height:auto; border:1px solid var(--line); border-radius:8px; }
figcaption { font-size:.87rem; color:var(--muted); margin-top:.5rem; }
.tablewrap { overflow-x:auto; border:1px solid var(--line); border-radius:8px; }
table { border-collapse:collapse; width:100%; font-size:.85rem; }
th, td { padding:.45rem .7rem; text-align:left; white-space:nowrap;
  border-bottom:1px solid var(--line); }
th { background:var(--card); font-weight:600; position:sticky; top:0; }
tbody tr:last-child td { border-bottom:none; }
td.gene { font-weight:600; font-family:ui-monospace,SFMono-Regular,Menlo,monospace; }
ul.plain { padding-left:1.1rem; }
pre.ai { background:var(--card); border:1px solid var(--line); border-left:3px solid var(--accent2);
  border-radius:8px; padding:1rem; white-space:pre-wrap; font-size:.9rem;
  font-family:inherit; }
footer { margin-top:3rem; padding-top:1rem; border-top:1px solid var(--line);
  color:var(--muted); font-size:.85rem; }
@media (prefers-color-scheme: dark) {
  :root { --ink:#e8eaed; --muted:#9aa5b1; --line:#333a42; --bg:#16191d; --card:#1e2227; }
  figure img { background:#fff; }
}
@media print { body { padding:0; } figure { page-break-inside:avoid; } }
"
}

#' Write a single-file HTML report for a run
#'
#' Collects the numbers, the figures and the annotation tables from a finished
#' run into one self-contained HTML file you can email, attach to a lab
#' notebook, or open months later without re-running anything. Figures are
#' embedded directly in the file when the knitr package is available.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param file Character. Output path. Defaults to `goiExplorer_report_<GOI>.html`
#'   inside the run's output directory.
#' @param n_genes Integer. Number of DE genes to include in the table.
#' @param ai Logical. Also ask the LLM for an interpretation and include it.
#'   Requires an API key; see [ai_agent()]. Failures are reported in the report
#'   rather than raised.
#' @param ... Passed to [ai_interpret()] when `ai = TRUE`.
#' @return The path to the written file, invisibly.
#' @examples
#' \dontrun{
#' res <- run_pipeline(input = "counts_dir", dataType = "counts", goi = "TP53")
#' write_report(res)
#' write_report(res, ai = TRUE) # with an LLM interpretation
#' }
#' @export
write_report <- function(res_output, file = NULL, n_genes = 50, ai = FALSE, ...) {
  stopifnot(base::is.list(res_output))
  goi <- if (base::is.null(res_output$goi)) "GOI" else res_output$goi

  if (base::is.null(file)) {
    outdir <- res_output$directory
    if (base::is.null(outdir) || !base::nzchar(outdir)) outdir <- base::tempdir()
    file <- base::file.path(outdir, base::paste0("goiExplorer_report_", goi, ".html"))
  }

  res.df <- res_output$res.df
  degs <- res_output$degs
  html <- base::character(0)

  ## --- header + headline numbers -----------------------------------------
  n_tested <- if (!base::is.null(res.df)) base::nrow(res.df) else NA_integer_
  n_deg <- if (!base::is.null(degs)) base::nrow(degs) else NA_integer_
  n_up <- if (!base::is.null(degs)) {
    base::sum(degs$log2FoldChange > 0, na.rm = TRUE)
  } else {
    NA_integer_
  }
  n_dn <- if (!base::is.null(degs)) {
    base::sum(degs$log2FoldChange < 0, na.rm = TRUE)
  } else {
    NA_integer_
  }

  goi_row <- if (!base::is.null(res.df) && goi %in% base::rownames(res.df)) {
    res.df[goi, , drop = FALSE]
  } else {
    NULL
  }

  stats <- c(
    .html_stat("genes tested", if (base::is.na(n_tested)) "-" else base::format(n_tested, big.mark = ",")),
    .html_stat("significant", if (base::is.na(n_deg)) "-" else base::format(n_deg, big.mark = ",")),
    .html_stat("up-regulated", if (base::is.na(n_up)) "-" else n_up),
    .html_stat("down-regulated", if (base::is.na(n_dn)) "-" else n_dn)
  )
  if (!base::is.null(goi_row)) {
    stats <- c(
      stats,
      .html_stat(
        base::paste0(goi, " log2FC"),
        base::round(goi_row$log2FoldChange[1], 2),
        base::paste0("padj = ", base::signif(goi_row$padj[1], 3))
      )
    )
  }

  html <- c(
    html,
    "<div class=\"wrap\">",
    "<header>",
    base::paste0("<h1>", .html_escape(goi), " &mdash; RNA-Seq report</h1>"),
    base::paste0(
      "<p class=\"sub\">Generated by goiExplorer on ",
      .html_escape(base::format(base::Sys.time(), "%Y-%m-%d %H:%M")),
      if (!base::is.null(res_output$ensemblSpecies)) {
        base::paste0(" &middot; ", .html_escape(res_output$ensemblSpecies))
      } else {
        ""
      },
      "</p>"
    ),
    "</header>",
    "<div class=\"stats\">", stats, "</div>"
  )

  ## --- gene of interest ---------------------------------------------------
  html <- c(html, "<h2>Gene of interest</h2>")
  if (!base::is.null(goi_row)) {
    html <- c(html, .html_table(
      goi_row, n = 1,
      cols = c(
        "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
        "gene_biotype", "entrezgene_id", "chromosome_name"
      )
    ))
    desc <- goi_row$description
    if (!base::is.null(desc) && !base::is.na(desc[1])) {
      html <- c(html, base::paste0("<p class=\"muted\">", .html_escape(desc[1]), "</p>"))
    }
  } else {
    html <- c(html, base::paste0(
      "<p class=\"muted\">", .html_escape(goi),
      " was not found in the results table.</p>"
    ))
  }

  ## --- figures ------------------------------------------------------------
  figures <- base::list(
    list("BoxplotPath", "Expression of the gene of interest per group."),
    list("BarplotPath", "Group means with standard deviation and the adjusted p-value."),
    list("CountplotPath", "Normalised counts per sample (log scale)."),
    list("VolcanoLabelledPath", "Volcano plot with the strongest genes labelled."),
    list("VolcanoplotPath", "Volcano plot highlighting the gene of interest."),
    list("RankPlotPath", "Where the gene of interest ranks by fold change."),
    list("DEGHeatmapPath", "Top differentially expressed genes, row z-scores."),
    list("CorrelatedGenesPlotPath", "Genes co-expressed with the gene of interest."),
    list("PCAplotPath", "Principal component analysis of the samples."),
    list("SampleDistancePlotPath", "Sample-to-sample distances."),
    list("LibrarySizePlotPath", "Library sizes and size factors."),
    list("PvalueHistogramPath", "Distribution of raw p-values."),
    list("DispersionPlotPath", "Dispersion estimates."),
    list("plotMAPath", "MA plot."),
    list("DiseasePlotPath", "Diseases associated with the gene of interest.")
  )
  fig_html <- base::character(0)
  for (f in figures) {
    p <- res_output[[f[[1]]]]
    if (base::is.null(p) || !base::file.exists(p)) next
    fig_html <- c(fig_html, base::paste0(
      "<figure>", .html_img(p, f[[2]]),
      "<figcaption>", .html_escape(f[[2]]), "</figcaption></figure>"
    ))
  }
  if (base::length(fig_html) > 0) {
    html <- c(html, "<h2>Figures</h2>", fig_html)
  }

  ## --- DE table -----------------------------------------------------------
  if (!base::is.null(degs) && base::nrow(degs) > 0) {
    html <- c(
      html,
      base::paste0("<h2>Differentially expressed genes</h2>"),
      base::paste0(
        "<p class=\"sub\">Showing the top ",
        base::min(n_genes, base::nrow(degs)), " of ", base::nrow(degs),
        " genes that passed the cutoffs, ranked by adjusted p-value.</p>"
      ),
      .html_table(
        degs[base::order(degs$padj), , drop = FALSE],
        n = n_genes,
        cols = c("baseMean", "log2FoldChange", "pvalue", "padj", "gene_biotype")
      )
    )
  }

  ## --- diseases -----------------------------------------------------------
  dis <- res_output$DEGsToDiseases
  if (!base::is.null(dis) && base::nrow(dis) > 0) {
    html <- c(html, "<h2>Associated diseases</h2>", .html_table(dis, n = 20))
  }

  ## --- pathways -----------------------------------------------------------
  paths <- res_output$KEGGpaths
  if (!base::is.null(paths) && base::length(paths) > 0) {
    nm <- base::sub("\\.png$", "", base::basename(base::unlist(paths)))
    html <- c(
      html, "<h2>KEGG pathways</h2>",
      "<ul class=\"plain\">",
      base::paste0("<li>", .html_escape(base::unique(nm)), "</li>"),
      "</ul>"
    )
  }

  ## --- optional AI interpretation ----------------------------------------
  if (base::isTRUE(ai)) {
    txt <- base::tryCatch(
      ai_interpret(res_output, ...),
      error = function(e) base::paste("AI request failed:", base::conditionMessage(e))
    )
    html <- c(
      html, "<h2>Interpretation</h2>",
      "<p class=\"sub\">Generated by a language model from the summary above. Verify before citing.</p>",
      base::paste0("<pre class=\"ai\">", .html_escape(txt), "</pre>")
    )
  }

  ## --- run summary --------------------------------------------------------
  html <- c(
    html,
    "<h2>Run summary</h2>",
    base::paste0("<pre class=\"ai\">", .html_escape(summarise_run(res_output)), "</pre>"),
    "<footer>Produced by the goiExplorer R package.</footer>",
    "</div>"
  )

  doc <- base::paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>goiExplorer &mdash; ", .html_escape(goi), "</title>\n<style>",
    .report_css(), "</style>\n</head>\n<body>\n",
    base::paste(html, collapse = "\n"),
    "\n</body>\n</html>\n"
  )

  base::writeLines(doc, file, useBytes = TRUE)
  base::message("Report written to: ", file)
  base::invisible(file)
}
