## ---------------------------------------------------------------------------
## Additional exploratory / QC plots for a goiExplorer run.
##
## Every function here takes the `res_output` list produced by `pipeline()` /
## `run_pipeline()` and returns a ggplot object (or `NULL` when the pieces it
## needs are missing). Nothing here ever stops the pipeline: a plot that cannot
## be drawn simply warns and returns `NULL`.
## ---------------------------------------------------------------------------

# Row-wise variance without pulling in matrixStats.
.row_vars <- function(x) {
  n <- ncol(x)
  if (is.null(n) || n < 2) {
    return(rep(NA_real_, nrow(x)))
  }
  base::rowSums((x - base::rowMeans(x))^2) / (n - 1)
}

# Variance-stabilised expression matrix (genes x samples) from res_output.
.vsd_matrix <- function(res_output) {
  vsd <- res_output$vsd
  if (base::is.null(vsd)) {
    return(NULL)
  }
  if (base::is.matrix(vsd)) {
    return(vsd)
  }
  if (base::requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    return(base::as.matrix(SummarizedExperiment::assay(vsd)))
  }
  NULL
}

# colData of the DESeqDataSet, falling back to the per-gene count table.
.coldata <- function(res_output) {
  dds <- res_output$dds
  if (!base::is.null(dds) &&
    base::requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    cd <- base::tryCatch(
      base::as.data.frame(SummarizedExperiment::colData(dds)),
      error = function(e) NULL
    )
    if (!base::is.null(cd) && base::nrow(cd) > 0) {
      return(cd)
    }
  }
  gc <- res_output$geneCounts
  if (!base::is.null(gc)) {
    return(gc[, base::setdiff(base::names(gc), "count"), drop = FALSE])
  }
  NULL
}

# mcols() of a DESeqDataSet without hard-depending on S4Vectors.
.mcols <- function(x) {
  if (base::requireNamespace("S4Vectors", quietly = TRUE)) {
    return(base::tryCatch(S4Vectors::mcols(x), error = function(e) NULL))
  }
  if (base::requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    return(base::tryCatch(SummarizedExperiment::mcols(x), error = function(e) NULL))
  }
  NULL
}

# hclust ordering that degrades gracefully for tiny matrices.
.cluster_order <- function(m, by = c("row", "col")) {
  by <- base::match.arg(by)
  if (by == "col") m <- base::t(m)
  n <- base::nrow(m)
  if (base::is.null(n) || n < 3) {
    return(base::seq_len(base::max(n, 0)))
  }
  base::tryCatch(
    stats::hclust(stats::dist(m))$order,
    error = function(e) base::seq_len(n)
  )
}

#' A shared minimal theme for goiExplorer plots
#'
#' Thin wrapper around [ggplot2::theme_minimal()] that centres and bolds the
#' title, so every plot the package produces looks like it belongs to the same
#' report.
#'
#' @param base_size Numeric. Base font size passed to [ggplot2::theme_minimal()].
#' @return A ggplot2 theme object.
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   ggtitle("Demo") +
#'   theme_goi()
#' }
#' @export
theme_goi <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey30"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Save a ggplot and return the path (NULL when outdir is not supplied).
.save_plot <- function(g, outdir, name, width = 7, height = 5) {
  if (base::is.null(g) || base::is.null(outdir) || !base::nzchar(outdir)) {
    return(NULL)
  }
  p <- base::file.path(outdir, name)
  base::tryCatch(
    {
      # bg = "white": a transparent PNG turns into an unreadable black-on-black
      # figure the moment it lands in a dark-themed report or slide
      ggplot2::ggsave(p, g, dpi = 300, width = width, height = height, bg = "white")
      p
    },
    error = function(e) {
      base::warning("Could not save ", name, ": ", e$message)
      NULL
    }
  )
}


#' PCA of the samples
#'
#' Principal component analysis on the variance-stabilised counts, using the
#' `ntop` most variable genes. The percentage of variance explained is printed
#' on each axis, which makes it easy to see whether the experimental groups are
#' actually the main source of variation in the data.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param ntop Integer. Number of most variable genes to use (default 500).
#' @param intgroup Character. Column of `colData` used for colouring.
#' @param palette Character vector of fill colours, one per group.
#' @param label_samples Logical. Draw sample names next to the points.
#' @return A ggplot object, or `NULL` if no variance-stabilised data is available.
#' @examples
#' \dontrun{
#' res <- run_pipeline(input = "counts_dir", dataType = "counts", goi = "TP53")
#' plot_pca(res)
#' }
#' @importFrom ggplot2 ggplot aes labs .data
#' @export
plot_pca <- function(res_output,
                     ntop = 500,
                     intgroup = "condition",
                     palette = c("#4f8832", "#f79c18"),
                     label_samples = TRUE) {
  mat <- .vsd_matrix(res_output)
  if (base::is.null(mat) || base::ncol(mat) < 2) {
    base::warning("plot_pca(): no variance-stabilised matrix available.")
    return(NULL)
  }

  rv <- .row_vars(mat)
  sel <- base::order(rv, decreasing = TRUE)[base::seq_len(base::min(ntop, base::nrow(mat)))]
  pca <- stats::prcomp(base::t(mat[sel, , drop = FALSE]))
  pct <- base::round(100 * pca$sdev^2 / base::sum(pca$sdev^2), 1)

  df <- base::data.frame(
    PC1 = pca$x[, 1],
    PC2 = if (base::ncol(pca$x) >= 2) pca$x[, 2] else 0,
    sample = base::colnames(mat),
    stringsAsFactors = FALSE
  )

  cd <- .coldata(res_output)
  if (!base::is.null(cd) && intgroup %in% base::names(cd)) {
    df$group <- base::as.factor(cd[[intgroup]][base::match(df$sample, base::rownames(cd))])
  } else {
    df$group <- base::factor("all")
  }

  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$group)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.9) +
    ggplot2::labs(
      title = "Sample PCA",
      subtitle = base::paste0("top ", base::length(sel), " most variable genes"),
      x = base::paste0("PC1: ", pct[1], "% variance"),
      y = base::paste0("PC2: ", if (base::length(pct) >= 2) pct[2] else 0, "% variance"),
      colour = intgroup
    ) +
    theme_goi()

  if (base::length(base::levels(df$group)) <= base::length(palette)) {
    g <- g + ggplot2::scale_colour_manual(values = palette)
  }
  if (label_samples) {
    g <- g + ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$sample),
      size = 3, show.legend = FALSE, max.overlaps = Inf
    )
  }
  g
}


#' Sample-to-sample distance heatmap
#'
#' Euclidean distances between samples on the variance-stabilised scale,
#' reordered by hierarchical clustering. Replicates of the same group should sit
#' next to each other; a sample that clusters with the wrong group is usually a
#' swap or an outlier worth checking before trusting the DE results.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param ntop Integer. Number of most variable genes used for the distance.
#'   Use `Inf` for all genes.
#' @return A ggplot object, or `NULL` when no expression matrix is available.
#' @examples
#' \dontrun{
#' plot_sample_distances(res)
#' }
#' @export
plot_sample_distances <- function(res_output, ntop = Inf) {
  mat <- .vsd_matrix(res_output)
  if (base::is.null(mat) || base::ncol(mat) < 2) {
    base::warning("plot_sample_distances(): no variance-stabilised matrix available.")
    return(NULL)
  }
  if (base::is.finite(ntop) && ntop < base::nrow(mat)) {
    rv <- .row_vars(mat)
    mat <- mat[base::order(rv, decreasing = TRUE)[base::seq_len(ntop)], , drop = FALSE]
  }

  dm <- base::as.matrix(stats::dist(base::t(mat)))
  ord <- .cluster_order(dm, by = "row")
  dm <- dm[ord, ord, drop = FALSE]
  labs <- base::rownames(dm)

  long <- base::data.frame(
    s1 = base::factor(base::rep(labs, times = base::length(labs)), levels = labs),
    s2 = base::factor(base::rep(labs, each = base::length(labs)), levels = base::rev(labs)),
    distance = base::as.vector(dm),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(long, ggplot2::aes(x = .data$s1, y = .data$s2, fill = .data$distance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "#2c3e50", high = "white", name = "distance") +
    ggplot2::labs(
      title = "Sample-to-sample distances",
      subtitle = "Euclidean distance on variance-stabilised counts",
      x = NULL, y = NULL
    ) +
    theme_goi() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )
}


#' Heatmap of the top differentially expressed genes
#'
#' Row z-scores of the variance-stabilised counts for the `n` most significant
#' genes, clustered on both axes. The gene of interest is always included and
#' its label is highlighted, so you can see immediately whether it behaves like
#' the rest of the signature.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param n Integer. Number of top genes (ranked by adjusted p-value).
#' @param goi Character. Gene to highlight; defaults to the pipeline's GOI.
#' @param highlight_colour Colour used for the GOI label.
#' @return A ggplot object, or `NULL` when the required pieces are missing.
#' @examples
#' \dontrun{
#' plot_deg_heatmap(res, n = 30)
#' }
#' @export
plot_deg_heatmap <- function(res_output,
                             n = 40,
                             goi = res_output$goi,
                             highlight_colour = "#c0392b") {
  mat <- .vsd_matrix(res_output)
  res.df <- res_output$res.df
  if (base::is.null(mat) || base::is.null(res.df) || base::nrow(res.df) == 0) {
    base::warning("plot_deg_heatmap(): need both `vsd` and `res.df`.")
    return(NULL)
  }

  ranked <- res.df[base::order(res.df$padj, -base::abs(res.df$log2FoldChange)), , drop = FALSE]
  genes <- base::rownames(ranked)[base::seq_len(base::min(n, base::nrow(ranked)))]
  if (!base::is.null(goi) && goi %in% base::rownames(mat)) {
    genes <- base::unique(c(genes, goi))
  }
  genes <- base::intersect(genes, base::rownames(mat))
  if (base::length(genes) < 2) {
    base::warning("plot_deg_heatmap(): fewer than two genes to plot.")
    return(NULL)
  }

  z <- base::t(base::scale(base::t(mat[genes, , drop = FALSE])))
  z[!base::is.finite(z)] <- 0

  gord <- .cluster_order(z, by = "row")
  sord <- .cluster_order(z, by = "col")
  z <- z[gord, sord, drop = FALSE]

  long <- base::data.frame(
    gene = base::factor(base::rep(base::rownames(z), times = base::ncol(z)),
      levels = base::rownames(z)
    ),
    sample = base::factor(base::rep(base::colnames(z), each = base::nrow(z)),
      levels = base::colnames(z)
    ),
    z = base::as.vector(z),
    stringsAsFactors = FALSE
  )

  goi_row <- base::match(goi, base::rownames(z))
  sub <- "row z-scores of variance-stabilised counts"

  g <- ggplot2::ggplot(long, ggplot2::aes(x = .data$sample, y = .data$gene, fill = .data$z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, name = "z-score"
    ) +
    theme_goi() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 7),
      panel.grid = ggplot2::element_blank()
    )

  # outline the GOI's row rather than recolouring its axis label: ggplot2 only
  # tolerates vectorised element_text() by accident, and warns about it
  if (!base::is.na(goi_row)) {
    g <- g + ggplot2::annotate("rect",
      xmin = 0.5, xmax = base::ncol(z) + 0.5,
      ymin = goi_row - 0.5, ymax = goi_row + 0.5,
      fill = NA, colour = highlight_colour, linewidth = 0.7
    )
    sub <- base::paste0(sub, "  |  ", goi, " outlined")
  }

  g + ggplot2::labs(
    title = base::paste0("Top ", base::length(genes), " differentially expressed genes"),
    subtitle = sub,
    x = NULL, y = NULL
  )
}


#' Histogram of raw p-values
#'
#' A sanity check on the differential expression test. A well-behaved
#' experiment gives a flat histogram with a peak near zero; a hill in the middle
#' or a spike at one usually points at a mis-specified design or heavy
#' filtering.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param bins Integer. Number of histogram bins.
#' @param fill Bar fill colour.
#' @return A ggplot object, or `NULL` when no p-values are available.
#' @examples
#' \dontrun{
#' plot_pvalue_histogram(res)
#' }
#' @export
plot_pvalue_histogram <- function(res_output, bins = 50, fill = "#4f8832") {
  pv <- res_output$pvalue_raw
  if (base::is.null(pv) && !base::is.null(res_output$res.df)) {
    pv <- res_output$res.df$pvalue
  }
  pv <- pv[base::is.finite(pv)]
  if (base::length(pv) == 0) {
    base::warning("plot_pvalue_histogram(): no p-values available.")
    return(NULL)
  }

  df <- base::data.frame(pvalue = pv)
  expected <- base::length(pv) / bins

  ggplot2::ggplot(df, ggplot2::aes(x = .data$pvalue)) +
    ggplot2::geom_histogram(bins = bins, fill = fill, colour = "white", linewidth = 0.2) +
    ggplot2::geom_hline(
      yintercept = expected, linetype = "dashed",
      colour = "grey40", linewidth = 0.4
    ) +
    ggplot2::labs(
      title = "Distribution of raw p-values",
      subtitle = "dashed line = expectation under the null hypothesis",
      x = "p-value", y = "number of genes"
    ) +
    theme_goi()
}


#' Dispersion estimates
#'
#' The ggplot equivalent of [DESeq2::plotDispEsts()]: per-gene estimates, the
#' fitted mean-dispersion trend, and the final shrunken values used by the test.
#' Gene-wise estimates should scatter around the fitted red line.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @return A ggplot object, or `NULL` when the dispersion columns are missing.
#' @examples
#' \dontrun{
#' plot_dispersion(res)
#' }
#' @export
plot_dispersion <- function(res_output) {
  dds <- res_output$dds
  md <- if (base::is.null(dds)) NULL else .mcols(dds)
  if (base::is.null(md) || !all(c("baseMean", "dispGeneEst", "dispFit") %in% base::colnames(md))) {
    base::warning("plot_dispersion(): dispersion estimates not found on `dds`.")
    return(NULL)
  }

  df <- base::data.frame(
    baseMean = base::as.numeric(md$baseMean),
    geneEst  = base::as.numeric(md$dispGeneEst),
    fitted   = base::as.numeric(md$dispFit),
    final    = base::as.numeric(md$dispersion)
  )
  df <- df[base::is.finite(df$baseMean) & df$baseMean > 0, , drop = FALSE]
  if (base::nrow(df) == 0) {
    base::warning("plot_dispersion(): nothing to plot.")
    return(NULL)
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$baseMean)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$geneEst, colour = "gene estimate"),
      size = 0.4, alpha = 0.35, na.rm = TRUE
    ) +
    ggplot2::geom_point(ggplot2::aes(y = .data$final, colour = "final (shrunken)"),
      size = 0.4, alpha = 0.35, na.rm = TRUE
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fitted, colour = "fitted trend"),
      linewidth = 0.8, na.rm = TRUE
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(
      values = c(
        "gene estimate"    = "grey50",
        "final (shrunken)" = "#2166ac",
        "fitted trend"     = "#b2182b"
      ),
      name = NULL
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(
      title = "Dispersion estimates",
      x = "mean of normalised counts", y = "dispersion"
    ) +
    theme_goi() +
    ggplot2::theme(legend.position = "bottom")
}


#' Where the gene of interest ranks among all genes
#'
#' Every gene ordered by log2 fold change, with the GOI marked. The subtitle
#' reports its rank and percentile, which answers the very first question people
#' ask about a GOI: "is this actually one of the strongest changes, or just one
#' that happens to be significant?"
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param goi Character. Gene to highlight; defaults to the pipeline's GOI.
#' @param lfcCutoff Numeric. Fold-change cutoff drawn as guide lines.
#' @param pCutoff Numeric. Adjusted p-value cutoff used to colour significance.
#' @return A ggplot object, or `NULL` when results are missing.
#' @examples
#' \dontrun{
#' plot_goi_rank(res)
#' }
#' @export
plot_goi_rank <- function(res_output,
                          goi = res_output$goi,
                          lfcCutoff = 1,
                          pCutoff = 0.05) {
  res.df <- res_output$res.df
  if (base::is.null(res.df) || base::nrow(res.df) == 0) {
    base::warning("plot_goi_rank(): `res.df` not available.")
    return(NULL)
  }

  df <- res.df[base::is.finite(res.df$log2FoldChange), , drop = FALSE]
  df <- df[base::order(df$log2FoldChange), , drop = FALSE]
  df$rank <- base::seq_len(base::nrow(df))
  df$significant <- base::ifelse(
    !base::is.na(df$padj) & df$padj < pCutoff & base::abs(df$log2FoldChange) > lfcCutoff,
    "significant", "not significant"
  )

  hit <- df[base::rownames(df) %in% goi, , drop = FALSE]
  sub <- if (base::nrow(hit) == 1) {
    base::paste0(
      goi, ": rank ", hit$rank[1], " of ", base::nrow(df),
      " (", base::round(100 * hit$rank[1] / base::nrow(df), 1), "th percentile), log2FC = ",
      base::round(hit$log2FoldChange[1], 2)
    )
  } else {
    base::paste0(goi, " not found in the results table")
  }

  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$log2FoldChange)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$significant), size = 0.4, alpha = 0.6) +
    ggplot2::scale_colour_manual(
      values = c("significant" = "#b2182b", "not significant" = "grey65"),
      name = NULL
    ) +
    ggplot2::geom_hline(
      yintercept = c(-lfcCutoff, 0, lfcCutoff),
      linetype = c("dashed", "solid", "dashed"),
      colour = "grey40", linewidth = 0.3
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(
      title = base::paste0("Fold-change ranking (", goi, ")"),
      subtitle = sub,
      x = "genes ranked by log2 fold change", y = "log2 fold change"
    ) +
    theme_goi() +
    ggplot2::theme(legend.position = "bottom")

  if (base::nrow(hit) == 1) {
    g <- g +
      ggplot2::geom_point(
        data = hit, ggplot2::aes(x = .data$rank, y = .data$log2FoldChange),
        colour = "#2c3e50", size = 3
      ) +
      ggrepel::geom_text_repel(
        data = hit,
        ggplot2::aes(x = .data$rank, y = .data$log2FoldChange, label = goi),
        colour = "#2c3e50", fontface = "bold", size = 4, min.segment.length = 0
      )
  }
  g
}


#' Genes co-expressed with the gene of interest
#'
#' Pearson correlation of every gene against the GOI across all samples, on the
#' variance-stabilised scale. Useful as a cheap, annotation-free way to find
#' candidate partners of the GOI in this particular dataset.
#'
#' Note that with few samples and two clearly separated groups the top hits are
#' largely genes that follow the same group split, so treat this as a
#' hypothesis generator rather than evidence of co-regulation.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param n Integer. Number of genes to show.
#' @param goi Character. Gene to correlate against; defaults to the pipeline's GOI.
#' @param min_var Numeric. Genes with variance at or below this value are dropped.
#' @return A ggplot object, or `NULL` when the GOI is absent from the matrix.
#'   The underlying table is attached as the `"data"` attribute.
#' @examples
#' \dontrun{
#' p <- plot_top_correlated(res, n = 25)
#' head(attr(p, "data"))
#' }
#' @export
plot_top_correlated <- function(res_output, n = 20, goi = res_output$goi, min_var = 0) {
  mat <- .vsd_matrix(res_output)
  if (base::is.null(mat) || base::is.null(goi) || !goi %in% base::rownames(mat)) {
    base::warning("plot_top_correlated(): GOI not present in the expression matrix.")
    return(NULL)
  }
  if (base::ncol(mat) < 3) {
    base::warning("plot_top_correlated(): at least three samples are needed.")
    return(NULL)
  }

  keep <- base::is.finite(.row_vars(mat)) & .row_vars(mat) > min_var
  mat <- mat[keep, , drop = FALSE]
  if (!goi %in% base::rownames(mat)) {
    base::warning("plot_top_correlated(): GOI has no variance across samples.")
    return(NULL)
  }

  r <- base::as.vector(stats::cor(base::t(mat), mat[goi, ]))
  base::names(r) <- base::rownames(mat)
  r <- r[base::names(r) != goi]
  r <- r[base::is.finite(r)]
  if (base::length(r) == 0) {
    base::warning("plot_top_correlated(): no finite correlations.")
    return(NULL)
  }

  top <- r[base::order(base::abs(r), decreasing = TRUE)][base::seq_len(base::min(n, base::length(r)))]
  df <- base::data.frame(
    gene = base::names(top),
    r = base::as.numeric(top),
    direction = base::ifelse(top >= 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )
  df$gene <- stats::reorder(base::factor(df$gene), df$r)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$r, y = .data$gene, fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("positive" = "#b2182b", "negative" = "#2166ac"), name = NULL
    ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::labs(
      title = base::paste0("Genes co-expressed with ", goi),
      subtitle = "Pearson correlation on variance-stabilised counts",
      x = "correlation with the gene of interest", y = NULL
    ) +
    theme_goi() +
    ggplot2::theme(legend.position = "bottom")

  base::attr(g, "data") <- df[base::order(-base::abs(df$r)), c("gene", "r", "direction")]
  g
}


#' Library sizes and size factors
#'
#' Total raw counts per sample with the DESeq2 size factor printed on each bar.
#' Wildly uneven library sizes, or size factors far from one, are worth knowing
#' about before interpreting a single gene.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param intgroup Character. Column of `colData` used for colouring.
#' @param palette Character vector of fill colours, one per group.
#' @return A ggplot object, or `NULL` when the count matrix is missing.
#' @examples
#' \dontrun{
#' plot_library_sizes(res)
#' }
#' @export
plot_library_sizes <- function(res_output,
                               intgroup = "condition",
                               palette = c("#4f8832", "#f79c18")) {
  counts <- res_output$counts
  if (base::is.null(counts) || base::ncol(counts) == 0) {
    base::warning("plot_library_sizes(): no count matrix available.")
    return(NULL)
  }

  df <- base::data.frame(
    sample = base::colnames(counts),
    libsize = base::colSums(counts) / 1e6,
    stringsAsFactors = FALSE
  )

  sf <- base::tryCatch(DESeq2::sizeFactors(res_output$dds), error = function(e) NULL)
  df$sizeFactor <- if (!base::is.null(sf)) {
    base::round(base::as.numeric(sf[base::match(df$sample, base::names(sf))]), 2)
  } else {
    NA_real_
  }

  cd <- .coldata(res_output)
  if (!base::is.null(cd) && intgroup %in% base::names(cd)) {
    df$group <- base::as.factor(cd[[intgroup]][base::match(df$sample, base::rownames(cd))])
  } else {
    df$group <- base::factor("all")
  }
  df$sample <- base::factor(df$sample, levels = df$sample)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$sample, y = .data$libsize, fill = .data$group)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::labs(
      title = "Library sizes",
      subtitle = "labels show the DESeq2 size factor",
      x = NULL, y = "million assigned reads", fill = intgroup
    ) +
    theme_goi() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (base::any(base::is.finite(df$sizeFactor))) {
    g <- g + ggplot2::geom_text(ggplot2::aes(label = .data$sizeFactor),
      vjust = -0.4, size = 3, colour = "grey25", na.rm = TRUE
    )
  }
  if (base::length(base::levels(df$group)) <= base::length(palette)) {
    g <- g + ggplot2::scale_fill_manual(values = palette)
  }
  g
}


#' Diseases associated with the gene of interest
#'
#' Lollipop chart of the `DEGsToDiseases` scores already computed by the
#' pipeline. The pipeline writes this table to disk; this turns it into the
#' figure the README promises.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param n Integer. Maximum number of diseases to show.
#' @return A ggplot object, or `NULL` when no disease table is available.
#' @examples
#' \dontrun{
#' plot_disease_associations(res)
#' }
#' @export
plot_disease_associations <- function(res_output, n = 10) {
  x <- res_output$DEGsToDiseases
  if (base::is.null(x) || base::nrow(x) == 0 || base::ncol(x) < 2) {
    base::warning("plot_disease_associations(): no disease table available.")
    return(NULL)
  }

  df <- base::data.frame(
    disease = base::as.character(x[[1]]),
    score = base::suppressWarnings(base::as.numeric(x[[2]])),
    stringsAsFactors = FALSE
  )
  df <- df[base::is.finite(df$score) & base::nzchar(df$disease), , drop = FALSE]
  if (base::nrow(df) == 0) {
    base::warning("plot_disease_associations(): no usable rows.")
    return(NULL)
  }
  df <- df[base::order(df$score, decreasing = TRUE), , drop = FALSE]
  df <- df[base::seq_len(base::min(n, base::nrow(df))), , drop = FALSE]
  df$disease <- stats::reorder(base::factor(df$disease), df$score)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$score, y = .data$disease)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$score, yend = .data$disease),
      colour = "grey70", linewidth = 0.6
    ) +
    ggplot2::geom_point(size = 3.5, colour = "#f79c18") +
    ggplot2::labs(
      title = base::paste0("Diseases associated with ", res_output$goi),
      subtitle = "Open Targets association score (via KnowSeq)",
      x = "association score", y = NULL
    ) +
    theme_goi()
}


#' Volcano plot labelling the strongest genes
#'
#' Companion to the volcano plot the pipeline already produces: instead of
#' labelling only the GOI, this one also names the most significant genes, so
#' the plot works as a standalone summary figure.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param n_labels Integer. Number of top genes (by adjusted p-value) to label.
#' @param goi Character. Gene to always label; defaults to the pipeline's GOI.
#' @param lfcCutoff Numeric. Fold-change cutoff.
#' @param pCutoff Numeric. Adjusted p-value cutoff.
#' @return A ggplot object, or `NULL` when results are missing.
#' @examples
#' \dontrun{
#' plot_volcano_top(res, n_labels = 20)
#' }
#' @export
plot_volcano_top <- function(res_output,
                             n_labels = 15,
                             goi = res_output$goi,
                             lfcCutoff = 1,
                             pCutoff = 0.05) {
  res.df <- res_output$res.df
  if (base::is.null(res.df) || base::nrow(res.df) == 0) {
    base::warning("plot_volcano_top(): `res.df` not available.")
    return(NULL)
  }

  df <- base::data.frame(
    gene = base::rownames(res.df),
    log2FoldChange = res.df$log2FoldChange,
    padj = res.df$padj,
    stringsAsFactors = FALSE
  )
  df <- df[base::is.finite(df$log2FoldChange) & base::is.finite(df$padj), , drop = FALSE]
  df$padj[df$padj <= 0] <- base::min(df$padj[df$padj > 0], na.rm = TRUE)
  df$status <- base::ifelse(
    df$padj >= pCutoff | base::abs(df$log2FoldChange) <= lfcCutoff, "not significant",
    base::ifelse(df$log2FoldChange > 0, "up", "down")
  )

  ranked <- df[base::order(df$padj), , drop = FALSE]
  ranked <- ranked[ranked$status != "not significant", , drop = FALSE]
  lab <- ranked[base::seq_len(base::min(n_labels, base::nrow(ranked))), , drop = FALSE]
  if (!base::is.null(goi) && goi %in% df$gene) {
    lab <- base::unique(base::rbind(lab, df[df$gene == goi, , drop = FALSE]))
  }

  g <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$log2FoldChange,
    y = -base::log10(.data$padj)
  )) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$status), size = 0.6, alpha = 0.6) +
    ggplot2::scale_colour_manual(
      values = c("up" = "#b2182b", "down" = "#2166ac", "not significant" = "grey70"),
      name = NULL
    ) +
    ggplot2::geom_vline(
      xintercept = c(-lfcCutoff, lfcCutoff),
      linetype = "dashed", colour = "grey40", linewidth = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = -base::log10(pCutoff),
      linetype = "dashed", colour = "grey40", linewidth = 0.3
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(
      title = "Volcano plot",
      subtitle = base::paste0(
        "|log2FC| > ", lfcCutoff, " and padj < ", pCutoff,
        "  |  ", base::sum(df$status == "up"), " up, ",
        base::sum(df$status == "down"), " down"
      ),
      x = "log2 fold change", y = "-log10 adjusted p-value"
    ) +
    theme_goi() +
    ggplot2::theme(legend.position = "bottom")

  if (base::nrow(lab) > 0) {
    g <- g + ggrepel::geom_text_repel(
      data = lab,
      ggplot2::aes(label = .data$gene),
      size = 3, max.overlaps = Inf, min.segment.length = 0,
      fontface = base::ifelse(lab$gene == goi, "bold", "plain"),
      show.legend = FALSE
    )
  }
  g
}


#' Build every exploratory plot for a finished run
#'
#' Runs the QC and exploration plots on an existing `res_output`, adds them to
#' the list, and optionally writes them to disk. Each plot is attempted
#' independently: one failure never aborts the rest.
#'
#' This is called for you at the end of [pipeline()], but it is also the cheap
#' way to redraw the figures (with different cutoffs, colours or gene counts)
#' from a run you already have in memory.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param outdir Character or `NULL`. Directory to write PNGs into. Defaults to
#'   `res_output$directory`; pass `NULL` to skip writing files.
#' @param lfcCutoff,pCutoff Numeric cutoffs passed through to the plots that use them.
#' @param palette Character vector of group colours.
#' @param n_heatmap Integer. Number of genes in the DEG heatmap.
#' @param n_correlated Integer. Number of genes in the co-expression plot.
#' @param verbose Logical. Report progress with `message()`.
#' @return `res_output` with the extra plots (and their `*Path` entries) added.
#' @examples
#' \dontrun{
#' res <- run_pipeline(input = "counts_dir", dataType = "counts", goi = "TP53")
#' res <- build_extra_plots(res, n_heatmap = 60)
#' res$PCAplot
#' }
#' @export
build_extra_plots <- function(res_output,
                              outdir = res_output$directory,
                              lfcCutoff = 1,
                              pCutoff = 0.05,
                              palette = c("#4f8832", "#f79c18"),
                              n_heatmap = 40,
                              n_correlated = 20,
                              verbose = TRUE) {
  stopifnot(base::is.list(res_output))
  goi <- res_output$goi
  suffix <- if (base::is.null(goi)) "" else base::paste0("_", goi)

  specs <- base::list(
    list(
      slot = "PCAplot", file = "pca.png", w = 7, h = 5.5,
      fun = function() plot_pca(res_output, palette = palette)
    ),
    list(
      slot = "SampleDistancePlot", file = "sampleDistances.png", w = 6.5, h = 5.5,
      fun = function() plot_sample_distances(res_output)
    ),
    list(
      slot = "DEGHeatmap", file = base::paste0("degHeatmap", suffix, ".png"), w = 7, h = 8,
      fun = function() plot_deg_heatmap(res_output, n = n_heatmap, goi = goi)
    ),
    list(
      slot = "PvalueHistogram", file = "pvalueHistogram.png", w = 6.5, h = 4.5,
      fun = function() plot_pvalue_histogram(res_output, fill = palette[1])
    ),
    list(
      slot = "DispersionPlot", file = "dispersion.png", w = 6.5, h = 5,
      fun = function() plot_dispersion(res_output)
    ),
    list(
      slot = "RankPlot", file = base::paste0("rankPlot", suffix, ".png"), w = 7, h = 5,
      fun = function() {
        plot_goi_rank(res_output, goi = goi, lfcCutoff = lfcCutoff, pCutoff = pCutoff)
      }
    ),
    list(
      slot = "CorrelatedGenesPlot", file = base::paste0("coexpressed", suffix, ".png"), w = 6.5, h = 6,
      fun = function() plot_top_correlated(res_output, n = n_correlated, goi = goi)
    ),
    list(
      slot = "LibrarySizePlot", file = "librarySizes.png", w = 7, h = 4.5,
      fun = function() plot_library_sizes(res_output, palette = palette)
    ),
    list(
      slot = "DiseasePlot", file = base::paste0("diseases", suffix, ".png"), w = 7, h = 5,
      fun = function() plot_disease_associations(res_output)
    ),
    list(
      slot = "VolcanoLabelled", file = base::paste0("volcanoLabelled", suffix, ".png"), w = 7.5, h = 6,
      fun = function() {
        plot_volcano_top(res_output, goi = goi, lfcCutoff = lfcCutoff, pCutoff = pCutoff)
      }
    )
  )

  for (s in specs) {
    g <- base::tryCatch(
      base::suppressWarnings(s$fun()),
      error = function(e) {
        base::warning("Could not build ", s$slot, ": ", e$message)
        NULL
      }
    )
    if (base::is.null(g)) {
      if (verbose) base::message("  - skipped ", s$slot)
      next
    }
    res_output[[s$slot]] <- g
    p <- .save_plot(g, outdir, s$file, width = s$w, height = s$h)
    if (!base::is.null(p)) res_output[[base::paste0(s$slot, "Path")]] <- p
    if (verbose) base::message("  - ", s$slot)
  }

  res_output
}
