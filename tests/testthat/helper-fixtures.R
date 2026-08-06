# A synthetic `res_output` that looks enough like a real run for the plotting,
# summary and report code, without needing DESeq2, biomaRt or the network.
fake_res_output <- function(n_genes = 200, n_samples = 8, goi = "GOI1") {
  set.seed(42)
  half <- n_samples %/% 2
  samples <- c(paste0("s", seq_len(half), "_H"), paste0("s", seq_len(half), "_D"))
  # MIRROR tracks the GOI exactly, ANTI is its exact opposite: the two anchors
  # the co-expression test relies on
  genes <- c(goi, "MIRROR", "ANTI", paste0("G", seq_len(n_genes - 3)))

  # variance-stabilised-looking matrix; the first 20 genes track the group
  mat <- matrix(stats::rnorm(n_genes * n_samples, mean = 8, sd = 1.5),
    nrow = n_genes, dimnames = list(genes, samples)
  )
  grp <- rep(c(0, 1), each = half)
  mat[1, ] <- mat[1, ] + 3 * grp
  mat["MIRROR", ] <- mat[goi, ]
  mat["ANTI", ] <- 20 - mat[goi, ]
  for (i in 4:20) mat[i, ] <- mat[i, ] + ifelse(i %% 2 == 0, 2.5, -2.5) * grp

  counts <- round(2^mat)
  storage.mode(counts) <- "integer"

  res.df <- data.frame(
    baseMean = rowMeans(counts),
    log2FoldChange = stats::rnorm(n_genes, sd = 0.3),
    lfcSE = stats::runif(n_genes, 0.1, 0.5),
    stat = stats::rnorm(n_genes),
    pvalue = stats::runif(n_genes, 0.05, 1),
    padj = stats::runif(n_genes, 0.2, 1),
    gene_biotype = "protein_coding",
    entrezgene_id = seq_len(n_genes),
    chromosome_name = "1",
    description = "a synthetic gene used in the tests",
    row.names = genes,
    stringsAsFactors = FALSE
  )
  res.df$log2FoldChange[1:20] <- c(3, rep(c(2.4, -2.4), length.out = 19))
  res.df$pvalue[1:20] <- 1e-8
  res.df$padj[1:20] <- 1e-6
  res.df <- res.df[order(res.df$padj), ]

  degs <- subset(res.df, padj < 0.05 & abs(log2FoldChange) > 1)

  geneCounts <- data.frame(
    count = counts[goi, ],
    condition = factor(rep(c("H", "D"), each = half), levels = c("H", "D")),
    row.names = samples
  )

  list(
    goi = goi,
    ensemblSpecies = "hsapiens_gene_ensembl",
    directory = NULL,
    counts = counts,
    vsd = mat,
    res.df = res.df,
    degs = degs,
    geneCounts = geneCounts,
    pvalue_raw = res.df$pvalue,
    DEGsToDiseases = data.frame(
      disease = c("disease A", "disease B", "disease C"),
      score = c(0.9, 0.6, 0.2),
      stringsAsFactors = FALSE
    ),
    KEGGpaths = c("/tmp/paths/hsa04064.NF-kappa B signaling pathway.png")
  )
}

# Every exported plotting function that takes only `res_output`.
plot_fun_names <- function() {
  c(
    "plot_pca", "plot_sample_distances", "plot_deg_heatmap",
    "plot_pvalue_histogram", "plot_dispersion", "plot_goi_rank",
    "plot_top_correlated", "plot_library_sizes",
    "plot_disease_associations", "plot_volcano_top"
  )
}
