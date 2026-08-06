test_that("every plot function returns a ggplot for a complete run", {
  skip_if_not_installed("ggplot2")
  res <- fake_res_output()

  # plot_dispersion is the one that needs a real DESeqDataSet
  for (fn in setdiff(plot_fun_names(), "plot_dispersion")) {
    g <- get(fn)(res)
    expect_s3_class(g, "ggplot")
  }
})

test_that("plot functions return NULL instead of erroring on missing pieces", {
  skip_if_not_installed("ggplot2")
  empty <- list(goi = "GOI1")

  for (fn in plot_fun_names()) {
    expect_warning(g <- get(fn)(empty))
    expect_null(g, info = fn)
  }
})

test_that("the PCA separates two clearly different groups", {
  skip_if_not_installed("ggplot2")
  res <- fake_res_output()
  g <- plot_pca(res)

  df <- ggplot2::ggplot_build(g)$plot$data
  expect_true(all(c("PC1", "PC2", "sample", "group") %in% names(df)))
  expect_equal(nrow(df), ncol(res$vsd))

  # the healthy and disease samples should sit on opposite sides of PC1
  h <- mean(df$PC1[grepl("_H$", df$sample)])
  d <- mean(df$PC1[grepl("_D$", df$sample)])
  expect_gt(abs(h - d), 0)
  expect_lt(h * d, 0)
})

test_that("the heatmap always includes the gene of interest", {
  skip_if_not_installed("ggplot2")
  res <- fake_res_output()
  # push the GOI far down the ranking so it is not in the top n by itself
  res$res.df[res$goi, "padj"] <- 0.99
  res$res.df <- res$res.df[order(res$res.df$padj), ]

  expect_silent(g <- plot_deg_heatmap(res, n = 5))
  expect_s3_class(g, "ggplot")
  expect_true(res$goi %in% levels(g$data$gene))
  expect_match(g$labels$subtitle, paste(res$goi, "outlined"), fixed = TRUE)
})

test_that("the rank plot reports where the gene of interest sits", {
  skip_if_not_installed("ggplot2")
  res <- fake_res_output()
  g <- plot_goi_rank(res)

  expect_match(g$labels$subtitle, res$goi, fixed = TRUE)
  expect_match(g$labels$subtitle, "percentile")
  # the GOI has the largest log2FC in the fixture, so it ranks last
  expect_match(g$labels$subtitle, paste("rank", nrow(g$data)))
})

test_that("co-expression ranks the genes that track the GOI", {
  skip_if_not_installed("ggplot2")
  res <- fake_res_output()
  g <- plot_top_correlated(res, n = 10)

  tbl <- attr(g, "data")
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), 10)
  expect_false(res$goi %in% as.character(tbl$gene))
  expect_true(all(abs(tbl$r) <= 1))
  # sorted by absolute correlation, strongest first
  expect_equal(order(abs(tbl$r), decreasing = TRUE), seq_len(nrow(tbl)))

  # the fixture's perfect copy and perfect opposite must come out on top,
  # with the right signs
  expect_true(all(c("MIRROR", "ANTI") %in% as.character(tbl$gene[1:2])))
  expect_equal(tbl$r[as.character(tbl$gene) == "MIRROR"], 1, tolerance = 1e-6)
  expect_equal(tbl$r[as.character(tbl$gene) == "ANTI"], -1, tolerance = 1e-6)
  expect_equal(tbl$direction[as.character(tbl$gene) == "MIRROR"], "positive")
  expect_equal(tbl$direction[as.character(tbl$gene) == "ANTI"], "negative")
})

test_that("the volcano labels the strongest genes and the GOI", {
  skip_if_not_installed("ggplot2")
  res <- fake_res_output()
  g <- plot_volcano_top(res, n_labels = 5)

  expect_s3_class(g, "ggplot")
  expect_match(g$labels$subtitle, "up")
  expect_match(g$labels$subtitle, "down")
})

test_that("build_extra_plots fills the result list and writes the PNGs", {
  skip_if_not_installed("ggplot2")
  res <- fake_res_output()
  td <- withr::local_tempdir()

  out <- build_extra_plots(res, outdir = td, verbose = FALSE)

  expect_true(is.list(out))
  # everything except the dispersion plot, which needs a DESeqDataSet
  for (slot in c(
    "PCAplot", "SampleDistancePlot", "DEGHeatmap", "PvalueHistogram",
    "RankPlot", "CorrelatedGenesPlot", "LibrarySizePlot", "DiseasePlot",
    "VolcanoLabelled"
  )) {
    expect_s3_class(out[[slot]], "ggplot")
    expect_true(file.exists(out[[paste0(slot, "Path")]]), info = slot)
  }
  expect_null(out$DispersionPlot)

  # the originals must survive
  expect_identical(out$res.df, res$res.df)
  expect_identical(out$goi, res$goi)
})

test_that("build_extra_plots skips writing when outdir is NULL", {
  skip_if_not_installed("ggplot2")
  out <- build_extra_plots(fake_res_output(), outdir = NULL, verbose = FALSE)

  expect_s3_class(out$PCAplot, "ggplot")
  expect_null(out$PCAplotPath)
})

test_that("theme_goi is a ggplot theme", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(theme_goi(), "theme")
})
