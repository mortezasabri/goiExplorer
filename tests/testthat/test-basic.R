test_that("run_pipeline returns a list with expected names", {
  skip_on_cran()
  # the full pipeline queries Ensembl, KEGG and Open Targets
  skip_if_offline()
  for (p in c("DESeq2", "biomaRt", "KnowSeq", "limma", "KEGGREST", "pathview")) {
    skip_if_not_installed(p)
  }

  input <- system.file("extdata", package = "goiExplorer")
  td <- withr::local_tempdir()
  res <- run_pipeline(input, "counts", parent_outdir = td, report = FALSE)

  expect_type(res, "list")
  for (slot in c("goi", "dds", "counts", "res", "vsd", "res.df", "degs", "geneCounts")) {
    expect_true(!is.null(res[[slot]]), info = slot)
  }
  expect_s3_class(res$Volcanoplot, "ggplot")
  expect_s3_class(res$PCAplot, "ggplot")
})

test_that("write_txt_xlsx writes both files", {
  skip_if_not_installed("writexl")
  td <- paste0(withr::local_tempdir(), .Platform$file.sep)
  df <- data.frame(gene = c("A", "B"), value = c(1.5, 2.5))

  write_txt_xlsx(df, td, prefix = "demo")

  expect_true(file.exists(file.path(td, "demo.txt")))
  expect_true(file.exists(file.path(td, "demo.xlsx")))
  expect_equal(read.table(file.path(td, "demo.txt"), header = TRUE, sep = "\t"), df)
})
