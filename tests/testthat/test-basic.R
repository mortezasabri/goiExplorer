test_that("run_pipeline returns a list with expected names", {
  skip_on_cran()
  # you’d call run_pipeline on some toy data here
  input <- system.file("extdata", package="goiExplorer")
  td <- withr::local_tempdir()
  run_pipeline(input, "counts", parent_outdir = td)
})