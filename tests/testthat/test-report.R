test_that("write_report produces a self-contained HTML file", {
  res <- fake_res_output()
  td <- withr::local_tempdir()
  f <- file.path(td, "report.html")

  expect_message(out <- write_report(res, file = f), "Report written")
  expect_equal(out, f)
  expect_true(file.exists(f))

  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "^<!DOCTYPE html>")
  expect_match(html, "</html>")
  expect_match(html, res$goi, fixed = TRUE)
  expect_match(html, "Differentially expressed genes")
  expect_match(html, "Associated diseases")
  expect_match(html, "Run summary")
  # the stylesheet must be inlined, not linked
  expect_match(html, "<style>")
  expect_no_match(html, "<link[^>]*stylesheet")
})

test_that("write_report embeds the figures it is given", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("knitr")

  res <- fake_res_output()
  td <- withr::local_tempdir()
  res <- build_extra_plots(res, outdir = td, verbose = FALSE)

  f <- write_report(res, file = file.path(td, "report.html"))
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")

  expect_match(html, "data:image/png;base64,")
  expect_match(html, "Principal component analysis")
})

test_that("write_report survives a run with almost nothing in it", {
  td <- withr::local_tempdir()
  f <- write_report(list(goi = "TP53"), file = file.path(td, "bare.html"))

  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "TP53", fixed = TRUE)
  expect_match(html, "was not found in the results table")
})

test_that("HTML special characters in the data are escaped", {
  res <- fake_res_output()
  res$res.df$description[1] <- "<script>alert(1)</script> & \"quoted\""
  res$DEGsToDiseases$disease[1] <- "<b>bold</b>"

  td <- withr::local_tempdir()
  f <- write_report(res, file = file.path(td, "escaped.html"))
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")

  expect_no_match(html, "<script>alert", fixed = TRUE)
  expect_match(html, "&lt;script&gt;alert", fixed = TRUE)
  expect_no_match(html, "<b>bold</b>", fixed = TRUE)
})
