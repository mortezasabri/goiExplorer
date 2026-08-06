test_that("summarise_run reports the numbers that matter", {
  res <- fake_res_output()
  txt <- summarise_run(res)

  expect_type(txt, "character")
  expect_length(txt, 1)

  expect_match(txt, "goiExplorer run summary")
  expect_match(txt, res$goi, fixed = TRUE)
  expect_match(txt, "H \\(n=4\\) vs D \\(n=4\\)")
  expect_match(txt, "Significant genes: 20 \\(11 up, 9 down\\)")
  expect_match(txt, "higher in the case group")
  expect_match(txt, "disease A")
  expect_match(txt, "NF-kappa B signaling pathway", fixed = TRUE)
})

test_that("summarise_run copes with an empty or partial run", {
  expect_match(summarise_run(NULL), "No goiExplorer results")
  expect_match(summarise_run(list()), "run summary")
  expect_match(summarise_run(list(goi = "TP53")), "TP53", fixed = TRUE)
})

test_that("summarise_run never leaks the count matrix", {
  res <- fake_res_output()
  txt <- summarise_run(res)

  # a handful of raw counts should not appear anywhere in the digest
  vals <- as.character(res$counts[5:10, 1])
  expect_false(any(vapply(vals, function(v) grepl(v, txt, fixed = TRUE), logical(1))))
})

test_that("ai_agent reports a missing key instead of stopping", {
  withr::local_envvar(c(ANTHROPIC_API_KEY = "", OPENAI_API_KEY = ""))

  out <- ai_agent("what happened?")
  expect_type(out, "character")
  expect_match(out, "^AI request failed:")
  expect_match(out, "ANTHROPIC_API_KEY")
})

test_that("ai_agent validates its provider", {
  out <- ai_agent("hello", provider = "not-a-provider", api_key = "x")
  expect_match(out, "^AI request failed:")
})

test_that("ai_agent rejects an empty prompt", {
  out <- ai_agent("", api_key = "x")
  expect_match(out, "^AI request failed:")
})

test_that("the default provider follows whichever key is set", {
  withr::with_envvar(
    c(ANTHROPIC_API_KEY = "a", OPENAI_API_KEY = "b"),
    expect_equal(goiExplorer:::.ai_default_provider(), "anthropic")
  )
  withr::with_envvar(
    c(ANTHROPIC_API_KEY = "", OPENAI_API_KEY = "b"),
    expect_equal(goiExplorer:::.ai_default_provider(), "openai")
  )
})

test_that("each provider has a default model", {
  expect_type(goiExplorer:::.ai_default_model("anthropic"), "character")
  expect_type(goiExplorer:::.ai_default_model("openai"), "character")
  expect_error(goiExplorer:::.ai_default_model("nope"))
})
