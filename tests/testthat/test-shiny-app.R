app_server_fun <- function() {
  app_dir <- system.file("shiny/goiExplorer_app", package = "goiExplorer")
  skip_if(app_dir == "", "shiny app not installed")
  source(file.path(app_dir, "server.R"), local = new.env())$value
}

test_that("the shipped app files load and expose ui and server objects", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyFiles")
  skip_if_not_installed("DT")

  app_dir <- system.file("shiny/goiExplorer_app", package = "goiExplorer")
  skip_if(app_dir == "", "shiny app not installed")

  ui <- source(file.path(app_dir, "ui.R"), local = new.env())$value
  srv <- source(file.path(app_dir, "server.R"), local = new.env())$value

  expect_s3_class(ui, "shiny.tag.list")
  expect_true(is.function(srv))
  expect_named(formals(srv), c("input", "output", "session"))
})

test_that("the UI declares every control the server listens to", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyFiles")
  skip_if_not_installed("DT")

  app_dir <- system.file("shiny/goiExplorer_app", package = "goiExplorer")
  skip_if(app_dir == "", "shiny app not installed")
  html <- as.character(source(file.path(app_dir, "ui.R"), local = new.env())$value)

  # the AI panel used to be referenced by the server but missing from the UI
  for (id in c(
    "ai_query", "ai_ask", "ai_clear", "ai_provider", "ai_key", "ai_model",
    "ai_chat", "summary_cards", "download_ui", "run_summary"
  )) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
})

test_that("the chat records the exchange and clears on request", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyFiles")
  skip_if_not_installed("fs")
  withr::local_envvar(c(ANTHROPIC_API_KEY = "", OPENAI_API_KEY = ""))

  shiny::testServer(app_server_fun(), {
    session$setInputs(
      ai_provider = "anthropic", ai_key = "", ai_model = "",
      ai_query = "Is my gene a strong hit?"
    )
    session$setInputs(ai_ask = 1)

    msgs <- chat()
    expect_length(msgs, 2)
    expect_equal(msgs[[1]]$role, "user")
    expect_equal(msgs[[1]]$content, "Is my gene a strong hit?")
    expect_equal(msgs[[2]]$role, "assistant")
    # no key is configured, so the failure must come back as a message
    expect_match(msgs[[2]]$content, "AI request failed")

    # a second question keeps the history and stays alternating
    session$setInputs(ai_query = "And what about the pathways?")
    session$setInputs(ai_ask = 2)
    expect_length(chat(), 4)
    expect_equal(vapply(chat(), function(m) m$role, character(1)),
      c("user", "assistant", "user", "assistant"))

    session$setInputs(ai_clear = 1)
    expect_length(chat(), 0)
  })
})

test_that("an empty question is ignored", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyFiles")
  skip_if_not_installed("fs")

  shiny::testServer(app_server_fun(), {
    session$setInputs(ai_provider = "anthropic", ai_key = "", ai_model = "", ai_query = "   ")
    expect_error(session$setInputs(ai_ask = 1), NA)
    expect_length(chat(), 0)
  })
})

test_that("status text reacts before any run", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinyFiles")
  skip_if_not_installed("fs")

  shiny::testServer(app_server_fun(), {
    expect_match(output$status, "Waiting for you")
    expect_equal(output$output_path, "")
  })
})
