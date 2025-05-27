if (!requireNamespace("openai", quietly = TRUE)) {
  stop("Please install the openai package to use ai_agent().")
}


#' Chat with your GOI Explorer results
#'
#' @param prompt A natural language question about your last pipeline run.
#' @param model OpenAI model to use.
#' @return A character response.
#' @export
ai_agent <- function(prompt, model = "gpt-4") {
  stopifnot(nzchar(Sys.getenv("OPENAI_API_KEY")))
  tryCatch({
    openai::create_chat_completion(
      model = model,
      messages = list(
        list(role = "system",
             content = "You are a bioinformatics assistant helping interpret info regarding all sort of analysis regarding gene-of-interest (GOI) in a dataset."),
        list(role = "user", content = prompt)
      )
    )$choices[[1]]$message$content
  }, error = function(e) paste("AI request failed:", e$message))
}
