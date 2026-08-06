## ---------------------------------------------------------------------------
## A small LLM client so you can ask questions about a finished run in plain
## English. Supports Anthropic (Claude) and OpenAI; the API key never leaves
## your machine except in the request to the provider you chose.
##
## Nothing here is loaded at package build time and no provider package is
## required to install goiExplorer.
## ---------------------------------------------------------------------------

.ai_default_model <- function(provider) {
  switch(provider,
    anthropic = "claude-sonnet-5",
    openai    = "gpt-4o-mini",
    base::stop("Unknown provider: ", provider)
  )
}

.ai_default_provider <- function() {
  if (base::nzchar(base::Sys.getenv("ANTHROPIC_API_KEY"))) {
    return("anthropic")
  }
  if (base::nzchar(base::Sys.getenv("OPENAI_API_KEY"))) {
    return("openai")
  }
  "anthropic"
}

.ai_api_key <- function(provider, api_key = NULL) {
  if (!base::is.null(api_key) && base::nzchar(api_key)) {
    return(api_key)
  }
  env <- if (provider == "anthropic") "ANTHROPIC_API_KEY" else "OPENAI_API_KEY"
  key <- base::Sys.getenv(env)
  if (!base::nzchar(key)) {
    base::stop(
      "No API key found. Set ", env, " (e.g. Sys.setenv(", env,
      " = \"...\")) or pass it as `api_key`."
    )
  }
  key
}

.ai_system_prompt <- function() {
  base::paste(
    "You are a bioinformatics analyst helping a researcher interpret an RNA-Seq",
    "differential expression run centred on one gene of interest (GOI).",
    "You are given a summary of the actual results; base every quantitative",
    "claim on it and say plainly when the summary does not contain the answer.",
    "Do not invent gene names, p-values or fold changes.",
    "Be concise and concrete, use the numbers you were given, and flag caveats",
    "(small sample size, multiple testing, correlation vs causation) where they",
    "genuinely matter."
  )
}

# Format one section of the results digest; returns character(0) when empty.
.ai_section <- function(title, body) {
  body <- body[base::nzchar(body)]
  if (base::length(body) == 0) {
    return(base::character(0))
  }
  c(base::paste0("## ", title), body, "")
}

#' Summarise a pipeline run as plain text
#'
#' Turns the `res_output` list into a compact, human-readable digest: the GOI's
#' statistics, how many genes moved, the strongest genes in each direction, and
#' the disease / pathway annotations. This is what [ai_agent()] sends to the
#' model as context, and it is a useful thing to read (or paste into a lab
#' notebook) on its own.
#'
#' Only the summary is sent to the LLM provider - never the count matrix.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param max_genes Integer. How many genes to list per direction.
#' @return A single character string.
#' @examples
#' \dontrun{
#' res <- run_pipeline(input = "counts_dir", dataType = "counts", goi = "TP53")
#' cat(summarise_run(res))
#' }
#' @export
summarise_run <- function(res_output, max_genes = 15) {
  if (base::is.null(res_output) || !base::is.list(res_output)) {
    return("No goiExplorer results are available yet.")
  }
  goi <- res_output$goi
  out <- c("# goiExplorer run summary", "")

  ## --- design -------------------------------------------------------------
  design <- base::character(0)
  if (!base::is.null(goi)) design <- c(design, base::paste0("Gene of interest: ", goi))
  if (!base::is.null(res_output$ensemblSpecies)) {
    design <- c(design, base::paste0("Species dataset: ", res_output$ensemblSpecies))
  }
  gc <- res_output$geneCounts
  if (!base::is.null(gc) && "condition" %in% base::names(gc)) {
    tab <- base::table(gc$condition)
    design <- c(design, base::paste0(
      "Groups: ",
      base::paste0(base::names(tab), " (n=", base::as.integer(tab), ")", collapse = " vs ")
    ))
  }
  if (!base::is.null(res_output$counts)) {
    design <- c(design, base::paste0(
      "Genes retained after filtering: ", base::nrow(res_output$counts)
    ))
  }
  out <- c(out, .ai_section("Design", design))

  ## --- the gene of interest ----------------------------------------------
  goi_lines <- base::character(0)
  res.df <- res_output$res.df
  if (!base::is.null(res.df) && !base::is.null(goi) && goi %in% base::rownames(res.df)) {
    row <- res.df[goi, , drop = FALSE]
    lfc <- row$log2FoldChange[1]
    goi_lines <- c(
      goi_lines,
      base::paste0("log2 fold change: ", base::round(lfc, 3)),
      base::paste0("adjusted p-value: ", base::signif(row$padj[1], 3)),
      base::paste0("raw p-value: ", base::signif(row$pvalue[1], 3)),
      base::paste0("base mean: ", base::round(row$baseMean[1], 1)),
      base::paste0(
        "direction: ",
        if (base::is.na(lfc)) "unknown" else if (lfc > 0) "higher in the case group" else "higher in the healthy group"
      )
    )
    for (col in c("gene_biotype", "entrezgene_id", "chromosome_name", "description")) {
      if (col %in% base::names(row) && !base::is.na(row[[col]][1])) {
        goi_lines <- c(goi_lines, base::paste0(col, ": ", row[[col]][1]))
      }
    }
  }
  if (!base::is.null(gc) && "condition" %in% base::names(gc) && "count" %in% base::names(gc)) {
    means <- base::tapply(gc$count, gc$condition, base::mean)
    goi_lines <- c(goi_lines, base::paste0(
      "mean normalised counts: ",
      base::paste0(base::names(means), " = ", base::round(base::as.numeric(means), 1), collapse = ", ")
    ))
  }
  out <- c(out, .ai_section(base::paste0("Gene of interest (", goi, ")"), goi_lines))

  ## --- differential expression -------------------------------------------
  de_lines <- base::character(0)
  degs <- res_output$degs
  if (!base::is.null(degs) && base::nrow(degs) > 0) {
    up <- degs[!base::is.na(degs$log2FoldChange) & degs$log2FoldChange > 0, , drop = FALSE]
    dn <- degs[!base::is.na(degs$log2FoldChange) & degs$log2FoldChange < 0, , drop = FALSE]
    de_lines <- c(de_lines, base::paste0(
      "Significant genes: ", base::nrow(degs),
      " (", base::nrow(up), " up, ", base::nrow(dn), " down)"
    ))
    fmt <- function(d, n) {
      if (base::nrow(d) == 0) {
        return("")
      }
      d <- d[base::order(d$padj), , drop = FALSE]
      d <- d[base::seq_len(base::min(n, base::nrow(d))), , drop = FALSE]
      base::paste(
        base::paste0(
          base::rownames(d), " (log2FC ", base::round(d$log2FoldChange, 2),
          ", padj ", base::signif(d$padj, 2), ")"
        ),
        collapse = "; "
      )
    }
    de_lines <- c(
      de_lines,
      base::paste0("Top up-regulated: ", fmt(up, max_genes)),
      base::paste0("Top down-regulated: ", fmt(dn, max_genes))
    )
  } else if (!base::is.null(res.df)) {
    de_lines <- c(de_lines, "No genes passed the significance cutoffs.")
  }
  out <- c(out, .ai_section("Differential expression", de_lines))

  ## --- co-expression ------------------------------------------------------
  if (!base::is.null(res_output$CorrelatedGenesPlot)) {
    cor_df <- base::attr(res_output$CorrelatedGenesPlot, "data")
    if (!base::is.null(cor_df) && base::nrow(cor_df) > 0) {
      cor_df <- cor_df[base::seq_len(base::min(max_genes, base::nrow(cor_df))), , drop = FALSE]
      out <- c(out, .ai_section(
        "Genes most correlated with the GOI",
        base::paste(
          base::paste0(base::as.character(cor_df$gene), " (r=", base::round(cor_df$r, 2), ")"),
          collapse = "; "
        )
      ))
    }
  }

  ## --- diseases -----------------------------------------------------------
  dis <- res_output$DEGsToDiseases
  if (!base::is.null(dis) && base::nrow(dis) > 0 && base::ncol(dis) >= 2) {
    out <- c(out, .ai_section(
      "Associated diseases (KnowSeq / Open Targets)",
      base::paste0(
        base::seq_len(base::nrow(dis)), ". ", base::as.character(dis[[1]]),
        " (score ", dis[[2]], ")"
      )
    ))
  }

  ## --- pathways -----------------------------------------------------------
  paths <- res_output$KEGGpaths
  if (!base::is.null(paths) && base::length(paths) > 0) {
    nm <- base::basename(base::unlist(paths))
    nm <- base::sub("\\.png$", "", nm)
    nm <- base::sub("^([a-z]{3}[0-9]+)\\.", "\\1: ", nm)
    out <- c(out, .ai_section(
      "KEGG pathways containing the GOI",
      base::unique(nm)
    ))
  }

  ## --- GWAS ---------------------------------------------------------------
  coloc <- res_output$colocalisationsForGene
  if (!base::is.null(coloc) && base::nrow(coloc) > 0 && "Trait_reported" %in% base::names(coloc)) {
    traits <- base::unique(base::as.character(coloc$Trait_reported))
    traits <- traits[base::seq_len(base::min(20, base::length(traits)))]
    out <- c(out, .ai_section("GWAS traits colocalising with the GOI", traits))
  }

  base::paste(out, collapse = "\n")
}


# --- provider back-ends -----------------------------------------------------

.ai_require_http <- function() {
  missing <- base::character(0)
  for (p in c("httr", "jsonlite")) {
    if (!base::requireNamespace(p, quietly = TRUE)) missing <- c(missing, p)
  }
  if (base::length(missing) > 0) {
    base::stop(
      "ai_agent() needs the ", base::paste(missing, collapse = " and "),
      " package(s): install.packages(c(", base::paste0("\"", missing, "\"", collapse = ", "), "))"
    )
  }
}

.ai_call_anthropic <- function(messages, system, model, key, max_tokens, temperature) {
  body <- base::list(
    model = model,
    max_tokens = max_tokens,
    temperature = temperature,
    system = system,
    messages = messages
  )
  resp <- httr::POST(
    url = "https://api.anthropic.com/v1/messages",
    httr::add_headers(
      "x-api-key" = key,
      "anthropic-version" = "2023-06-01",
      "content-type" = "application/json"
    ),
    body = jsonlite::toJSON(body, auto_unbox = TRUE, null = "null"),
    encode = "raw",
    httr::timeout(180)
  )
  parsed <- jsonlite::fromJSON(
    httr::content(resp, as = "text", encoding = "UTF-8"),
    simplifyVector = FALSE
  )
  if (httr::status_code(resp) >= 300) {
    base::stop(
      "Anthropic API error (", httr::status_code(resp), "): ",
      if (!base::is.null(parsed$error$message)) parsed$error$message else "unknown error"
    )
  }
  blocks <- parsed$content
  txt <- base::vapply(
    blocks,
    function(b) if (base::identical(b$type, "text")) b$text else "",
    base::character(1)
  )
  base::paste(txt[base::nzchar(txt)], collapse = "\n")
}

.ai_call_openai <- function(messages, system, model, key, max_tokens, temperature) {
  msgs <- base::c(base::list(base::list(role = "system", content = system)), messages)
  body <- base::list(
    model = model,
    max_tokens = max_tokens,
    temperature = temperature,
    messages = msgs
  )
  resp <- httr::POST(
    url = "https://api.openai.com/v1/chat/completions",
    httr::add_headers(
      "Authorization" = base::paste("Bearer", key),
      "content-type" = "application/json"
    ),
    body = jsonlite::toJSON(body, auto_unbox = TRUE, null = "null"),
    encode = "raw",
    httr::timeout(180)
  )
  parsed <- jsonlite::fromJSON(
    httr::content(resp, as = "text", encoding = "UTF-8"),
    simplifyVector = FALSE
  )
  if (httr::status_code(resp) >= 300) {
    base::stop(
      "OpenAI API error (", httr::status_code(resp), "): ",
      if (!base::is.null(parsed$error$message)) parsed$error$message else "unknown error"
    )
  }
  parsed$choices[[1]]$message$content
}


#' Ask a question about your GOI Explorer results
#'
#' Sends your question, plus a text summary of the run (see [summarise_run()]),
#' to a large language model and returns the answer. The model can therefore
#' talk about *your* numbers instead of guessing: "is my gene one of the
#' strongest changes?", "what do the top down-regulated genes have in common?",
#' "how should I phrase this for a figure legend?".
#'
#' Two providers are supported. The default is Anthropic (Claude) when
#' `ANTHROPIC_API_KEY` is set, otherwise OpenAI when `OPENAI_API_KEY` is set.
#'
#' Only the summary text is transmitted - the count matrix and the raw results
#' table never leave your machine. Treat the answer as a reading aid, not as a
#' result: verify anything you plan to publish against the tables the pipeline
#' wrote to disk.
#'
#' @param prompt Character. Your question.
#' @param res_output Optional. The list returned by [run_pipeline()]; when given,
#'   its summary is prepended as context.
#' @param history Optional list of previous turns, each a list with `role`
#'   (`"user"` or `"assistant"`) and `content`. Lets you hold a conversation.
#' @param provider `"anthropic"` or `"openai"`. Defaults to whichever key is set.
#' @param model Character. Model name; defaults to a sensible model per provider.
#' @param api_key Character. Overrides the environment variable.
#' @param max_tokens Integer. Maximum length of the reply.
#' @param temperature Numeric. Sampling temperature; 0 keeps answers reproducible.
#' @return A character string: the model's answer, or a message starting with
#'   `"AI request failed:"` when the call could not be made.
#' @examples
#' \dontrun{
#' Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")
#' res <- run_pipeline(input = "counts_dir", dataType = "counts", goi = "TP53")
#'
#' ai_agent("Is TP53 among the strongest changes here, or just significant?", res)
#'
#' # multi-turn
#' q1 <- "Summarise the result in three sentences."
#' a1 <- ai_agent(q1, res)
#' ai_agent("Now rewrite that for a general audience.",
#'   res,
#'   history = list(
#'     list(role = "user", content = q1),
#'     list(role = "assistant", content = a1)
#'   )
#' )
#' }
#' @export
ai_agent <- function(prompt,
                     res_output = NULL,
                     history = NULL,
                     provider = NULL,
                     model = NULL,
                     api_key = NULL,
                     max_tokens = 1024,
                     temperature = 0) {
  base::tryCatch(
    {
      stopifnot(base::is.character(prompt), base::length(prompt) == 1, base::nzchar(prompt))
      provider <- base::match.arg(
        if (base::is.null(provider)) .ai_default_provider() else provider,
        c("anthropic", "openai")
      )
      .ai_require_http()
      key <- .ai_api_key(provider, api_key)
      if (base::is.null(model)) model <- .ai_default_model(provider)

      user_msg <- prompt
      if (!base::is.null(res_output)) {
        user_msg <- base::paste0(
          "Here are the results of my run:\n\n",
          summarise_run(res_output),
          "\n\n---\n\nMy question: ", prompt
        )
      }

      messages <- base::c(
        if (base::is.null(history)) base::list() else history,
        base::list(base::list(role = "user", content = user_msg))
      )

      if (provider == "anthropic") {
        .ai_call_anthropic(
          messages, .ai_system_prompt(), model, key, max_tokens, temperature
        )
      } else {
        .ai_call_openai(
          messages, .ai_system_prompt(), model, key, max_tokens, temperature
        )
      }
    },
    error = function(e) base::paste("AI request failed:", base::conditionMessage(e))
  )
}


#' Ask the model to interpret a finished run
#'
#' Convenience wrapper around [ai_agent()] with a ready-made prompt: what the
#' result says about the gene of interest, what stands out in the rest of the
#' data, and what to check next.
#'
#' @param res_output The list returned by [run_pipeline()] / [pipeline()].
#' @param ... Passed to [ai_agent()] (`provider`, `model`, `api_key`, ...).
#' @return A character string.
#' @examples
#' \dontrun{
#' cat(ai_interpret(res))
#' }
#' @export
ai_interpret <- function(res_output, ...) {
  ai_agent(
    base::paste(
      "Interpret this run for me. Cover, in short paragraphs:",
      "(1) what happened to the gene of interest and how convincing it is;",
      "(2) what the overall differential expression looks like and whether the",
      "run seems technically sound;",
      "(3) the two or three most useful follow-up analyses or validations.",
      "Reference specific genes and numbers from the summary."
    ),
    res_output = res_output,
    ...
  )
}
