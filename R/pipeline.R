#' Run the counts-based RNA-Seq pipeline for the gene-of-interest
#'
#' Reads in a plain counts table plus metadata, performs DESeq2 analysis,
#' enrichment, and writes results to disk.
#'
#' @param dds a `dds` object from DESeq2
#' @param goi a gene of interest
#' @param parent_outdir Character. Directory where all outputs (tables, plots) will be written.
#' If it does not exist, it will be created recursively.
#' @param abr_healthy abbreviation for healthy group; for example H for healthy
#' @param abr_case abbreviation for case group; for example D for disease
#' @param ensemblSpecies e.g. hsapiens_gene_ensembl. You can get that from biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))$dataset
#' @param rowNamesOfCounts the type of the rownames in count table, e.g. external_gene_name
#' @param interestedTerms give priority to interested terms (e.g. cardio or inflammation) in otargen::colocalisationsForGene
#' @param intgroup the column name of the group, e.g. `condition`
#' @param lfcCutoff cutoff on logFC
#' @param pCutoff cutoff on adj pvalue
#' @param pAdjustMethod method to adjusting the pvalues. see ?p.adjust
#' @param palette palette to fill the plots
#' @return Files saved in parent_outdir:
#'   - `res.df`: a `DESeqDataSet` object
#'   - `degs`: a `data.frame` of annotated DE results
#'   - plots: boxplot, barplot, Count Plot, Volcano, MA plot, Top10 associated diseases, Pathways, colocalisationsForGene result for the GOI
#' @details
#' Internally this function:
#' 1. Reads the dds
#' 2. runs `DESeq()`
#' 3. Annotates and shrinks logfold changes with `get_res_df()`
#' @examples
#' \dontrun{
#' pipeline_counts(
#'   dds = "inst/extdata",
#'   outdir = "counts_output"
#' )
#' }
#' @importFrom DESeq2 resultsNames DESeq results counts plotMA
#' @importFrom biomaRt getBM
#' @importFrom ggplot2 ggplot ggsave aes labs scale_y_continuous rel geom_point
#' @importFrom magrittr %>%
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom ggrepel geom_text_repel
#' @importFrom rstatix add_xy_position t_test
#' @importFrom KnowSeq DEGsToDiseases
#' @importFrom ggpubr ggbarplot ggboxplot stat_pvalue_manual
#' @importFrom pathview pathview
#' @importFrom utils data
#' @export
pipeline <- function(dds,
                     goi,
                     parent_outdir,
                     abr_healthy,
                     abr_case,
                     ensemblSpecies,
                     rowNamesOfCounts,
                     interestedTerms,
                     intgroup,
                     lfcCutoff,
                     pCutoff,
                     pAdjustMethod,
                     palette) {
  options(width = 100, timeout = 600)
  stopifnot(inherits(dds, "DESeqDataSet"))
  stopifnot(inherits(goi, "character"))
  stopifnot(inherits(parent_outdir, "character"))
  stopifnot(inherits(abr_healthy, "character"))
  stopifnot(inherits(abr_case, "character"))
  stopifnot(inherits(ensemblSpecies, "character"))
  stopifnot(inherits(rowNamesOfCounts, "character"))
  stopifnot(inherits(interestedTerms, "character"))
  stopifnot(inherits(intgroup, "character"))
  stopifnot(inherits(lfcCutoff, "numeric"))
  stopifnot(inherits(pCutoff, "numeric"))
  stopifnot(inherits(pAdjustMethod, "character"))
  stopifnot(inherits(palette, "character"))

  if (ensemblSpecies == "hsapiens_gene_ensembl") {
    organism <- "human"
    ScientificName <- "Homo sapiens"
    abrScientificName <- "hsa"
  } else {
    stop('ensemblSpecies should be one of: \n
         biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))$dataset \n
         and for now the only valid value is "hsapiens_gene_ensembl"')
  }

  res_output <- list() # saving the directories and results in this as an output of run_pipeline()
  res_output$goi <- goi
  res_output$ensemblSpecies <- ensemblSpecies



  # ensure we have a real output directory
  if (missing(parent_outdir) || is.null(parent_outdir)) {
    parent_outdir <- file.path(tempdir(), paste0("goiExplorer_", Sys.time()))
  }
  if (!dir.exists(parent_outdir)) dir.create(parent_outdir, recursive = TRUE)

  # Ensure trailing slash (use .Platform$file.sep for cross-platform)
  if (substr(parent_outdir, nchar(parent_outdir), nchar(parent_outdir)) != .Platform$file.sep) {
    parent_outdir <- paste0(parent_outdir, .Platform$file.sep)
  }
  res_output$directory <- outdir <- parent_outdir

  # Try the *default* Ensembl site first, then fall back to a cached file
  mart_path <- system.file("extdata", "hsa_mart.rds", package = "goiExplorer")
  mart <- tryCatch(
    biomaRt::useEnsembl("ensembl", dataset = ensemblSpecies),
    error = function(e_live) {
      message("Default Ensembl failed, trying US mirror")
      warning(
        "Live Ensembl query failed: ", e_live$message,
        "\nFalling back to local cache in extdata."
      )
      # Fall back to your shipped RDS
      base::readRDS(mart_path)
    }
  )

  ## ----DESeq2
  # The condition of interest should go at the end
  message("Number of genes before omiting low expressed genes: ")
  nrow(dds)
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 1
  dds <- dds[keep, ]
  message("Number of genes after omiting low expressed genes: ")
  nrow(dds)

  dds <- suppressMessages(DESeq2::DESeq(dds))
  dds <- DESeq2::estimateSizeFactors(dds)
  counts <- DESeq2::counts(dds)

  res_output$dds <- dds
  res_output$counts <- counts

  # ----------- DE
  DESeq2::resultsNames(dds)
  message("Summary of differentially expressed genes: ")
  summary(res <- DESeq2::results(dds, alpha = pCutoff, pAdjustMethod = pAdjustMethod))
  res_output$res <- res


  # ----------- Normalization
  vsd <- DESeq2::vst(dds, blind = FALSE)
  res_output$vsd <- vsd


  ## ----DESeq2 out
  res.df <- get_res_df(dds, rowNamesOfCounts, mart, pCutoff, pAdjustMethod)
  res$pvalue <- ifelse(is.na(res$pvalue), 1, res$pvalue)
  res$padj <- ifelse(is.na(res$padj), 1, res$padj)
  res.df$pvalue <- ifelse(is.na(res.df$pvalue), 1, res.df$pvalue)
  res.df$padj <- ifelse(is.na(res.df$padj), 1, res.df$padj)

  # degs: only those passed the cutoffs
  degs <- subset(res.df, padj < pCutoff & abs(log2FoldChange) > lfcCutoff)
  write_txt_xlsx(degs, outdir,
    prefix = paste0("degs_", abr_case, "_Vs_", abr_healthy)
  )
  write_txt_xlsx(res.df, outdir,
    prefix = paste0("de_", abr_case, "_Vs_", abr_healthy)
  )

  res_output$degs <- degs
  res_output$res.df <- res.df

  ## ----abstract
  attributes <- c(
    rowNamesOfCounts, "external_gene_name",
    "ensembl_gene_id",
    "gene_biotype", "description"
  )
  des <- tryCatch(
    biomaRt::getBM(
      attributes = attributes,
      filters    = rowNamesOfCounts,
      values     = goi,
      mart       = mart
    ),
    error = function(e) {
      warning(
        "Could not fetch annotation from Ensembl: ", e$message,
        "  skipping annotation."
      )
      # Return an empty data.frame with the right columns
      base::data.frame(
        base::matrix(
          ncol = base::length(attributes), nrow = 0,
          dimnames = base::list(NULL, attributes)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  if (rowNamesOfCounts == "external_gene_name") {
    des <- des[, -2]
  }
  res_output$des <- des

  entrezgene_id <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
    keys = goi,
    keytype = "SYMBOL", column = "ENTREZID"
  )
  res_output$entrezgene_id <- entrezgene_id

  link <- unique(paste0(
    "https://useast.ensembl.org/",
    sub(" ", "_", ScientificName),
    "/Gene/Summary?g=", des$ensembl_gene_id
  ))


  df <- geneCounts <-
    DESeq2::plotCounts(
      dds,
      gene = goi,
      intgroup = intgroup,
      returnData = TRUE
    )
  res_output$geneCounts <- geneCounts


  ## ----Barplot
  stat.test <- df %>% rstatix::t_test(count ~ condition)
  stat.test$p <- base::round(res.df[goi, ]$padj, 3)
  stat.test$statistic <- base::round(res.df[goi, ]$stat, 3)
  stat.test <- rstatix::add_significance(stat.test)
  stat.test.mean <- stat.test %>% add_xy_position(fun = "mean_sd", x = intgroup)

  g <- ggbarplot(df,
    x = intgroup, y = "count", add = "mean_sd",
    fill = intgroup,
    palette = palette
  ) +
    stat_pvalue_manual(stat.test.mean,
      label = "padj: {p} {p.signif}",
      tip.length = 0.01
    ) +
    labs(title = paste0("Barplot for ", goi), y = "Normalized counts")
  p <- paste0(outdir, "barplot_", goi, ".png")
  ggsave(p,
    g,
    dpi = 300
  )
  res_output$Barplot <- g
  res_output$BarplotPath <- p


  ## ----Boxplot
  stat.test.max <- stat.test %>% add_xy_position(x = intgroup)

  rm(g)
  g <- ggpubr::ggboxplot(df,
    x = intgroup, y = "count",
    fill = intgroup,
    palette = palette
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.75),
      ggplot2::aes(group = condition),
      color = "steelblue4"
    ) +
    ggpubr::stat_pvalue_manual(stat.test.max,
      label = "padj: {p} {p.signif}",
      vjust = -1, bracket.nudge.y = 1
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::labs(title = paste0("Boxplot for ", goi), y = "Normalized counts")
  p <- paste0(outdir, "boxplot_", goi, ".png")
  ggplot2::ggsave(p,
    g,
    dpi = 300
  )
  res_output$Boxplot <- g
  res_output$BoxplotPath <- p


  ## ----Count Plot
  rm(g)
  g <- ggplot2::ggplot(
    geneCounts,
    ggplot2::aes(x = condition, y = count, color = condition)
  ) +
    ggplot2::scale_y_log10() +
    ggbeeswarm::geom_beeswarm(cex = 3) +
    ggplot2::ggtitle(goi) +
    ggplot2::theme(legend.position = "none")
  p <- paste0(outdir, "plotCounts_", goi, ".png")
  ggplot2::ggsave(p,
    g,
    dpi = 300
  )
  res_output$Countplot <- g
  res_output$CountplotPath <- p

  ## ----Volcano
  res.df$threshold <- res.df$padj < pCutoff & abs(res.df$log2FoldChange) > lfcCutoff
  res.df$genelabels <- base::rownames(res.df) %in% goi # name of genes want to display
  base::options(ggrepel.max.overlaps = Inf) # in case of ggrepel: 1 unlabeled data points (too many overlaps). Consider increasing max.overlaps


  rm(g)
  g <- ggplot2::ggplot(res.df) +
    ggplot2::geom_point(ggplot2::aes(x = log2FoldChange, y = -log10(padj), colour = threshold),
      size = 0.4
    ) +
    ggrepel::geom_text_repel(aes(
      x = log2FoldChange, y = -log10(padj),
      label = ifelse(genelabels == T,
        base::rownames(res.df), ""
      )
    )) +
    ggplot2::ggtitle(paste0("Volcano Plot (", goi, ")")) +
    ggplot2::xlab("log2 fold change") +
    ggplot2::ylab("-log10 adjusted p-value") +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = ggplot2::rel(1.5), hjust = 0.5),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.25))
    )
  p <- paste0(outdir, "volcanoPlot_", goi, ".png")
  ggplot2::ggsave(p,
    g,
    dpi = 300
  )
  res_output$Volcanoplot <- g
  res_output$VolcanoplotPath <- p

  ## ----MA plot
  # df is your results table, with at least:
  #   baseMean, log2FoldChange, padj, genelabels
  df <- as.data.frame(res)
  # 1. Define significance
  sig_cutoff <- df$padj < pCutoff & abs(df$log2FoldChange) > lfcCutoff
  df$significant <- ifelse(sig_cutoff, "yes", "no")

  # 2. Build the MA‐plot
  g <- ggplot2::ggplot(df, ggplot2::aes(x = baseMean, y = log2FoldChange)) +
    # log‐scale the x axis
    ggplot2::scale_x_log10() +

    # all points
    geom_point(ggplot2::aes(color = df$significant), alpha = 0.5, size = 0.5) +

    # highlight significant genes in a distinct color
    ggplot2::scale_color_manual(
      values = c(
        no = "gray60",
        yes = "red" # or palette[["sig"]]
      ),
      guide = ggplot2::guide_legend(title = "Significant")
    ) +

    # optional horizontal line at LFC = 0
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +

    # annotate only the top N genes (or all significant ones)
    ggrepel::geom_text_repel(
      data = base::subset(df, base::rownames(df) %in% goi),
      ggplot2::aes(label = goi),
      size = 6,
      max.overlaps = 10,
      box.padding = 0.3,
      segment.size = 0.2
    ) +
    ggplot2::labs(
      x = "Mean of normalized counts (log10)",
      y = "log2 Fold Change",
      title = base::paste("MA plot for", goi)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title      = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
  # Save the plot
  p <- base::paste0(outdir, "MAplot_", goi, ".png")
  ggplot2::ggsave(p,
    plot = g,
    dpi = 300,
    width = 6,
    height = 5
  )
  res_output$plotMA <- g
  res_output$test <- base::subset(df[, c("log2FoldChange", "padj", "significant")], base::rownames(df) %in% goi)
  ## ----DEGsToDiseases
  x <- base::data.frame(KnowSeq::DEGsToDiseases(goi, size = 10, getEvidences = TRUE))[, 1:2]
  x[, 2] <- base::as.numeric(x[, 2])
  x[, 2] <- base::round(x[, 2], 3)
  write_txt_xlsx(x, outdir,
    prefix = paste0("DEGsToDiseases_", goi)
  )
  res_output$DEGsToDiseases <- x



  ## ----Pathways
  outKEGG <- paste0(outdir, "KEGG/")
  if (!dir.exists(outKEGG)) dir.create(outKEGG, recursive = TRUE)
  tmp <- getwd()
  setwd(outKEGG)

  geneList <- degs$log2FoldChange
  names(geneList) <- degs$entrezgene_id
  geneList <- geneList[stats::complete.cases(names(geneList))]
  geneList <- geneList[order(geneList, decreasing = T)]

  z <- limma::getGeneKEGGLinks(abrScientificName)
  path.ids <- z[z$GeneID == entrezgene_id, ]$PathwayID

  for (i in 1:length(path.ids)) {
    names(path.ids)[i] <- KEGGREST::keggGet(path.ids[i])[[1]]$PATHWAY
  }

  for (i in 1:length(path.ids)) {
    tryCatch(
      {
        var2 <- path.ids[i]
        nm <- paste0(var2, ": ", names(var2))
        message(paste0("\n *** Initiating loop ", i, "/", length(path.ids), ": ", nm))

        geneOnThePathID <- KEGGREST::keggGet(path.ids[i])[[1]]$GENE
        geneOnThePathID <- geneOnThePathID[seq(1, length(geneOnThePathID), by = 2)]
        geneOnThePathID <- geneList[which(names(geneList) %in% geneOnThePathID)]
        ceil <- base::round(max(abs(geneOnThePathID)), 1)
        if (requireNamespace("pathview", quietly = TRUE)) {
          utils::data("bods", package = "pathview", envir = environment())
          tryCatch(
            {
              pathview(
                gene.data = geneList,
                pathway.id = path.ids[i],
                species = abrScientificName,
                out.suffix = names(var2),
                kegg.native = T,
                limit = list(gene = ceil, cpd = 1)
              )
            },
            error = function(e) {
              warning("Could not generate KEGG plot for ", names(var2), ": ", e$message)
            }
          )
        } else {
          warning("pathview package not installed; skipping KEGG plots")
        }
      },
      error = function(e) {
      }
    )
  }

  dir.create("paths")
  filesstrings::file.move(
    list.files(pattern = paste0(
      abrScientificName,
      "[0-9]*\\..*\\.png"
    )),
    "paths"
  ) # hsa for human
  setwd(tmp)
  l <- list.files(file.path(outKEGG, "paths"), full.names = T)
  res_output$KEGGpaths <- l



  x <- otargen::colocalisationsForGene(goi)
  x <- as.data.frame(x)
  if (nrow(x) > 0) {
    if (length(x) > 10) {
      ind <- grep(paste(interestedTerms, collapse = "|"),
        x$Trait_reported,
        ignore.case = TRUE
      )
      message("Here are the most interested ones")
      write_txt_xlsx(x[ind, ], outdir,
        prefix = paste0("colocalisationsForGene_", goi)
      )
    } else {
      message("Here are the all traits associated with this gene")
      write_txt_xlsx(x, outdir,
        prefix = paste0("colocalisationsForGene_", goi)
      )
    }
    res_output$colocalisationsForGene <- x
  }




  if (length(x$Study) == 1) {
    grDevices::png(paste0(outdir, "manhattan1_", goi, ".png"),
      units = "in", height = 10, width = 14, res = 200
    )
    otargen::plot_manhattan(otargen::manhattan(x$Study[1]))
    grDevices::dev.off()
  } else if (length(x$Study) == 2) {
    grDevices::png(paste0(outdir, "manhattan1_", goi, ".png"),
      units = "in", height = 10, width = 14, res = 200
    )
    otargen::plot_manhattan(otargen::manhattan(x$Study[1]))
    grDevices::dev.off()

    grDevices::png(paste0(outdir, "manhattan2_", goi, ".png"),
      units = "in", height = 10, width = 14, res = 200
    )
    otargen::plot_manhattan(otargen::manhattan(x$Study[2]))
    grDevices::dev.off()
  }


  return(res_output)
}
