utils::globalVariables(c(
  "padj", "log2FoldChange", "ScientificName", "condition",
  "count", "threshold", "genelabels", "baseMean",
  "abrScientificName", "excludedSample", "txi", "rel"
))


#' Launch the GOI Explorer Shiny app
#' @export
run_app <- function() {
  app_dir <- base::system.file("shiny/goiExplorer_app", package = "goiExplorer")
  if (app_dir == "") {
    stop("Could not find Shiny app. Try re-installing the package.")
  }
  shiny::runApp(app_dir)
}


#' Export a data.frame to TXT and XLSX
#'
#' Writes a tab-delimited TXT and an Excel file side-by-side.
#'
#' @param data A data.frame or tibble to save.
#' @param outdir Path to an existing directory.
#' @param prefix File name prefix (defaults to the name of `data`).
#' @return Invisible `NULL`. Files are written to `outdir`.
#' @details
#' This will call `utils::write.table()` and `writexl::write_xlsx()`.
#' If `prefix` does not exist it will be created based on the data name.
#' @examples
#' \dontrun{
#' df <- matrix(20, 10, 5)
#' write_txt_xlsx(df, tempdir())
#'}
#' @importFrom writexl write_xlsx
#' @export
write_txt_xlsx <- function(data, outdir, prefix = base::deparse(base::substitute(data)) ) {
  stopifnot(is.data.frame(data))
  utils::write.table(as.data.frame(data), 
                     file = base::paste0(outdir, prefix, ".txt"), 
                     quote = F, sep = "\t", row.names = FALSE)
  writexl::write_xlsx(base::data.frame(data), base::paste0(outdir, prefix,".xlsx"))
}










#' Get annotated DESeq2 results
#'
#' Fetches `results()` from a `DESeqDataSet`, adds a “significant” flag,
#' and joins in external gene annotations via biomaRt.
#'
#' @param dds A `DESeqDataSet` object after running `DESeq()`  
#' @param rowNamesOfCounts Character. The column in your annotation dataset (e.g. `"external_gene_name"`)  
#'   used to match against the DESeq results row names.  
#' @param mart A `biomaRt` Mart object (created with `useMart()`) for annotation lookup.  
#' @param pCutoff Numeric scalar. The `alpha` threshold for adjusted p-values (default `0.05`).  
#' @param pAdjustMethod Character. Method for p-value adjustment (default `"BH"`).  
#' @return A `data.frame` with one row per gene, including all `results()` columns  
#'   plus annotation fields (`external_gene_name`, `entrezgene_id`, etc.), sorted by `padj`.  
#' @examples
#' \dontrun{
#' library(DESeq2); library(biomaRt)
#' mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#' dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
#' dds <- DESeq(dds)
#' df   <- get_res_df(dds, row_key = "ensembl_gene_id", mart = mart)
#' head(df)
#' }
#' @importFrom DESeq2 results
#' @importFrom biomaRt getBM
#' @export
get_res_df <- function(dds, 
                       rowNamesOfCounts, 
                       mart, 
                       pCutoff = 0.05, 
                       pAdjustMethod = "BH") {
  df <- DESeq2::results(dds, alpha = pCutoff, pAdjustMethod = pAdjustMethod, tidy = TRUE)
  
  # Set a boolean column for significance
  df$significant <- ifelse(df$padj < pCutoff, "Significant", NA)
  
  df2 <- biomaRt::getBM(attributes=c(rowNamesOfCounts, "external_gene_name", 
                                     'ensembl_gene_id', 
                            "entrezgene_id", 'gene_biotype', "description", 
                            'chromosome_name', 
                            'start_position', 'end_position', 'strand'), 
               filters = rowNamesOfCounts, 
               values = df$row, 
               mart = mart)
  
  colnames(df)[1] <- rowNamesOfCounts
  df <- base::merge(df, df2, all.x=TRUE)
  df <- df[!duplicated(df[,1]), ]
  res.df <- df[order(df$padj), ]
  rownames(res.df) <- res.df[,1]
  return(res.df)
}
