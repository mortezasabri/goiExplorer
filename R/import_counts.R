#' Import the counts-based RNA-Seq dataset
#'
#' Reads in a plain counts table plus metadata (or both via a `list`), create dds object,
#'
#' @param input Character. Path to a folder containing  
#'   - `counts.txt`: gene × sample matrix  
#'   - `gr.txt`: one group/condition per line  
#' @param abr_healthy abbreviation for healthy group; for example H for healthy
#' @param abr_case abbreviation for case group; for example D for disease
#' @return A `dds`: a `DESeqDataSet` object  
#' @details
#' Internally this function:  
#' 1. Reads the two input files via `read.table()`: counts.txt and gr.txt or list of these two 
#' 2. Constructs a `DESeqDataSet` and return `dds` 
#' @examples
#' \dontrun{
#' pipeline_counts(
#'   input = "inst/extdata",
#'   outdir   = "counts_output"
#' )
#' }
#' @importFrom utils read.table
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom stats relevel
#' @export
import_counts <- function(input,
                          abr_healthy,
                          abr_case) {
  
  if (is.character(input)) {
    gr <- read.table(file.path(input, "gr.txt"), header = F)[,1]
    counts <- read.table(file.path(input, "counts.txt"), 
                         header = T, 
                         row.names = 1, 
                         check.names = F)
  } else if (is.list(input)) { counts <- input$counts; gr <- input$gr }
  ## ----coldata object
  healthy <- grep(abr_healthy, colnames(counts))
  case <- grep(abr_case, colnames(counts))
  
  gr <- as.factor(gr)
  gr <- relevel(gr, abr_healthy)
  coldata <- data.frame(
    name = colnames(counts),
    condition = gr,
    row.names = colnames(counts)
  )

  
  
  
  ## ----DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ condition)
  stopifnot(inherits(dds, "DESeqDataSet"))
  return(dds)
}
