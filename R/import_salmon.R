#' Run the Salmon-based RNA-Seq pipeline for the gene-of-interest
#'
#' Uses Salmon quantifications plus a tx2gene mapping to build counts,
#'
#' @param input Character vector of paths to Salmon output directories.
#' @param abr_healthy abbreviation for healthy group; for example H for healthy
#' @param abr_case abbreviation for case group; for example D for disease
#' @param tx2gene   Path to an `.rds` file or data.frame mapping transcripts -> genes.
#' @return A `dds`: a `DESeqDataSet` object  
#' @details
#' Internally this function:  
#' 1. Imports counts via `tximport::tximport()` and the provided `tx2gene` mapping  
#' 2. Constructs a `DESeqDataSet` and return `dds` 
#' @examples
#' \dontrun{
#' pipeline_salmon(
#'   salmonDirs = list.dirs("inst/extdata/salmon", recursive = FALSE),
#'   tx2gene    = system.file("extdata/tx2gene.rds", package="goiExplorer"),
#'   outdir     = "salmon_output"
#' )
#' }
#' @importFrom tximport tximport
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom stats relevel
#' @export
import_salmon <- function(input,
                          abr_healthy,
                          abr_case,
                          tx2gene) {

  
  ## ----preparing for importing the data
  samples <- read.table(file.path(input, "ID.txt"), header = F)[,1]
  files <- file.path(input, samples, "quant.sf")
  
  sampIDs <- gsub( paste0(input, "/(.+)/quant.sf"), "\\1", files)
  gr <- read.table(file.path(input, "gr.txt"), header = F)[,1]
  
  names(files) <- paste0(sampIDs, "_", gr)
  stopifnot( all(file.exists(files)) )
  
  # reordering H to D
  healthy <- grep(abr_healthy, gr)
  case <- grep(abr_case, gr)
  files <- files[c(healthy, case)]; gr <- gr[c(healthy, case)]
  
  ## ----importing the data
  # txi object
  if (is.na("tx2gene")) {
    tx2gene <- readRDS(system.file("extdata", "tx2gene.rds", package = "goiExplorer"))
  }
  
  if (!is.na(excludedSample)) {
    ind <- match(excludedSample, names(files))
    gr <- gr[-ind]
    files <- files[-ind]
    txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene) 
    # use --gencode on salmon index
  } else {
    txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene) 
    # use --gencode on salmon index
  }
  ## ----coldata object
  healthy <- grep(abr_healthy, colnames(txi$counts))
  case <- grep(abr_case, colnames(txi$counts))
  
  gr <- as.factor(gr)
  gr <- relevel(gr, abr_healthy)
  coldata <- data.frame(
    name = colnames(txi$counts),
    condition = gr,
    row.names = colnames(txi$counts)
  )

  
  
  
  ## ----DESeq2
  dds <- DESeq2::DESeqDataSetFromTximport(txi = txi,
                                  colData = coldata,
                                  design = ~ condition)
  stopifnot(inherits(dds, "DESeqDataSet"))
  return(dds)
  
}
