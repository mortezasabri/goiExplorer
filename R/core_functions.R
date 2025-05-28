#' Execute the full RNA-Seq analysis pipeline
#'
#' Dispatches either the “counts” or “salmon” workflow based on `dataType`,
#' then runs differential expression, annotation, KEGG, and plots.
#'
#' @param input Character. Path (or list of two) to the directory containing input files:
#'   - for `"counts"`: `counts.txt`, `gr.txt` or a list of two slots: `counts` (table) and `gr` (vector)
#'   - for `"salmon"`: `ID.txt`, `ID.txt`, `gr.txt`, plus any Salmon `quant.sf` dirs
#' @param dataType Character, one of `"counts"` or `"salmon"` indicating
#'   which underlying pipeline to use.
#' @param goi a gene of interest
#' @param parent_outdir    Character. Directory where all outputs (tables, plots) will be written.
#'   If it doesn’t exist, it will be created recursively.
#' @param abr_healthy abbreviation for healthy group; for example H for healthy
#' @param abr_case abbreviation for case group; for example D for disease
#' @param ensemblSpecies e.g. hsapiens_gene_ensembl. You can get that from biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))$dataset
#' @param rowNamesOfCounts the type of the rownames in count table, e.g. external_gene_name
#' @param interestedTerms give priority to interested terms (e.g. cardio or inflammation) in otargen::colocalisationsForGene
#' @param intgroup the column name of the group, e.g. `condition`
#' @param excludedSample optional: the samples you want to exclude. should be the names of the corresponded column
#' @param lfcCutoff cutoff on logFC
#' @param pCutoff cutoff on adj pvalue
#' @param pAdjustMethod method to adjusting the pvalues. see ?p.adjust
#' @param palette palette to fill the plots
#' @param tx2gene   Path to an `.rds` file or data.frame mapping transcripts -> genes.
#'   a directory also will be created, icluding all outputs.
#' @return A `list` with elements:
#'   - `dds`: the `DESeqDataSet` object  
#'   - `results`: a `data.frame` of annotated DE results  
#'   - `enrichment`: the `enrichResult` object from ORA  
#' @examples
#' \dontrun{
#' # counts-based run
#' meta <- data.frame(
#'   sample = readLines(system.file("extdata","ID.txt",package="goiExplorer")),
#'   group  = readLines(system.file("extdata","gr.txt",package="goiExplorer"))
#' )
#' run_pipeline(
#'   input     = system.file("extdata", package = "goiExplorer"),
#'   dataType  = "counts",
#'   parent_outdir    = "counts_output/"
#' )
#' }
#' @export
run_pipeline <- function(input, 
                         dataType, 
                         goi = "CYLD",
                         parent_outdir,
                         abr_healthy = "H",
                         abr_case = "D",
                         ensemblSpecies = "hsapiens_gene_ensembl",
                         rowNamesOfCounts = "external_gene_name",
                         interestedTerms = c("cardi", "inflammat"),
                         intgroup = c("condition"),
                         excludedSample = NA,
                         lfcCutoff = 1,
                         pCutoff = 0.05,
                         pAdjustMethod = "fdr", 
                         palette = c("#4f8832", "#f79c18"),
                         tx2gene = NA) {
  lfcCutoff <- as.numeric(lfcCutoff)
  pCutoff <- as.numeric(pCutoff)
  if (!exists("parent_outdir")) {
    parent_outdir <- file.path(tempdir(check = T), "goiExplorer_output")
  }
  if(!dir.exists(parent_outdir)) dir.create(parent_outdir, recursive = TRUE)
  
  if (!exists("input")) {
    stop("Please define the input directory")
  }
  
  if (!dataType %in% c("counts", "salmon")) {
    stop("dataType must be either 'counts' or 'salmon'")
  }
  
  
  
  
  if (dataType == "salmon") {
    message("Reading gr.txt and ID.txt inside the quant (input) directory")
    
    dds <- import_salmon(input,
                         abr_healthy,
                         abr_case,
                         tx2gene)

  } else if (dataType == "counts") {
    
    dds <- import_counts(input,
                  abr_healthy,
                  abr_case)
  }
  
  res_output <- pipeline(dds, 
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
                         palette)
  return(res_output)
  message("Outputs are being saved in: ", parent_outdir)
}

