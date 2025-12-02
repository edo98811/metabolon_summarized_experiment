#' Convert Results to Metabolon Format.
#' 
#' This function converts a results table into a format compatible with Metabolon data analysis.
#' 
#' @param table A data frame containing the results to be converted.
#' @param format A character string specifying the format of the input results.
#'   Options are `"DESeqResults"` (default) or `"topTable"`.
#' @param omics A character string specifying the type of omics data.
#'   Options are `"Transcriptomics"` (default) or `"Proteomics"`.
#' @param custom_colnames An optional character vector specifying custom column names for the output.
#'   If `NULL`, default column names will be used.
#' @param output_file A character string specifying the path to save the output file.
#'   If `NULL` and `save_file` is `TRUE`, a default filename will be generated. Needs to be csv!
#'  @param save_file A logical value indicating whether to save the output to a file.
#'   Default is `TRUE`.
#' @return A data frame formatted according to the standard required by Metabolon.
#' 
#' @examples
#' \dontrun{
#' # Convert DESeqResults to Metabolon format
#' results_metabolon <- results_to_metabolon(
#'   table = deseq_results,
#'   format = "DESeqResults",
#'   omics = "Transcriptomics",
#'   input_features = "gene_symbol",
#'   save_file = TRUE
#' )
#' }
#' @importFrom utils write.table
#' @export
results_to_metabolon <- function(table,
                                format = "topTable",
                                omics = "Transcriptomics",
                                custom_colnames = NULL,
                                output_file = NULL, 
                                save_file = T) {

  # Make sure output_file is set if save_file is TRUE
  if (is.null(output_file) && save_file) output_file <- paste0("results_to_metabolon_", Sys.Date(), ".csv")

  # Create the results for metabolon
  results_for_metabolon <- make_results_for_metabolon(table, format, omics, custom_colnames)

  # Write the output to a file
  if (save_file) {
    message("Saving results to: ", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE)
    write.table(results_for_metabolon, output_file, sep = ",", row.names = FALSE)
  } 
  return(results_for_metabolon)
}