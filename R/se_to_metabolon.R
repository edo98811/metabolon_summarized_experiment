#' Convert SummarizedExperiment to Metabolon Format
#'
#' This function converts a `SummarizedExperiment` object into a format compatible
#' with Metabolon data analysis.
#'
#' @param se A `SummarizedExperiment` object containing the assay data.
#' @param cdt A character string specifying the path to the client data table path.
#' @param output_file A character string specifying the path to save the output file.
#'   If `NULL` and `save_file` is `TRUE`, a default filename will be generated and it will be saved in the working directory.
#' @param save_file A logical value indicating whether to save the output to a file. Needs to be csv!
#'   Default is `TRUE`.
#' @param sample_id_column A character string specifying the column name in the colData
#'   that contains the sample IDs matching the PARENT_SAMPLE_NAME names. Default is `NULL`.
#' @param metabolon_sample_column A character string specifying the column name in the client data table that needs
#' to be matched with the one in the SummarizedExperiment. The default is "CLIENT_SAMPLE_ID".
#'
#' @examples
#' \dontrun{
#' se <- SummarizedExperiment::SummarizedExperiment()
#' result <- se_to_metabolon(se, cdt = "cdt.xlsx", input_features = "gene_symbol", save_file = TRUE)
#' }
#'
#' @importFrom readxl read_excel
#' @importFrom utils write.table
#' @importFrom SummarizedExperiment colData assay
#' @export
se_to_metabolon <- function(se,
                            cdt,
                            output_file = NULL,
                            save_file = TRUE,
                            se_sample_column = NULL, 
                            metabolon_sample_column = "CLIENT_SAMPLE_ID") {

  # Check correctness of input
  if (!inherits(se, "SummarizedExperiment")) stop("The input object is not a SummarizedExperiment.")
  if (is.null(output_file) && save_file) output_file <- paste0("se_to_metabolon_", Sys.Date(), ".csv")

  # Read metadata
  metadata_metabolon <- make_metadata(readxl::read_excel(cdt, sheet = 3))
  metadata_se <- as.data.frame(colData(se))


  ids_to_match <- if (!is.null(se_sample_column)) {
    if (!(se_sample_column %in% colnames(metadata_se))) {
      stop("The specified se_sample_column does not exist in the colData of the SummarizedExperiment.")
    }
    metadata_se[[se_sample_column]]
  } else {
    colnames(se)
  }

  # Read assay
  assay <- assay(se)

  # Add PARENT_SAMPLE_NAME column (from documentation)
  assay_dataframe <- data.frame(t(assay))
  assay_dataframe$PARENT_SAMPLE_NAME <- metadata_metabolon$PARENT_SAMPLE_NAME[match(ids_to_match, metadata_metabolon[[metabolon_sample_column]])]
  assay_dataframe <- assay_dataframe[!is.na(assay_dataframe$PARENT_SAMPLE_NAME),, drop = FALSE] # if no matching then dont do anything
  if (nrow(assay_dataframe) == 0) stop("Sample ids not matching with sample ids from metabolon")
  assay_dataframe <- assay_dataframe[c("PARENT_SAMPLE_NAME", setdiff(names(assay_dataframe), "PARENT_SAMPLE_NAME"))] # Put the col PARENT_SAMPLE_NAME as first
  
  
  # Write the output to a file
  if (save_file) {
    message("Saving results to: ", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE)
    write.table(assay_dataframe, output_file, sep = ",", row.names = FALSE)
  }
  return(assay_dataframe)
}
