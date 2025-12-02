#' Convert SummarizedExperiment to Metabolon Format
#'
#' This function converts a `SummarizedExperiment` object into a format compatible
#' with Metabolon data analysis.
#'
#' @param se A `SummarizedExperiment` object containing the assay data.
#' @param cdt A character string specifying the path to the client data table path.
#' @param output_file A character string specifying the path to save the output file.
#'   If `NULL` and `save_file` is `TRUE`, a default filename will be generated and it will be saved in the working directory.
#' @param annotation_db An AnnotationDbi OrgDb object (e.g., org.Hs.eg.db::org.Hs.eg.db) for mapping input features to ENSEMBL IDs.
#' @param save_file A logical value indicating whether to save the output to a file. Needs to be csv!
#'   Default is `TRUE`.
#' @param sample_id_column A character string specifying the column name in the colData
#'   that contains the sample IDs matching the PARENT_SAMPLE_NAME names. Default is `NULL`.
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
                            sample_id_column = NULL) {

  # Check correctness of input
  if (!inherits(se, "SummarizedExperiment")) stop("The input object is not a SummarizedExperiment.")
  if (is.null(output_file) && save_file) output_file <- paste0("se_to_metabolon_", Sys.Date(), ".csv")

  # Read metadata
  metadata_metabolon <- make_metadata(readxl::read_excel(cdt, sheet = 3))
  metadata_se <- as.data.frame(colData(se))

  ids_to_match <- if (!is.null(sample_id_column)) {
    if (!(sample_id_column %in% colnames(metadata_se))) {
      stop("The specified sample_id_column does not exist in the colData of the SummarizedExperiment.")
    }
    metadata_se[[sample_id_column]]
  } else {
    colnames(se)
  }

  # Read assay
  assay <- assay(se)

  # Add PARENT_SAMPLE_NAME column (from documentation)
  assay_dataframe <- data.frame(t(assay))
  assay_dataframe$PARENT_SAMPLE_NAME <- rownames(metadata_metabolon)[match(ids_to_match, metadata_metabolon$PARENT_SAMPLE_NAME)]
  assay_dataframe <- assay_dataframe[, c("PARENT_SAMPLE_NAME", setdiff(names(assay_dataframe), "PARENT_SAMPLE_NAME"))]

  # Write the output to a file
  if (save_file) {
    message("Saving results to: ", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE)
    write.table(assay_dataframe, output_file, sep = ",", row.names = FALSE)
  }
  return(assay_dataframe)
}
