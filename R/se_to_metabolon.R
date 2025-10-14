#' Convert SummarizedExperiment to Metabolon Format
#'
#' This function converts a `SummarizedExperiment` object into a format compatible
#' with Metabolon data analysis. It allows for customization of input features
#' and provides options to save the output to a file.
#'
#' @param se A `SummarizedExperiment` (or derived) object. The input data to be converted.
#' @param cdt The client data table from metabolon.
#' @param input_features A character string specifying the type of input features.
#'   Options are `"ensembl_id"` (default), `"gene_symbol"`, or `"uniprot_id"`.
#' @param output_file A character string specifying the path to save the output file.
#'   If `NULL` and `save_file` is `TRUE`, a default filename will be generated.
#' @param organism A character string specifying the organism. Default is `"Hs"`
#'   (Homo sapiens), Options are: "Mm".
#' @param save_file A logical value indicating whether to save the output to a file.
#'   Default is `TRUE`.
#'
#' @return A transposed assay matrix with row names and column names formatted according to the standard required by metabolon.
#'
#' @details
#' - The function checks if the input object is a `SummarizedExperiment`.
#' - Metadata is read from the CDT.
#' - Annotations are created based on the specified organism and input features.
#' - The assay matrix is transposed and formatted according to the metabolon standard.
#' - If `save_file` is `TRUE`, the output is written to a CSV file.
#'
#' @examples
#' \dontrun{
#' se <- SummarizedExperiment::SummarizedExperiment()
#' result <- se_to_metabolon(se, cdt = "cdt.xlsx", input_features = "gene_symbol", save_file = TRUE)
#' }
#'
#' @export
se_to_metabolon <- function(se,
                            cdt,
                            input_features = "ensembl_id",
                            output_file = NULL,
                            organism = "Hs",
                            save_file = T) {
  # Check correctness of input
  if (!inherits(se, "SummarizedExperiment")) stop("The input object is not a SummarizedExperiment.")
  if (is.null(output_file) && save_file) output_file <- paste0("se_to_metabolon_", Sys.Date(), ".csv")
  # Check if organism is valid
  if (!organism %in% c("Mm", "Hs")) {
    stop("Invalid organism. Choose either 'Mm' (Mus musculus) or 'Hs' (Homo sapiens).")
  }

  # Read metadata
  metadata <- openxlsx2::read_xlsx(cdt, sheet = 3, rowNames = TRUE)

  # Get transposed matrix
  assay <- as.data.frame(assay(se))
  assay <- map_genes(rownames(se), assay, input_features, organism)

  # Add PARENT_SAMPLE_NAME column
  assay_transposed <- t(assay)
  assay_transposed$PARENT_SAMPLE_NAME <- rownames(metadata)[match(colnames(assay), metadata$CLIENT_SAMPLE_ID)]

  # Write the output to a file
  if (save_file) {
    message("Saving results to: ", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE)
    write.table(assay_transposed, output_file, sep = ",", row.names = TRUE)
  }
  return(assay_transposed)
}
