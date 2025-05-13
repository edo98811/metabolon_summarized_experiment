#' Convert SummarizedExperiment to Metabolon Format
#'
#' This function converts a `SummarizedExperiment` object into a format compatible 
#' with Metabolon data analysis. It allows for customization of input features 
#' and provides options to save the output to a file.
#'
#' @param se A `SummarizedExperiment` object. The input data to be converted.
#' @param cdt A parameter (not used in the current implementation).
#' @param input_features A character string specifying the type of input features. 
#'   Options are `"ensemble_id"` (default), `"gene_symbol"`, or `"uniprot_id"`.
#' @param output_file A character string specifying the path to save the output file. 
#'   If `NULL` and `save_file` is `TRUE`, a default filename will be generated.
#' @param organism A character string specifying the organism. Default is `"Hs"` 
#'   (Homo sapiens).
#' @param save_file A logical value indicating whether to save the output to a file. 
#'   Default is `FALSE`.
#'
#' @return A transposed assay matrix with row names and column names formatted 
#'   according to the specified `input_features`.
#'
#' @details
#' - The function checks if the input object is a `SummarizedExperiment`.
#' - Metadata is read from an external file (assumed to be in Excel format).
#' - Annotations are created based on the specified organism and input features.
#' - The assay matrix is transposed and formatted according to the input features.
#' - If `save_file` is `TRUE`, the output is written to a CSV file.
#'
#' @examples
#' \dontrun{
#' se <- SummarizedExperiment::SummarizedExperiment()
#' result <- se_to_metabolon(se, cdt = NULL, input_features = "gene_symbol", save_file = TRUE)
#' }
#'
#' @export
se_to_metabolon <- function(se, 
                            cdt, 
                            input_features = "ensembl_id",
                            output_file = NULL, 
                            organism = "Hs",
                            save_file = F) {

  # Check correctness of input
  if (!inherits(se, "SummarizedExperiment")) stop("The input object is not a SummarizedExperiment.")
  if (is.null(output_file) && save_file) output_file <- paste0("se_to_metabolon_", Sys.Date(), ".csv")

  # Read metadata
  metadata <- openxlsx2::read_xlsx(cdt, sheet = 3, rowNames = TRUE)

  # Get transposed matrix
  assay_transposed <- as.data.frame(t(assay(se)))

  # Add PARENT_SAMPLE_NAME column
  assay_transposed$PARENT_SAMPLE_NAME <- rownames(metadata)[match(rownames(assay_transposed), metadata$CLIENT_SAMPLE_ID)]
  
  # anns <- create_annotations(se, organism = organism, annotation_type = input_features)

  # Make the rownames conform to the metabolon format and transpose the assay
  assay_transposed <- switch(input_features,
    "gene_symbol" = {
      gene_ids <- mapIds(
        org.Hs.eg.db,
        keys = rownames(se),
        column = "ENSEMBL",
        keytype = "SYMBOL",
        multiVals = "first"
      )
      unmapped <- is.na(uniprot_ids)
      if (any(unmapped)) {
        warning("Some rows could not be mapped to UniProt IDs and will be removed.")
        assay_transposed <- assay_transposed[!unmapped, , drop = FALSE]
        uniprot_ids <- uniprot_ids[!unmapped]
      }
      colnames(assay_transposed) <- make.names(uniprot_ids, unique = TRUE)
      assay_transposed
    },
    
    "uniprot_id" = {
      uniprot_ids <- mapIds(
      org.Hs.eg.db,
      keys = rownames(se),
      column = "ENSEMBL",
      keytype = "UNIPROT",
      multiVals = "first"
      )
      unmapped <- is.na(uniprot_ids)
      if (any(unmapped)) {
        warning("Some rows could not be mapped to UniProt IDs and will be removed.")
        assay_transposed <- assay_transposed[!unmapped, , drop = FALSE]
        uniprot_ids <- uniprot_ids[!unmapped]
      }
      colnames(assay_transposed) <- make.names(uniprot_ids, unique = TRUE)
      assay_transposed
    },
    "ensembl_id" = {
      # No mapping needed; retain original Ensembl IDs
      assay_transposed
    },

    stop("Invalid input_features. Choose from 'gene_symbol', 'uniprot_id', or 'ensembl_id'")
  )

  # Write the output to a file
  if (save_file) {
    dir.create(dirname(output_file), showWarnings = FALSE)
    write.table(assay_transposed, output_file, sep = ",", row.names = TRUE)
  }

  return(assay_transposed)
}

