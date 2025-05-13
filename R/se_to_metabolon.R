se_to_metabolon <- function(se, 
                            cdt, 
                            input_features = "ensemble_id",
                            output_file = NULL, 
                            organism = "Hs",
                            save_file = F) {

  # Check correctness of input
  if (!inherits(se, "SummarizedExperiment")) stop("The input object is not a SummarizedExperiment.")
  if (is.null(output_file) && save_file) output_file <- paste0("se_to_metabolon_", Sys.Date(), ".csv")

  # Read metadata
  metadata <- openxlsx2::read_xlsx(metadata_file, sheet = 3, rowNames = TRUE)
  
  # Add PARENT_SAMPLE_NAME column
  assay_transposed$PARENT_SAMPLE_NAME <- rownames(metadata)[match(rownames(assay_transposed), metadata$CLIENT_SAMPLE_ID)]
  
  anns <- create_annotations(se, organism = organism, annotation_type = input_features)

  # Make the rownames conform to the metabolon format and transpose the assay
  assay_transposed <- switch(input_features,
    "gene_symbol" = {
      assay_transposed <- t(assay(se))
      colnames(assay_transposed) <- rownames_to_ensembl(rowData(se), anns, from = "uniprot_id")
      return(assay_transposed)
    },
    "uniprot_id" = {
      assay_transposed <- t(assay(se))
      colnames(assay_transposed) <- rownames_to_ensembl(rowData(se), anns, from = "uniprot_id")
      return(assay_transposed)
    },
    "ensemble_id" = {
      return(t(assay(se)))
    },
    stop("Invalid input_features. Choose from 'gene_symbol', 'uniprot_id', or 'ensemble_id'")
  )

  # Write the output to a file
  if (save_file) write.table(assay_transposed, output_file, sep = ",", row.names = TRUE)

  return(assay_transposed)
}

