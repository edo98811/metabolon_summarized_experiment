#' @title Convert Results to Metabolon Format
#' @description This function converts differential expression analysis results 
#' (from limma, DESeq2, or custom formats) into a format compatible with Metabolon.
#' @param table A data frame containing the results of differential expression analysis.
#' @param format A character string specifying the format of the input table. 
#' Options are "topTable" (limma), "DESeqResults" (DESeq2), "edgeR" (not yet supported), 
#' or "Custom". Default is "topTable".
#' @param omics A character string specifying the type of omics data. Default is "transcriptomics".
#' @param custom_colnames A character vector specifying custom column names for "Custom" format. 
#' Must include at least "p-value", "FDR", and "FC". Default is NULL.
#' @param output_file A character string specifying the output file path. If NULL and `save_file` is TRUE, 
#' a default file name will be generated. Default is NULL.
#' @param input_features A character string specifying the type of input features. 
#' Options are "gene_symbol", "uniprot_id", or "ensembl_id". Default is "ensembl_id".
#' @param organism A character string specifying the organism. Default is "Hs" (Homo sapiens).
#' @param save_file A logical value indicating whether to save the output to a file. Default is FALSE.
#' @return A data frame formatted for Metabolon, with columns: Feature, Omics, p-value, FDR, and FC.
#' @details The function checks the correctness of the input format and ensures that the required 
#' columns are present in the input table. It supports workflows for "topTable" (limma), 
#' "DESeqResults" (DESeq2), and "Custom" formats. The "edgeR" format is not yet supported.
#' @examples
#' # Example usage with limma::topTable results
#' results <- limma::topTable(fit)
#' metabolon_results <- results_to_metabolon(results, format = "topTable")
#' 
#' # Example usage with DESeq2::results
#' results <- DESeq2::results(dds)
#' metabolon_results <- results_to_metabolon(results, format = "DESeqResults")
#' 
#' # Example usage with custom format
#' custom_results <- data.frame(
#'   p_value = c(0.01, 0.05),
#'   FDR = c(0.02, 0.06),
#'   FC = c(1.5, -1.2)
#' )
#' metabolon_results <- results_to_metabolon(custom_results, format = "Custom", 
#'                                           custom_colnames = c("p_value", "FDR", "FC"))
#' @export
results_to_metabolon <- function(table,
                                format = "topTable",
                                omics = "Transcriptomics",
                                custom_colnames = NULL,
#                                 use_adj_pvalue = TRUE, 
                                output_file = NULL, 
                                input_features = "ensembl_id",
                                organism = "Hs",
                                save_file = T) {
  
  # Check correctness of input
  if (!format %in% c("topTable", "DESeqResults", "edgeR", "Custom")) {
    stop("Invalid format. Choose from 'topTable', 'DESeqResults', 'Custom' or 'edgeR'.")
  }
  # Check if organism is valid
  if (!organism %in% c("Mm", "Hs")) {
    stop("Invalid organism. Choose either 'Mm' (Mus musculus) or 'Hs' (Homo sapiens).")
  }
  if (is.null(output_file) && save_file) output_file <- paste0("results_to_metabolon_", Sys.Date(), ".csv")
  if (format == "edgeR") stop("edgeR format is not supported yet.")

  # Create the results for metabolon
  results_for_metabolon <- make_results_for_metabolon(table, format, omics, custom_colnames)

  # Make the rownames conform to the metabolon format and transpose the assay
  results_for_metabolon <- map_genes(results_for_metabolon$Feature, results_for_metabolon, input_features, organism)
  results_for_metabolon$Feature <- rownames(results_for_metabolon)

  # Write the output to a file
  if (save_file) {
    message("Saving results to: ", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE)
    write.table(results_for_metabolon, output_file, sep = ",", row.names = FALSE)
  } else return(results_for_metabolon)
}