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
#' @param use_adj_pvalue A logical value indicating whether to use adjusted p-values. Default is TRUE.
#' @param output_file A character string specifying the output file path. If NULL and `save_file` is TRUE, 
#' a default file name will be generated. Default is NULL.
#' @param input_features A character string specifying the type of input features. 
#' Options are "gene_symbol", "uniprot_id", or "ensemble_id". Default is "ensemble_id".
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
                                omics = "transcriptomics",
                                custom_colnames = NULL,
                                use_adj_pvalue = TRUE, 
                                output_file = NULL, 
                                input_features = "ensemble_id",
                                organism = "Hs",
                                save_file = F) {
  
  # Check correctness of input
  if (!format %in% c("topTable", "DESeqResults", "edgeR", "Custom")) {
    stop("Invalid format. Choose from 'topTable', 'DESeqResults', 'Custom' or 'edgeR'.")
  }
  if (is.null(output_file) && save_file) output_file <- paste0("results_to_metabolon_", Sys.Date(), ".csv")
  if (format == "edgeR") stop("edgeR format is not supported yet.")




  # Implements the workflows
  if(format == "topTable") {

    # Check if the input inherits from data frame
    if (!inherits(table, "data.frame")) stop("Input 'table' must inherit from class 'data.frame' (or output of limma::topTable).")

    # Check if the input has the required columns
    required_columns <- c("logFC", "P.Value", "adj.P.Val")
    
    if (length(setdiff(required_columns, colnames(table))) > 0) {
      stop(paste("The input table is missing the following required columns:", paste(missing_columns, collapse = ", ")))
    }

    results_for_metabolon <- data.frame(
      Feature = rownames(table),
      Omics = omics,
      "p-value" = table$P.Value,
      FDR = table$adj.P.Val,
      FC = table$logFC,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else if(format == "DESeqResults") {

    # Check if the input inherits from data frame
    if (!inherits(table, "DESeqResults")) stop("Input 'table' must inherit from class DESeqResults.")
    table <- as.data.frame(table)

    # Check if the input has the required columns
    required_columns <- c("log2FoldChange", "pvalue", "padj")
    
    if (length(setdiff(required_columns, colnames(table))) > 0) {
      stop(paste("The input table is missing the following required columns:", paste(setdiff(required_columns, colnames(table)), collapse = ", ")))
    }

    results_for_metabolon <- data.frame(
      Feature = rownames(table),
      Omics = omics,
      "p-value" = table$pvalue,
      FDR = table$padj,
      FC = table$log2FoldChange,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else if(format == "Custom") {

    # # Check if the input has the required columns
    # if (!is.character(custom_colnames) || length(custom_colnames) < 3) {
    #   stop("'custom_colnames' must be a character vector of at least 3 elements.")
    # }
    # required_columns <- names(custom_colnames)
    required_columns <- c("p-value", "FDR", "FC")

    if (length(setdiff(colnames()(colnames), colnames(table))) > 0) {
      stop(paste("The input table is missing the following required columns:", paste(setdiff(required_columns, colnames(table)), collapse = ", ")))
    }

    results_for_metabolon <- data.frame(
      Feature = rownames(table),
      Omics = omics,
      "p-value" = table$p-value,
      FDR = table$FDR,
      FC = table$FC,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  # Make the rownames conform to the metabolon format and transpose the assay
  results_for_metabolon <- switch(input_features,
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
        results_for_metabolon <- results_for_metabolon[!unmapped, , drop = FALSE]
        uniprot_ids <- uniprot_ids[!unmapped]
      }
      colnames(results_for_metabolon) <- make.names(uniprot_ids, unique = TRUE)
      results_for_metabolon
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
        results_for_metabolon <- results_for_metabolon[!unmapped, , drop = FALSE]
        uniprot_ids <- uniprot_ids[!unmapped]
      }
      colnames(results_for_metabolon) <- make.names(uniprot_ids, unique = TRUE)
      results_for_metabolon
    },
    "ensembl_id" = {
      # No mapping needed; retain original Ensembl IDs
      results_for_metabolon
    },

    stop("Invalid input_features. Choose from 'gene_symbol', 'uniprot_id', or 'ensembl_id'")
  )
  
  # Write the output to a file
  if (save_file) {
    dir.create(dirname(output_file), showWarnings = FALSE)
    write.table(results_for_metabolon, output_file, sep = ",", row.names = TRUE)
  }
  return(results_for_metabolon)
}