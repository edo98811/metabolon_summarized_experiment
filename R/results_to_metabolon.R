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


  # Check if the input inherits from data frame
  if (!inherits(table, "data.frame")) stop("Input 'table' must inherit from class 'data.frame' (output of limma::topTable or DESeq2::results).")

  # Implements the workflows
  if(format == "topTable") {
    # Check if the input has the required columns
    required_columns <- c("logFC", "P.Value", "adj.P.Val")
    
    if (length(setdiff(required_columns, colnames(table))) > 0) {
      stop(paste("The input table is missing the following required columns:", paste(missing_columns, collapse = ", ")))
    }

    results_for_metabolon <- data.frame(
      Feature = rownames(table),
      Omics = omics,
      p-value = table$P.Value,
      FDR = table$adj.P.Val,
      FC = table$logFC,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } if(format == "DESeqResults") {
    # Check if the input has the required columns
    required_columns <- c("log2FoldChange", "pvalue", "padj")
    
    if (length(setdiff(required_columns, colnames(table))) > 0) {
      stop(paste("The input table is missing the following required columns:", paste(setdiff(required_columns, colnames(table)), collapse = ", ")))
    }

    results_for_metabolon <- data.frame(
      Feature = rownames(table),
      Omics = omics,
      p-value = table$pvalue,
      FDR = table$padj,
      FC = table$log2FoldChange,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } if(format == "Custom") {

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
      p-value = table$p-value,
      FDR = table$FDR,
      FC = table$FC,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  # Make the rownames conform to the metabolon format
  results_for_metabolon <- switch(input_features,
    "gene_symbol" = {
      rownames(results_for_metabolon) <- rownames_to_ensembl(results_for_metabolon, anns, from = "uniprot_id")
      return(results_for_metabolon)
    },
    "uniprot_id" = {
      rownames(results_for_metabolon) <- rownames_to_ensembl(results_for_metabolon, anns, from = "uniprot_id")
      return(results_for_metabolon)
    },
    "ensemble_id" = {
      return(results_for_metabolon)
    },
    stop("Invalid input_features. Choose from 'gene_symbol', 'uniprot_id', or 'ensemble_id'")
  )
  
  # Write the output to a file
  if (save_file) write.table(results_for_metabolon, output_file, sep = ",", row.names = TRUE)
  return(results_for_metabolon)
}