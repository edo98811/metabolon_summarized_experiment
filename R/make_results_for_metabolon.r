make_results_for_metabolon <- function(table, format, omics, custom_colnames = NULL) {
  # Implements the workflows
  if (format == "topTable") {
    # Check if the input inherits from data frame
    if (!inherits(table, "data.frame")) stop("Input 'table' must inherit from class 'data.frame' (or output of limma::topTable).")

    # Check if the input has the required columns
    required_columns <- c("logFC", "P.Value", "adj.P.Val")
    missing_columns <- setdiff(required_columns, colnames(table))
    if (length(missing_columns) > 0) {
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
  } else if (format == "DESeqResults") {
    # Check if the input inherits from data frame
    if (!(inherits(table, "DESeqResults") || inherits(table, "data.frame"))) stop("Input 'table' must inherit from class DESeqResults.")
    table <- as.data.frame(table)

    # Check if the input has the required columns
    required_columns <- c("log2FoldChange", "pvalue", "padj")
    missing_columns <- setdiff(required_columns, colnames(table))
    if (length(missing_columns) > 0) {
      stop(paste("The input table is missing the following required columns:", paste(missing_columns, collapse = ", ")))
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
  } else if (format == "Custom") {
    # Check if the input has the required columns
    required_columns <- c("p-value", "FDR", "FC")
    missing_columns <- setdiff(required_columns, colnames(table))
    if (length(missing_columns) > 0) {
      stop(paste("The input table is missing the following required columns:", paste(missing_columns, collapse = ", ")))
    }

    results_for_metabolon <- data.frame(
      Feature = rownames(table),
      Omics = omics,
      "p-value" = table[["p-value"]],
      FDR = table$FDR,
      FC = table$FC,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else {
    stop("Invalid format specified. Must be one of 'topTable', 'DESeqResults', or 'Custom'.")
  }

  return(results_for_metabolon)
}
