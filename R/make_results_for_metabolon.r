#' Prepare results for Metabolon upload
#'
#' This function transforms results tables (from limma, DESeq2, or custom formats)
#' into a format compatible with Metabolon multi-omics module.
#'
#' @param table A data frame containing the results of differential expression analysis.
#' @param format A character string specifying the format of the input table.
#' Options are "topTable" (limma), "DESeqResults" (DESeq2), "edgeR" (not yet supported),
#' or "Custom". Default is "topTable".
#' @param omics A character string specifying the type of omics data. Default is "transcriptomics".
#' @param custom_colnames A character vector specifying custom column names for "Custom" format.
#' Must include at least "p-value", "FDR", and "FC". Default is NULL.
#' @return A data frame formatted for Metabolon, with columns: Feature, Omics, p-value, FDR, and FC.
make_results_for_metabolon <- function(table, format, omics, custom_colnames = NULL) {
  format <- match.arg(format, c("topTable", "DESeqResults", "Custom"))

  # Check if the table has the required columns
  table <- check_table(table, format)
  table <- check_columns(table, format)

  # Implement the 
  results_for_metabolon <- switch(format,
    topTable = {
      data.frame(
        Feature = rownames(table),
        Omics = omics,
        "p-value" = table$P.Value,
        FDR = table$adj.P.Val,
        FC = table$logFC,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    },
    DESeqResults = {
      data.frame(
        Feature = rownames(table),
        Omics = omics,
        "p-value" = table$pvalue,
        FDR = table$padj,
        FC = table$log2FoldChange,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    },
    Custom = {
      data.frame(
        Feature = rownames(table),
        Omics = omics,
        "p-value" = table[["p-value"]],
        FDR = table$FDR,
        FC = table$FC,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }
  )

  return(results_for_metabolon)
}

check_columns <- function(table, format) {
  # Define the required columns for each format
  required_columns <- switch(format,
    topTable = c("logFC", "P.Value", "adj.P.Val"),
    DESeqResults = c("log2FoldChange", "pvalue", "padj"),
    Custom = c("p-value", "FDR", "FC")
  )

  missing_columns <- setdiff(required_columns, colnames(table))
  if (length(missing_columns) > 0) {
    stop(paste("The input table is missing the following required columns:", paste(missing_columns, collapse = ", ")))
  }

  return(table)
}


check_table <- function(table, format) {
  table <- switch(format,
    topTable = {
      if (!inherits(table, "data.frame")) stop("Input 'table' must inherit from class 'data.frame' (or output of limma::topTable).")
      table
    },
    DESeqResults = {
      if (!(inherits(table, "DESeqResults") || inherits(table, "data.frame"))) stop("Input 'table' must inherit from class DESeqResults.")
      as.data.frame(table)
    },
    Custom = {
      if (!inherits(table, "data.frame")) stop("Input 'table' must inherit from class 'data.frame'.")
      table
    }
  )

  return(table)
}
