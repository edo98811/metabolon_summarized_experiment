# Convert Results to Metabolon Format

This function converts differential expression analysis results (from
limma, DESeq2, or custom formats) into a format compatible with
Metabolon multi-omics module.

## Usage

``` r
results_to_metabolon(
  table,
  format = "topTable",
  omics = "Transcriptomics",
  custom_colnames = NULL,
  output_file = NULL,
  input_features = "ensembl_id",
  organism = "Hs",
  save_file = T
)
```

## Arguments

- table:

  A data frame containing the results of differential expression
  analysis.

- format:

  A character string specifying the format of the input table. Options
  are "topTable" (limma), "DESeqResults" (DESeq2), "edgeR" (not yet
  supported), or "Custom". Default is "topTable".

- omics:

  A character string specifying the type of omics data. Default is
  "transcriptomics".

- custom_colnames:

  A character vector specifying custom column names for "Custom" format.
  Must include at least "p-value", "FDR", and "FC". Default is NULL.

- output_file:

  A character string specifying the output file path. If NULL and
  `save_file` is TRUE, a default file name will be generated. Default is
  NULL.

- input_features:

  A character string specifying the type of input features. Options are
  "gene_symbol", "uniprot_id", or "ensembl_id". Default is "ensembl_id".

- organism:

  A character string specifying the organism. Can be "Mm" or "Hs".
  Default is "Hs" (Homo sapiens).

- save_file:

  A logical value indicating whether to save the output to a file.
  Default is TRUE. If a path is not provided, a default filename will be
  generated based on the current date and saved in the working
  directory.

## Value

A data frame formatted for Metabolon, with columns: Feature, Omics,
p-value, FDR, and FC.

## Details

The function checks the correctness of the input format and ensures that
the required columns are present in the input table. It supports
workflows for "topTable" (limma), "DESeqResults" (DESeq2), and "Custom"
formats. The "edgeR" format is not yet supported.

## Examples

``` r
# Example usage with limma::topTable results
if (FALSE) { # \dontrun{
results <- limma::topTable(fit)
metabolon_results <- results_to_metabolon(results, format = "topTable")

# Example usage with DESeq2::results
results <- DESeq2::results(dds)
metabolon_results <- results_to_metabolon(results, format = "DESeqResults")

# Example usage with custom format
custom_results <- data.frame(
  p_value = c(0.01, 0.05),
  FDR = c(0.02, 0.06),
  FC = c(1.5, -1.2)
)
metabolon_results <- results_to_metabolon(custom_results, format = "Custom", 
                                          custom_colnames = c("p_value", "FDR", "FC"))
} # }
```
