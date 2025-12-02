# Convert Results to Metabolon Format.

This function converts a results table into a format compatible with
Metabolon data analysis.

## Usage

``` r
results_to_metabolon(
  table,
  format = "topTable",
  omics = "Transcriptomics",
  custom_colnames = NULL,
  output_file = NULL,
  save_file = TRUE
)
```

## Arguments

- table:

  A data frame containing the results to be converted.

- format:

  A character string specifying the format of the input results. Options
  are `"DESeqResults"` (default) or `"topTable"`.

- omics:

  A character string specifying the type of omics data. Options are
  `"Transcriptomics"` (default) or `"Proteomics"`.

- custom_colnames:

  An optional character vector specifying custom column names for the
  output. If `NULL`, default column names will be used.

- output_file:

  A character string specifying the path to save the output file. If
  `NULL` and `save_file` is `TRUE`, a default filename will be
  generated. Needs to be csv!

- save_file:

  A logical value indicating whether to save the output to a file.
  Default is `TRUE`.

## Value

A data frame formatted according to the standard required by Metabolon.

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert DESeqResults to Metabolon format
results_metabolon <- results_to_metabolon(
  table = deseq_results,
  format = "DESeqResults",
  omics = "Transcriptomics",
  input_features = "gene_symbol",
  save_file = TRUE
)
} # }
```
