# Convert SummarizedExperiment to Metabolon Format

This function converts a `SummarizedExperiment` object into a format
compatible with Metabolon data analysis.

## Usage

``` r
se_to_metabolon(
  se,
  cdt,
  output_file = NULL,
  save_file = TRUE,
  sample_id_column = NULL
)
```

## Arguments

- se:

  A `SummarizedExperiment` object containing the assay data.

- cdt:

  A character string specifying the path to the client data table path.

- output_file:

  A character string specifying the path to save the output file. If
  `NULL` and `save_file` is `TRUE`, a default filename will be generated
  and it will be saved in the working directory.

- save_file:

  A logical value indicating whether to save the output to a file. Needs
  to be csv! Default is `TRUE`.

- sample_id_column:

  A character string specifying the column name in the colData that
  contains the sample IDs matching the PARENT_SAMPLE_NAME names. Default
  is `NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
se <- SummarizedExperiment::SummarizedExperiment()
result <- se_to_metabolon(se, cdt = "cdt.xlsx", input_features = "gene_symbol", save_file = TRUE)
} # }
```
