# Convert SummarizedExperiment to Metabolon Format

This function converts a `SummarizedExperiment` object into a format
compatible with Metabolon data analysis.

## Usage

``` r
se_to_metabolon(
  se,
  cdt,
  input_features = "ensembl_id",
  output_file = NULL,
  organism = "Hs",
  save_file = TRUE,
  sample_id_column = "CLIENT_SAMPLE_ID"
)
```

## Arguments

- se:

  A `SummarizedExperiment` (or derived) object. The input data to be
  converted.

- cdt:

  The client data table from metabolon.

- input_features:

  A character string specifying the type of input features. Options are
  `"ensembl_id"` (default), `"gene_symbol"`, or `"uniprot_id"`.

- output_file:

  A character string specifying the path to save the output file. If
  `NULL` and `save_file` is `TRUE`, a default filename will be
  generated.

- organism:

  A character string specifying the organism. Default is `"Hs"` (Homo
  sapiens), Options are: "Mm".

- save_file:

  A logical vawrlue indicating whether to save the output to a file.
  Default is `TRUE`.

- sample_id_column:

  The column in the metadata that needs to be matched to the sample id
  in the gse.

## Value

A transposed assay matrix with row names and column names formatted
according to the standard required by metabolon.

## Details

- The function checks if the input object is a `SummarizedExperiment`.

- Metadata is read from the CDT.

- Annotations are created based on the specified organism and input
  features.

- The assay matrix is transposed and formatted according to the
  metabolon standard.

- If `save_file` is `TRUE`, the output is written to a CSV file.

## Examples

``` r
if (FALSE) { # \dontrun{
se <- SummarizedExperiment::SummarizedExperiment()
result <- se_to_metabolon(se, cdt = "cdt.xlsx", input_features = "gene_symbol", save_file = TRUE)
} # }
```
