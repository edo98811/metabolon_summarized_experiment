# Convert CDT File to SummarizedExperiment Object

This function reads a CDT (Client Data Table) Excel file and converts it
into a `SummarizedExperiment` object. The function allows selecting
specific data types and optionally saves the resulting object to a file.

## Usage

``` r
cdt_to_se(
  cdt,
  output_file = NULL,
  save_file = T,
  data_type = "batch_norm_imputed",
  rowdata_key = "INCHIKEY"
)
```

## Arguments

- cdt:

  Character. Path to the CDT Excel file. Must exist.

- output_file:

  Character. Path to save the resulting `SummarizedExperiment` object.
  If `NULL` and `save_file` is `TRUE`, a default filename will be
  generated.

- save_file:

  Logical. Whether to save the resulting `SummarizedExperiment` object
  to a file. Default is `TRUE`. If `FALSE`, the object will be returned.
  If a path is not provided, a default filename will be generated based
  on the current date and saved in the working directory.

- data_type:

  Character. The type of data to extract from the CDT file. Options are:

  - `"peak_area"`: Extracts peak area data (sheet 4).

  - `"batch_norm"`: Extracts batch normalized data (sheet 5).

  - `"batch_norm_imputed"`: Extracts batch normalized and imputed data
    (sheet 6).

  - `"mass_extracted"`: Extracts mass extracted data (sheet 7).

  - `"log_transformed"`: Extracts log-transformed data (sheet 8).

  Default is `"batch_norm_imputed"`.

- rowdata_key:

  Character. The column name in the row data (sheet 2) to use as row
  names.

## Value

A `SummarizedExperiment` object containing:

- `assays`: The assay data extracted from the CDT file.

- `rowData`: Metadata for the rows (e.g., compounds).

- `colData`: Metadata for the columns (e.g., samples).

## Details

The function performs the following steps:

1.  Validates the existence of the provided CDT.

2.  Reads metadata (sheet 3) and row data (sheet 2) from the CDT file.

3.  Extracts assay data based on the specified `data_type`.

4.  Cleans and formats the metadata and row data. The componds kept all
    the ones that have a non-NA value in the selected row data
    (annotation) column.

5.  Ensures consistency between metadata and assay data dimensions.

6.  Constructs a `SummarizedExperiment` object.

7.  Optionally saves the resulting object to a file if `save_file` is
    `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert a CDT file to a SummarizedExperiment object
se <- cdt_to_se(
  cdt = "path/to/cdt_file.xlsx",
  save_file = FALSE,
  data_type = "log_transformed"
)
} # }
```
