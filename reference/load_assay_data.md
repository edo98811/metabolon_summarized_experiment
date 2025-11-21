# Load assay data based on the specified data type

Load assay data based on the specified data type

## Usage

``` r
load_assay_data(cdt, data_type)
```

## Arguments

- cdt:

  Character. Path to the CDT Excel file.

- data_type:

  Character. The type of data to extract. Options are: "peak_area",
  "batch_norm", "batch_norm_imputed", "mass_extracted",
  "log_transformed".

## Value

A data frame containing the assay data.
