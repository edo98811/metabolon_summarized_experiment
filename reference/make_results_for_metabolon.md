# Prepare results for Metabolon upload

This function transforms results tables (from limma, DESeq2, or custom
formats) into a format compatible with Metabolon multi-omics module.

## Usage

``` r
make_results_for_metabolon(table, format, omics, custom_colnames = NULL)
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

## Value

A data frame formatted for Metabolon, with columns: Feature, Omics,
p-value, FDR, and FC.
