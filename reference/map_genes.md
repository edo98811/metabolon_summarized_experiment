# Map Gene Identifiers to Ensembl IDs This function maps gene identifiers (gene symbols or UniProt IDs) to Ensembl Gene IDs using the appropriate organism database.

Map Gene Identifiers to Ensembl IDs

This function maps gene identifiers (gene symbols or UniProt IDs) to
Ensembl Gene IDs using the appropriate organism database.

## Usage

``` r
map_genes(keys_list, matrix, input_features, organism)
```

## Arguments

- keys_list:

  A character vector of gene identifiers (gene symbols, UniProt IDs, or
  Ensembl IDs).

- matrix:

  A data frame or matrix containing the assay data, with row names
  corresponding to the identifiers in `keys_list`.

- input_features:

  A character string specifying the type of input features. Options are
  `"gene_symbol"`, `"uniprot_id"`, or `"ensembl_id"`.

- organism:

  A character string specifying the organism. Can be `"Mm"` (Mus
  musculus) or `"Hs"` (Homo sapiens). Default is `"Hs"`.

## Value

matrix with row names updated to Ensembl IDs, suitable for downstream
analysis.
