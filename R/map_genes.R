  map_genes <- function(keys_list, matrix, input_features) {

    switch(input_features,
      "gene_symbol" = {
        gene_ids <- mapIds(
          org.Hs.eg.db,
          keys = keys_list,
          column = "ENSEMBL",
          keytype = "SYMBOL",
          multiVals = "first"
        )
        unmapped <- is.na(gene_ids)
        if (any(unmapped)) {
          warning("Some rows could not be mapped to UniProt IDs and will be removed.")
          matrix <- matrix[!unmapped, , drop = FALSE]
          gene_ids <- gene_ids[!unmapped]
        }
        colnames(matrix) <- make.names(gene_ids, unique = TRUE)
        matrix
      },
      
      "uniprot_id" = {
        uniprot_ids <- mapIds(
          org.Hs.eg.db,
          keys = keys_list,
          column = "ENSEMBL",
          keytype = "UNIPROT",
          multiVals = "first"
        )
        unmapped <- is.na(uniprot_ids)
        if (any(unmapped)) {
          warning("Some rows could not be mapped to UniProt IDs and will be removed.")
          matrix <- matrix[!unmapped, , drop = FALSE]
          uniprot_ids <- uniprot_ids[!unmapped]
        }
        colnames(matrix) <- make.names(uniprot_ids, unique = TRUE)
        matrix
      },
      
      "ensembl_id" = {
        # No mapping needed; retain original Ensembl IDs
        matrix
      },
      
      stop("Invalid input_features. Choose from 'gene_symbol', 'uniprot_id', or 'ensembl_id'")
    )
  }
