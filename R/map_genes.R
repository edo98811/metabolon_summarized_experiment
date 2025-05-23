map_genes <- function(keys_list, matrix, input_features, organism) {

  # Choose the correct OrgDb object using if/else (not ifelse)
  orgdb <- if (organism == "Mm") {
    org.Mm.eg.db::org.Mm.eg.db
  } else {
    org.Hs.eg.db::org.Hs.eg.db
  }

  switch(input_features,
    "gene_symbol" = {
      gene_ids <- AnnotationDbi::mapIds(
        orgdb,
        keys = keys_list,
        column = "ENSEMBL",
        keytype = "SYMBOL",
        multiVals = "first"
      )
      unmapped <- is.na(gene_ids)
      if (any(unmapped)) {
        warning("Some rows could not be mapped to Ensembl IDs and will be removed.")
        matrix <- matrix[!unmapped, , drop = FALSE]
        gene_ids <- gene_ids[!unmapped]
      }
      rownames(matrix) <- make.names(gene_ids, unique = TRUE)
      matrix
    },
    "uniprot_id" = {
      uniprot_ids <- AnnotationDbi::mapIds(
        orgdb,
        keys = keys_list,
        column = "ENSEMBL",
        keytype = "UNIPROT",
        multiVals = "first"
      )
      unmapped <- is.na(uniprot_ids)
      if (any(unmapped)) {
        warning("Some rows could not be mapped to Ensembl IDs and will be removed.")
        matrix <- matrix[!unmapped, , drop = FALSE]
        uniprot_ids <- uniprot_ids[!unmapped]
      }
      rownames(matrix) <- make.names(uniprot_ids, unique = TRUE)
      matrix
    },
    
    "ensembl_id" = {
      tryCatch({
      # No mapping needed; retain original Ensembl IDs
      rownames(matrix) <- make.names(keys_list, unique = TRUE)
      }, error = function(e) {
      stop("An error occurred while creating the results table IDs: ", e$message)
      })
      matrix
    },
    
    stop("Invalid input_features. Choose from 'gene_symbol', 'uniprot_id', or 'ensembl_id'")
  )
}
