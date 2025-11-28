#' Map Gene Identifiers to Ensembl IDs (Generic Across Any OrgDb)
#'
#' @param keys_list Character vector of gene identifiers (any keytype supported by OrgDb).
#' @param matrix A data frame or matrix whose row names correspond to `keys_list`.
#' @param input_keytype Character string specifying the keytype of the input identifiers
#'   (must be one of `keytypes(orgdb)`).
#' @param orgdb An AnnotationDbi OrgDb object (e.g., org.Hs.eg.db::org.Hs.eg.db).
#' @param remove_unmapped Logical; if TRUE (default), remove rows that cannot be mapped.
#'
#' @return The matrix/data.frame with row names replaced by Ensembl Gene IDs.
#' @importFrom AnnotationDbi mapIds keytypes columns
map_to_ensembl <- function(keys_list,
                           matrix,
                           input_keytype,
                           orgdb,
                           remove_unmapped = TRUE) {

  # Validate OrgDb 
  if (!inherits(orgdb, "OrgDb")) {
    stop("`orgdb` must be an AnnotationDbi OrgDb object (e.g., org.Hs.eg.db::org.Hs.eg.db).")
  }

  # Validate keytype 
  valid_keytypes <- AnnotationDbi::keytypes(orgdb)
  if (!(input_keytype %in% valid_keytypes)) {
    stop("Invalid input_keytype. Must be one of: ",
         paste(valid_keytypes, collapse = ", "))
  }

  # Validate ENSEMBL availability 
  if (!("ENSEMBL" %in% AnnotationDbi::columns(orgdb))) {
    stop("The supplied OrgDb does not support ENSEMBL mappings.")
  }

  # Perform mapping 
  mapped <- AnnotationDbi::mapIds(
    orgdb,
    keys       = keys_list,
    keytype    = input_keytype,
    column     = "ENSEMBL",
    multiVals  = "first"
  )

  # Handle unmapped 
  unmapped <- is.na(mapped)

  if (any(unmapped) && remove_unmapped) {
    warning(sum(unmapped), " identifiers could not be mapped to ENSEMBL and were removed.")
    matrix <- matrix[!unmapped, , drop = FALSE]
    mapped <- mapped[!unmapped]
  }

  # Update rownames 
  rownames(matrix) <- make.names(mapped, unique = TRUE)

  return(matrix)
}
