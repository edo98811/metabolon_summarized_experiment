create_annotations <- function(se,
                               organism = "Hs",
                               annotation_type = "ensembl") {

  if (!inherits(se, "SummarizedExperiment")) stop("se does not inherit from SummarizedExperiment object")

  # "www" → Main server (https://www.ensembl.org)
  # "useast" → US East (https://useast.ensembl.org)
  # "uswest" → US West (https://uswest.ensembl.org)
  # "asia" → Asia (https://asia.ensembl.org)
  if (organism == "Hs") {
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
  } else if (organism == "Mm") {
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "https://www.ensembl.org")
  } else {
    stop("Invalid organism")
  }

  # https://www.rdocumentation.org/packages/biomaRt/versions/2.28.0/topics/getBM
  anns <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "uniprotswissprot", "description"),
    filters = switch(annotation_type,
      "ensembl_id" = "ensembl_gene_id",
      "gene_symbol" = "external_gene_name",
      "uniprot_id" = "uniprotswissprot",
      stop("Invalid annotation_type Choose from 'ensembl_id', 'gene_symbol', 'uniprot_id'")
    ),
    values = rownames(se),
    mart = mart
  )

  colnames(anns) <- c("ensembl_gene_id", "gene_symbol", "uniprot_id", "description")

  anns <- anns[match(rownames(se), anns$ensembl_gene_id), ]
  return(anns)
}
