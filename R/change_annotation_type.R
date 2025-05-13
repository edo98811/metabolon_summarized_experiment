
rownames_to_ensembl <- function(objective_table, anns, from = "uniprot_id") {

    if (!"ensembl_gene_id" %in% colnames(assay_transposed)) {
        objective_table <- add_gene_info(objective_table, anns)
    }

    objective_table <- objective_table[, c("ensembl_gene_id", setdiff(names(objective_table), "ensembl_gene_id"))]
    na_genes <- objective_table$uniprot_id[is.na(objective_table$ensembl_gene_id)]

    # warning if NA or duplicated gene names
    if (length(na_genes) > 0) {
        warning("There are NA gene names (will be deleted): ", paste(na_genes, collapse = ", "))
    }

    # Warning if duplicated gene names
    duplicated_genes <- objective_table$uniprot_id[duplicated(objective_table$ensembl_gene_id)]
    if (length(duplicated_genes) > 0) {
        warning("There are duplicated gene name (will be deleted): ", paste(duplicated_genes, collapse = ", "))
    }

    # Remove duplicated and NA gene names
    objective_table <- objective_table[na_genes & duplicated_genes, ]

    # Assign the ensembl gene id as rownames and remove the column
    rownames(objective_table) <- objective_table$ensembl_gene_id
    objective_table <- objective_table[, -which(names(objective_table) == "ensembl_gene_id")]

    return(rownames(objective_table))
}


add_gene_info <- function(objective_table, anns) {

    if (!"uniprot_id" %in% colnames(objective_table)) {
        objective_table$uniprot_id <- rownames(objective_table)
    }  
    
    objective_table$gene_symbol <- anns[match(objective_table$uniprot_id, anns$unprot_id), "gene_symbol"]
    objective_table$ensembl_gene_id <- anns[match(objective_table$uniprot_id, anns$unprot_id), "ensembl_gene_id"]

    return(objective_table)
}

# to_prot_id <- function(objective_table, anns) {


#     objective_table$prot_ID <- rownames(anns[match(rownames(objective_table), anns$gene_symbol), ]) # rownames are protID (magari da cambiare per interoperability)
#     objective_table <- objective_table[, c("prot_ID", setdiff(names(objective_table), "prot_ID"))]
#     objective_table <- objective_table[!duplicated(objective_table$prot_ID) & !is.na(objective_table$prot_ID), ]
#     rownames(objective_table) <- objective_table$prot_ID
#     objective_table <- objective_table[, -which(names(objective_table) == "prot_ID")]

#     return(objective_table)
# }
