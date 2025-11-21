
# --- Tests for se_to_metabolon ---
test_that("se_to_metabolon returns a data.frame with expected columns", {
  table1 <- se_to_metabolon(
    gse,
    cdt = cdt_path,
    input_features = "ensembl_id",
    organism = "Hs",
    save_file = FALSE,
    sample_id_column = "CLIENT_SAMPLE_ID"
  )
  
  expect_s3_class(table1, "data.frame")
})

# --- Tests for results_to_metabolon with DESeqResults ---
test_that("results_to_metabolon handles DESeqResults correctly", {
  res1 <- results_to_metabolon(
    res_macrophage_IFNg_vs_naive_dds,
    format = "DESeqResults",
    omics = "transcriptomics",
    custom_colnames = NULL,
    input_features = "ensembl_id",
    organism = "Hs",
    save_file = FALSE
  )
  
  expect_s3_class(res1, "data.frame")
})

# --- Tests for results_to_metabolon with limma topTable ---
test_that("results_to_metabolon handles topTable correctly", {
  res2 <- results_to_metabolon(
    res_macrophage_IFNg_vs_naive_limma,
    format = "topTable",
    omics = "transcriptomics",
    custom_colnames = NULL,
    input_features = "ensembl_id",
    organism = "Hs",
    save_file = FALSE
  )
  
  expect_s3_class(res2, "data.frame")
})

