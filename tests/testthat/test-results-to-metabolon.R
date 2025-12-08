# --- Tests for se_to_metabolon ---
test_that("se_to_metabolon returns a data.frame with expected columns", {
  table1 <- se_to_metabolon(
    gse,
    cdt = cdt_path,
    save_file = FALSE
  )

  expect_s3_class(table1, "data.frame")
  expect_true("PARENT_SAMPLE_NAME" %in% colnames(table1))
  expect_true(all(sapply(table1[2:ncol(table1)], is.numeric)))
})

# --- Tests for results_to_metabolon with DESeqResults ---
test_that("results_to_metabolon handles DESeqResults correctly", {
  res1 <- results_to_metabolon(
    res_macrophage_IFNg_vs_naive_dds,
    format = "DESeqResults",
    omics = "Transcriptomics",
    custom_colnames = NULL,
    save_file = FALSE
  )

  expect_s3_class(res1, "data.frame")
  expect_equal(colnames(res1), c("Feature", "Omics", "p-value", "FDR", "FC"))
  expect_true(all(sapply(res1[c("p-value", "FDR", "FC")], is.numeric)))
  expect_true(all(sapply(res1[c("Feature", "Omics")], is.character)))
})

# --- Tests for results_to_metabolon with limma topTable ---
test_that("results_to_metabolon handles topTable correctly", {
  res2 <- results_to_metabolon(
    res_macrophage_IFNg_vs_naive_limma,
    format = "topTable",
    omics = "Transcriptomics",
    custom_colnames = NULL,
    save_file = FALSE
  )

  expect_s3_class(res2, "data.frame")
  expect_equal(colnames(res2), c("Feature", "Omics", "p-value", "FDR", "FC"))
  expect_true(all(sapply(res2[c("p-value", "FDR", "FC")], is.numeric)))
  expect_true(all(sapply(res2[c("Feature", "Omics")], is.character)))
})
