# --- Tests for cdt_to_se ---
test_that("cdt_to_se returns a SummarizedExperiment with expected assays", {
  se <- cdt_to_se(
    cdt = cdt_path,
    output_file = NULL,
    save_file = FALSE,
    data_type = "batch_norm_imputed"
  )
  
  expect_s4_class(se, "SummarizedExperiment")
  expect_true("counts" %in% assayNames(se))
  expect_true(is.numeric(assay(se)))
  expect_true(all(!is.na(assay(se))))
})