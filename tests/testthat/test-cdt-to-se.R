test_that("cdt_to_se works", {
  se <- cdt_to_se(
    cdt = cdt_path,
    output_file = "test/cdt_to_se_bni.rds",
    save_file = F,
    data_type = "batch_norm_imputed"
  )
  expect_s4_class(se, "SummarizedExperiment")
  # expect_equal(nrow(se), 100)
  # expect_equal(ncol(se), 6)
})
