testthat("cdt_to_se works", {
  se <- cdt_to_se(
    system.file("extdata", "macrophage_example.cdt", package = "metabolonize"),
    input_features = "gene_symbol",
    organism = "Mm",
    save_file = FALSE
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(nrow(se), 100)
  expect_equal(ncol(se), 6)
})

