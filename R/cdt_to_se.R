cdt_to_se <- function(cdt_path,
                      output_file = NULL,
                      save_file = F,
                      data_type = "match_norm_imputed") {

  # Define the path to the Excel file
  if (!file.exists(cdt_path)) stop("The provided cdt_path does not exist.")
  if (is.null(output_file) && save_file) output_file <- paste0("se_from_metabolon_", Sys.Date(), ".csv")

  # Load the metadata sheet
  metadata <- read_xlsx(cdt_path, sheet = 3)
  rownames(metadata) <- metadata$PARENT_SAMPLE_NAME
  metadata <- metadata[, !colnames(metadata) %in% "PARENT_SAMPLE_NAME"]

  rowdata <- read_xlsx(cdt_path, sheet = 2)
  rownames(rowdata) <- rowdata$COMP_ID

  # Load the assay data
  assay_data <- switch(data_type,
    "peak_area" = {
      t(read_xlsx(cdt_path, sheet = 4, rowNames = TRUE))
    },
    "batch_norm" = {
      t(read_xlsx(cdt_path, sheet = 5, rowNames = TRUE))
    },
    "match_norm_imputed" = {
      t(read_xlsx(cdt_path, sheet = 6, rowNames = TRUE))
    },
    "mass_extracted" = {
      t(read_xlsx(cdt_path, sheet = 7, rowNames = TRUE))
    },
    "log_transformed" = {
      t(read_xlsx(cdt_path, sheet = 8, rowNames = TRUE))
    },
  )

  # Check that the metadata dataframes are correct
  stopifnot(all(rownames(rowdata) == rownames(assay_data)))
  stopifnot(all(rownames(metadata) == colnames(assay_data)))

  # Create the SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = list(counts = assay_data),
    rowData = rowdata,
    colData = metadata
  )

  # Save the SummarizedExperiment object
  if (save_file) saveRDS(se, output_file)
  return(se)
  }
