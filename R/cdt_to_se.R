#' Convert CDT File to SummarizedExperiment Object
#'
#' This function reads a CDT (Comprehensive Data Table) Excel file and converts it into a
#' `SummarizedExperiment` object. The function allows selecting specific data types and
#' optionally saves the resulting object to a file.
#'
#' @param cdt Character. Path to the CDT Excel file. Must exist.
#' @param output_file Character. Path to save the resulting `SummarizedExperiment` object.
#'   If `NULL` and `save_file` is `TRUE`, a default filename will be generated.
#' @param save_file Logical. Whether to save the resulting `SummarizedExperiment` object
#'   to a file. Default is `TRUE`. If `FALSE`, the object will be returned. If a path is not provided,
#'   a default filename will be generated based on the current date and saved in the working directory.
#' @param data_type Character. The type of data to extract from the CDT file. Options are:
#'   \itemize{
#'     \item `"peak_area"`: Extracts peak area data (sheet 4).
#'     \item `"batch_norm"`: Extracts batch normalized data (sheet 5).
#'     \item `"batch_norm_imputed"`: Extracts batch normalized and imputed data (sheet 6).
#'     \item `"mass_extracted"`: Extracts mass extracted data (sheet 7).
#'     \item `"log_transformed"`: Extracts log-transformed data (sheet 8).
#'   }
#'   Default is `"mbtch_norm_imputed"`.
#'
#' @return A `SummarizedExperiment` object containing:
#'   \itemize{
#'     \item `assays`: The assay data extracted from the CDT file.
#'     \item `rowData`: Metadata for the rows (e.g., compounds).
#'     \item `colData`: Metadata for the columns (e.g., samples).
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates the existence of the provided `cdt`.
#'   \item Reads metadata (sheet 3) and row data (sheet 2) from the CDT file.
#'   \item Extracts assay data based on the specified `data_type`.
#'   \item Ensures consistency between metadata and assay data dimensions.
#'   \item Constructs a `SummarizedExperiment` object.
#'   \item Optionally saves the resulting object to a file if `save_file` is `TRUE`.
#' }
#'
#' @examples
#' \dontrun{
#' # Convert a CDT file to a SummarizedExperiment object
#' se <- cdt_to_se(
#'   cdt = "path/to/cdt_file.xlsx",
#'   save_file = FALSE,
#'   data_type = "log_transformed"
#' )
#' }
#'
#' @importFrom openxlsx2 wb_to_df
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
cdt_to_se <- function(cdt,
                      output_file = NULL,
                      save_file = T,
                      data_type = "batch_norm_imputed",
                      rowdata_key = "INCHIKEY") {

  # Define the path to the Excel file
  if (!file.exists(cdt)) stop("The provided cdt does not exist.")
  if (is.null(output_file) && save_file) output_file <- paste0("se_from_metabolon_", Sys.Date(), ".csv")

  # Load the metadata sheet
  metadata <- make_metadata(openxlsx2::wb_to_df(cdt, sheet = 3))

  # Load the row data sheet
  rowdata <- make_rowdata(openxlsx2::wb_to_df(cdt, sheet = 2))

  # Load the assay data
  assay_data <- load_assay_data(cdt, data_type)

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
  if (save_file) {
    message("Saving results to: ", output_file)
    dir.create(dirname(output_file), showWarnings = FALSE)
    saveRDS(se, output_file)
  }
  return(se)
}

make_rowdata <- function(rowdata, rowdata_key) {

  # Check if the specified rowdata_key exists in the rowdata columns
  if (!rowdata_key %in% colnames(rowdata)) {
    stop("The provided rowdata_key is not a column in the rowdata.")
  }

  # Remove rows where all values are NA
  rowdata <- rowdata[!rowSums(is.na(rowdata)) == ncol(rowdata), ] 
  
  # Check for duplicated values in the rowdata_key column
  if (anyDuplicated(rowdata[[rowdata_key]])) {
    stop("Duplicated values found in the rowdata_key column.")
  }

  # Check for missing values in the rowdata_key column
  if (any(is.na(rowdata[[rowdata_key]]))) {
    stop("Missing values found in the rowdata_key column.")
  }

  # Set the row names of the rowdata to the values in the INCHIKEY column
  rownames(rowdata) <- rowdata[[rowdata_key]]
  
  return(rowdata)
}

make_metadata <- function(metadata) {

  metadata <- metadata[!rowSums(is.na(metadata)) == ncol(metadata), ]
  rownames(metadata) <- metadata$PARENT_SAMPLE_NAME
  
  # control_samples <- rownames(metadata[is.na(metadata[,"CLIENT_SAMPLE_ID"]), ])
  control_samples <- is.na(metadata[, "CLIENT_SAMPLE_ID"])
  # metadata <- metadata[, !rownames(metadata) %in% "CLIENT_SAMPLE_ID"]
  metadata <- metadata[!control_samples, ]

  return(metadata)
}
#

load_assay_data <- function(cdt, data_type) {
  assay_data <- switch(data_type,
    "peak_area" = {
      assay_data <- t(openxlsx2::wb_to_df(cdt, sheet = 4, row_names = TRUE))
      assay_data <- assay_data[!rowSums(is.na(assay_data)) == ncol(assay_data), ]
    },
    "batch_norm" = {
      assay_data <- t(openxlsx2::wb_to_df(cdt, sheet = 5, row_names = TRUE))
      assay_data <- assay_data[!rowSums(is.na(assay_data)) == ncol(assay_data), ]
    },
    "batch_norm_imputed" = {
      assay_data <- t(openxlsx2::wb_to_df(cdt, sheet = 6, row_names = TRUE))
      assay_data <- assay_data[!rowSums(is.na(assay_data)) == ncol(assay_data), ]
    },
    "mass_extracted" = {
      assay_data <- t(openxlsx2::wb_to_df(cdt, sheet = 7, row_names = TRUE))
      assay_data <- assay_data[!rowSums(is.na(assay_data)) == ncol(assay_data), ]
    },
    "log_transformed" = {
      assay_data <- t(openxlsx2::wb_to_df(cdt, sheet = 8, row_names = TRUE))
      assay_data <- assay_data[!rowSums(is.na(assay_data)) == ncol(assay_data), ]
    },
    stop("Invalid data type. Choose from 'peak_area', 'batch_norm', 'batch_norm_imputed', 'mass_extracted', or 'log_transformed'")
  )
}