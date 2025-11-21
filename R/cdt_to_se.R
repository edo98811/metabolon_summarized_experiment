#' Convert CDT File to SummarizedExperiment Object
#'
#' This function reads a CDT (Client Data Table) Excel file and converts it into a
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
#'   Default is `"batch_norm_imputed"`.
#' @param rowdata_key Character. The column name in the row data (sheet 2) to use as row names.
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
#'   \item Validates the existence of the provided CDT.
#'   \item Reads metadata (sheet 3) and row data (sheet 2) from the CDT file.
#'   \item Extracts assay data based on the specified `data_type`.
#'   \item Cleans and formats the metadata and row data. The componds kept all the ones that have a non-NA value in the selected row data (annotation) column.
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
#' @importFrom readxl read_excel
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
  metadata <- make_metadata(readxl::read_excel(cdt, sheet = 3))

  # Load the assay data
  assay_data <- load_assay_data(cdt, data_type)

  # Load the row data sheet
  rowdata <- make_rowdata(readxl::read_excel(cdt, sheet = 2), rowdata_key = rowdata_key)

  # Ensure that the rownames of the assay_data match the rownames of the rowdata (delete absent rows)
  assay_data <- assay_data[rownames(assay_data) %in% rowdata$CHEM_ID, ]
  
  # Sort assay_data to match the order of rowdata
  assay_data <- assay_data[match(rowdata$CHEM_ID, rownames(assay_data)), ]
  rownames(assay_data) <- rownames(rowdata)

  # Check that the metadata dataframes are correct
  stopifnot(all(rownames(rowdata) == rownames(assay_data)))
  stopifnot(all(rownames(metadata) == colnames(assay_data)))

  # Create the SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
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

#' Clean and format the row data from the CDT file
#'
#' @param rowdata A data frame containing the row data from the CDT file.
#' @param rowdata_key A character string specifying the column name to use as row names.
#' @return A cleaned and formatted row data data frame.
make_rowdata <- function(rowdata, rowdata_key) {
  # Check if the specified rowdata_key exists in the rowdata columns
  if (!rowdata_key %in% colnames(rowdata)) {
    stop("The provided rowdata_key is not a column in the rowdata.")
  }

  # Remove rows where all values are NA
  rowdata <- rowdata[!rowSums(is.na(rowdata)) == ncol(rowdata), ]

  # Check for missing values in the rowdata_key column and delete those rows
  na_keys <- is.na(rowdata[[rowdata_key]])
  if (any(na_keys)) {
    warning(
      "Missing values found in the rowdata_key column. The following entities will be removed: ",
      paste(rowdata[is.na(rowdata[[rowdata_key]]), ]$CHEM_ID, collapse = ", ")
    )
    rowdata <- rowdata[!na_keys, ]
  }

  # Check for duplicated values in the rowdata_key column. Kept only the first occurrence. This is to ensure that the rownames are unique
  duplicated_rows <- rowdata[duplicated(rowdata[[rowdata_key]]), ]
  if (nrow(duplicated_rows) > 0) {
    warning(sprintf(
      "Duplicated values found in the %s (rowdata_key) column. The following duplicates will be removed: %s",
      rowdata_key,
      paste(duplicated_rows[[rowdata_key]], collapse = ", ")
    ))
    rowdata <- rowdata[!duplicated_rows, ]
  }

  # Set the row names of the rowdata to the values in the INCHIKEY column
  rowdata <- data.frame(rowdata, row.names = rowdata[[rowdata_key]])

  return(rowdata)
}

#' Clean and format the metadata from the CDT file
#'
#' @param metadata A data frame containing the metadata from the CDT file.
#' @return A cleaned and formatted metadata data frame.
make_metadata <- function(metadata) {
  metadata <- metadata[!rowSums(is.na(metadata)) == ncol(metadata), ]

  control_samples <- is.na(metadata[, "CLIENT_SAMPLE_ID"])
  metadata <- metadata[!control_samples, ]

  metadata <- data.frame(metadata, row.names = metadata$PARENT_SAMPLE_NAME)
  return(metadata)
}

#' Load assay data based on the specified data type
#'
#' @param cdt Character. Path to the CDT Excel file.
#' @param data_type Character. The type of data to extract. Options are:
#'   "peak_area", "batch_norm", "batch_norm_imputed", "mass_extracted", "log_transformed".
#' @importFrom readxl read_excel
#' @return A data frame containing the assay data.
load_assay_data <- function(cdt, data_type) {
  # Determine which sheet to read
  sheet <- switch(data_type,
    "peak_area"           = 4,
    "batch_norm"          = 5,
    "batch_norm_imputed"  = 6,
    "mass_extracted"      = 7,
    "log_transformed"     = 8,
    stop("Invalid data type. Choose from 'peak_area', 'batch_norm', 'batch_norm_imputed', 'mass_extracted', or 'log_transformed'")
  )

  # Read and transpose the sheet
  raw_data <- t(readxl::read_excel(cdt, sheet = sheet, .name_repair = "minimal"))

  # Remove rows with all NAs
  raw_data <- raw_data[!apply(raw_data, 1, function(x) all(is.na(x))), ]

  # First row contains sample names
  sample_names <- raw_data[1, ]

  # Remaining rows are numeric data
  numeric_data <- raw_data[-1, , drop = FALSE]

  # Convert all values to numeric
  numeric_data <- apply(numeric_data, 2, as.numeric)

  # Assign column names and row names
  colnames(numeric_data) <- sample_names
  rownames(numeric_data) <- rownames(raw_data)[-1]

  return(numeric_data)
}
