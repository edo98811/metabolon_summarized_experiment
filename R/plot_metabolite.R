#' Plot Metabolite Abundance
#'
#' This function generates a ggplot visualization of the abundance of a specified metabolite across samples in a SummarizedExperiment object.
#'
#' @param metabolite Character. The name of the metabolite to plot. Must be present in the rowData of \code{se_object}.
#' @param se_object SummarizedExperiment. The object containing metabolite data and sample metadata.
#' @param group_col Character (optional). The column name in colData to use for grouping samples in the plot. If NULL, all samples will be grouped together.
#' @param title Character (optional). The title for the plot.
#' @param rm.na Logical (default TRUE). Whether to remove samples with NA values in the metabolite abundance or grouping column.
#'
#' @return A ggplot object visualizing the abundance of the specified metabolite.
#'
#' @examples
#' # plot_metabolite("00112", se_object, group_col = "Condition", title = "Abundance")
#'
#' @import ggplot2
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom ggforce geom_sina
#' @importFrom ggrepel geom_text_repel
#'
#' @export
plot_metabolite <- function(metabolite, se_object, group_col = NULL, title = NULL, rm.na = TRUE) {
  if (!inherits(se_object, "SummarizedExperiment")) {
    stop("se_object must be a SummarizedExperiment.")
  }
  if (!(metabolite %in% rownames(SummarizedExperiment::rowData(se_object)))) {
    stop("Metabolite not found in rowData of se_object.")
  }
  if (!is.null(group_col) && !(group_col %in% colnames(SummarizedExperiment::colData(se_object)))) {
    stop("group_col not found in colData of se_object.")
  }

  # Extract the data for the specified metabolite
  metabolite_data <- data.frame(abundance = SummarizedExperiment::assay(se_object)[metabolite, ])
  
  # Handle group_col properly - create appropriate metadata dataframe
  if (!is.null(group_col)) {
    sample_metadata <- data.frame(group = SummarizedExperiment::colData(se_object)[[group_col]])
    names(sample_metadata) <- group_col
  } else {
    sample_metadata <- data.frame(group = rep("All Samples", ncol(se_object)))
    names(sample_metadata) <- "group"
    group_col <- "group"
  }

  # Remove NA values if requested
  if (rm.na) {
    if (!is.null(group_col) && group_col != "group") {
      valid_idx <- !is.na(metabolite_data[["abundance"]]) & !is.na(sample_metadata[[group_col]])
    } else {
      valid_idx <- !is.na(metabolite_data[["abundance"]])
    }
    metabolite_data <- metabolite_data[valid_idx, , drop = FALSE]
    sample_metadata <- sample_metadata[valid_idx, , drop = FALSE]
  }
  
  # Check if we have any data left
  if (nrow(metabolite_data) == 0) {
    stop("No valid data points remaining after NA removal.")
  }
  
  plot_data <- cbind(metabolite_data, sample_metadata)
  plot_data$label <- rownames(plot_data)
  
  p <- ggplot(plot_data, aes(x = .data[[group_col]], y = abundance, label = label, color = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, fill = "#A3C4DC", color = scales::alpha("#375A7F", 0.3), alpha = 0.3) +
    ggforce::geom_sina(size = 2, alpha = 0.8) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 10) +
    labs(
      title = ifelse(!is.null(title), title, paste0("Abundance of ", metabolite)),
      x = ifelse(group_col == "group", "Samples", group_col),
      y = "Abundance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    scale_color_brewer(palette = "Dark2")

  return(p)
  
}
