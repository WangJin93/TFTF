#' @title Visualization of correlation results
#' @description Presenting correlation analysis results using heat maps based on ggplot2.
#' @import dplyr psych ggplot2 reshape2 stringr
#' @param r The correlation coefficient matrix r of the correlation analysis results obtained from the functions cor_gcas_genelist(), cor_gcas_TIL(), and cor_gcas_drug().
#' @param p The P-value matrix p of the correlation analysis results obtained from the functions cor_gcas_genelist(), cor_gcas_TIL(), and cor_gcas_drug().
#' @examples
#' \dontrun{
#' genelist <- c("SIRPA","CTLA4","TIGIT","LAG3","VSIR","LILRB2","SIGLEC7","HAVCR2","LILRB4","PDCD1","BTLA")
#' dataset <-  c("GSE27262","GSE7670","GSE19188","GSE19804","GSE30219","GSE31210","GSE32665","GSE32863","GSE43458","GSE46539","GSE75037","GSE10072","GSE74706","GSE18842","GSE62113")
#' df <- get_expr_data(genes = "TNS1", datasets = dataset)
#' geneset_data <- get_expr_data(genes = genelist ,datasets = dataset)
#' result <- cor_gcas_genelist(df, geneset_data, sample_type = c("Tumor"))
#' viz_cor_heatmap(result$r, result$p)
#' }
#' @export
viz_cor_heatmap <- function(r, p) {
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    stop("The 'ggtree' package is required for this function. Please install it with: BiocManager::install('ggtree')")
  }
  if (!requireNamespace("aplot", quietly = TRUE)) {
    stop("The 'aplot' package is required for this function. Please install it with: install.packages('aplot')")
  }
  # Replace NA values in the correlation matrix with 0
  r[is.na(r)] <- 0

  # Perform hierarchical clustering if there are multiple rows
  if (nrow(r) > 1) {
    gg <- hclust(dist(r))    # Hierarchical clustering on rows
    r <- r[gg$order, ]       # Reorder rows based on clustering results
  }

  # Melt the correlation and p-value matrices
  melt_corr <- reshape2::melt(as.matrix(r))
  melt_p <- reshape2::melt(as.matrix(p))

  # Merge melted data frames on common columns
  melt_data <- merge(melt_corr, melt_p, by = c("Var1", "Var2"), sort = FALSE)
  colnames(melt_data) <- c("type", "cell_type", "corr", "PValue")

  # Assign significance stars based on PValue
  melt_data$text <- dplyr::case_when(
    melt_data$PValue < 0.001 ~ "***",
    melt_data$PValue < 0.01 ~ "**",
    melt_data$PValue < 0.05 ~ "*",
    melt_data$PValue > 0.05 ~ ""
  )

  # Reorder factor levels for cell_type
  melt_data$cell_type <- factor(melt_data$cell_type, levels = unique(melt_data$cell_type))
  melt_data$cell_type <- stringr::str_remove(melt_data$cell_type, "_XCELL")

  # Plot heatmap
  heat <- ggplot(melt_data, aes(cell_type, type)) +
    geom_tile(aes(fill = corr), colour = "grey", linewidth = 1) +
    scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.0, size = 14),
      axis.text.y = element_text(size = 14),
      plot.margin = margin(l = 30, unit = "pt")
    ) +
    geom_text(aes(label = text), col = "black", size = 7) +
    labs(fill = paste0("\n", "\n", "\n", "\n", "\n", "  *   p < 0.05", "\n", " **  p < 0.01", "\n",
                       "*** p < 0.001", "\n\n", "Correlation")) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "right")

  # Add hierarchical clustering dendrogram if there are multiple rows
  if (nrow(r) > 1) {
    h <- ggtree::ggtree(gg, layout = "rectangular", branch.length = "none")
    heat <- heat %>% aplot::insert_left(h, width = 0.1)
  }

  return(heat)
}
