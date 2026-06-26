#' @title Plot Venn Diagram for Differentially Expressed Genes
#' @description
#' This function plots a Venn diagram for lists of differentially expressed genes (DEGs) across multiple datasets.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid grid.newpage grid.draw
#' @param results List of character vectors. Each vector contains DEGs for a specific dataset.
#' @param fill_colors Character vector. Colors to fill the Venn diagram circles. Default is NULL, which uses a palette.
#' @param palette Character. Name of the RColorBrewer palette to use if fill_colors is not specified. Default is "Set1".
#' @param lty Numeric. Line type for the circles in the Venn diagram. Default is 2 (dashed line).
#' @param ... Additional arguments passed to venn.diagram function.
#' @return A list of intersected DEGs.
#' @examples
#' \dontrun{
#' df1 <- get_OSF_data(table = "GSE31210", action = "geo_data")
#' results1 <- DEGs_analysis(df1)
#' df2 <- get_OSF_data(table = "GSE19188", action = "geo_data")
#' results2 <- DEGs_analysis(df2)
#' DEGs_lists <- list("GSE31210" = results1, "GSE19188" = results2)
#' results <- get_DEGs_list(DEGs_lists)
#' plot_venn(results$DEG_up, palette = "Set1")
#' plot_venn(results$DEG_up, fill_colors = c("red", "green", "blue"), alpha = 0.5, cex = 1.5)
#' }
#' @export
plot_venn <- function(results, fill_colors = NULL, palette = "Set1", lty = 2, ...) {
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    stop("The 'VennDiagram' package is required for this function. Please install it with: install.packages('VennDiagram')")
  }

  # 确定数据集的数量
  intersects <- Reduce(intersect, results)
  grid::grid.newpage()
  num_sets <- length(results)

  # 使用RColorBrewer选择配色方案（如果未提供fill_colors参数）
  if (is.null(fill_colors)) {
    fill_colors <- RColorBrewer::brewer.pal(min(num_sets, RColorBrewer::brewer.pal.info[palette, "maxcolors"]), palette)
  }

  # 绘制Venn图，并传入其他参数
  VD <- VennDiagram::venn.diagram(results, filename = NULL, fill = fill_colors[1:num_sets],
                                  lty = rep(lty, num_sets), ...)

  # 绘制Venn图
  grid::grid.draw(VD)

  return(intersects)
}
