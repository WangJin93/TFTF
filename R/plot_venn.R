#' @title Plot venn diagram and folower plot
#' @description
#'  Plot venn diagram and folower plot to visualize the intersection of predicted results.
#' @import dplyr VennDiagram RColorBrewer grid
#' @param inter_data Intersection result obtained fron intersections() function.
#' @param color_panel Color palette in RColorBrewer, default "Set1".
#' @param font_size Default 1.
#' @param lwd Line width, default 2.
#' @param linetype "solid" = 1,"dashed" = 2, "dotted" = 3,"dot-dashed" = 4, "longdashed" = 5, "blank" = 0
#' @examples
#' \dontrun{
#' results <- predict_target(datasets=c("hTFtarget","KnockTF","FIMO_JASPAR","PWMEnrich_JASPAR"), cor_DB = c("TCGA","GTEx"),tf = "STAT3")
#' results_inter <- intersections(results)
#' plot_venn(results_inter)
#' }
#' @export
#'
plot_venn <- function(inter_data,
                      color_panel = "Set1",
                      font_size = 1,
                      lwd =2,
                      linetype = 2){
  dd <- which(names(inter_data)=="intersection")
  if (length(dd)!=0){
    inter_data <- inter_data[-which(names(inter_data)=="intersection")]
  }
  # keep only non-empty sets and avoid duplicated names breaking VennDiagram
  inter_data <- inter_data[vapply(inter_data, function(x) length(x) > 0, logical(1))]
  n <- length(inter_data)
  if (n == 0) {
    stop("No predicted results are available for plotting.")
  }
  if (n < 6){
    # VennDiagram can plot up to 5 sets; RColorBrewer palettes require >= 3
    # colors, so pad the palette for 1-2 sets as done in the TFTF app.
    if (n < 3) {
      fill <- RColorBrewer::brewer.pal(3, color_panel)[seq_len(n)]
    } else {
      fill <- RColorBrewer::brewer.pal(n, color_panel)
    }
    VD <- VennDiagram::venn.diagram(inter_data, filename = NULL, fill = fill,
                                    cex = font_size, margin = 0.2,
                                    cat.cex = font_size, lwd = lwd,
                                    lty = rep(as.numeric(linetype), n))
    grid::grid.newpage()
    grid::grid.draw(VD)
  } else {
    flowerplot(inter_data,
               ellipse_col_pal = color_panel,
               circle_col = "white",
               label_text_cex = font_size)
  }
  invisible(inter_data)
}
