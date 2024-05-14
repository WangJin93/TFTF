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

    if (length(dev.list()) > 0) {
    dev.off()
    }
  dd <- which(names(inter_data)=="intersection")
  if (length(dd)!=0){
    inter_data <- inter_data[-which(names(inter_data)=="intersection")]
  }
  if (length(inter_data) < 6){
    VD <- venn.diagram(inter_data, filename = NULL, fill = RColorBrewer::brewer.pal(length(inter_data), color_panel), cex = font_size,
                       margin = 0.2, cat.cex = font_size, lwd = lwd,
                       lty = rep(as.numeric(linetype), length(inter_data)))
    grid.draw(VD)

  }else{
    flowerplot(inter_data,
               ellipse_col_pal = color_panel,
               circle_col = "white",
               label_text_cex = font_size)


  }


}
