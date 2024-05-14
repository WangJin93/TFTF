#' @title Get intersections of prediction result
#' @description
#'  Get the intersections of the prediction result of predict_target() and predict_TF().
#'
#' @import dplyr
#' @param results The result of predict_target() and predict_TF().
#' @param datasets_index datasets_index of user inputted datasets. Default NULL, return the intersection of all datasets.
#' @examples
#' \dontrun{
#' results <- predict_target(datasets=c("hTFtarget","KnockTF","FIMO_JASPAR","PWMEnrich_JASPAR"), cor_DB = c("TCGA","GTEx"),tf = "STAT3")
#' results_inter <- intersections(results)
#' }
#' @export
#'
intersections <- function(results,datasets_index=NULL){
  if (is.null(datasets_index)) datasets_index = 1:(length(results)-1)
  results <- results$results[datasets_index]
  intersection <- Reduce(intersect,results)
  if(is.null(intersection)){
    intersection ="None"
  }
  results$intersection <- intersection
  return(results )
}
