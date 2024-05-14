#' @title Get TF-target data
#' @description
#' Get the TF-target targeting results by using the api. All results saved in MySQL database.
#' @import jsonlite
#' @param table One of the dataset in "hTFtarget","KnockTF","FIMO_JASPAR","PWMEnrich_JASPAR","ENCODE","CHEA","TRRUST","GTRD","ChIP_Atlas". For correlation analysis, "cor_"+tissue type was the table name, i.e. cor_TCGA, cor_Lung.
#' @param searchType TF or Target.
#' @param gene Any TF names for searchType == TF, and any gene symbols for searchType == Target
#' @examples
#' \dontrun{
#' results <- get_data(table = "ENCODE", searchType = "Target", gene = "GAPDH")
#' }
#' @export
#'
get_data <- function(table = "ENCODE",
                     searchType = "Target",
                     gene = "GAPDH"){
  if (length(gene)>1) gene <- paste0(gene,collapse = ",")
  url <- paste0("https://www.jingege.wang/TFTF/api.php?table=",table,"&searchType=",searchType,"&searchContent=",gene)
  res <- jsonlite::fromJSON(url)
  return(res)
}
