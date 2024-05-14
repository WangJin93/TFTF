#' @title  Tissue type information
#' @description
#' Get the tissue types in "TCGA" or "GTEx" database
#' @import dplyr
#' @param db "TCGA" or "GTEx".
#' @examples
#' \dontrun{
#' tissue_type("GTEx")
#' }
#' @export
#'
tissue_type <- function(db = "TCGA"){
  if (db == "TCGA"){
    print(tissue[,"TCGA"] %>% as.character())
  }else  if (db == "GTEx"){
    print(tissue[,"GTEx"] %>% na.omit() %>% as.character())
  }else{
    cat("Only TCGA and GTEx supported!")
  }

}
