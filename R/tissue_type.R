#' @title  Tissue type information
#' @description
#' Get the tissue types in "TCGA" or "GTEx" database
#' @import dplyr
#' @param db "TCGA" or "GTEx".
#' @return A character vector of tissue/cancer type names (invisibly); the
#' vector is also printed.
#' @examples
#' \dontrun{
#' tissue_type("GTEx")
#' }
#' @export
#'
tissue_type <- function(db = "TCGA"){
  if (db == "TCGA"){
    res <- as.character(tissue[["TCGA"]])
  }else  if (db == "GTEx"){
    res <- as.character(tissue[["GTEx"]])
    res <- res[!is.na(res)]
  }else{
    cat("Only TCGA and GTEx supported!")
    return(invisible(NULL))
  }
  print(res)
  invisible(res)
}
