#' @title Get dataset information
#' @description
#' Get the data information of specific Transcription Factor dataset.
#' @import dplyr
#' @param dataset Transcription Factor datasets used in this package, e.g. "ENCODE".
#' Both the internal database names (e.g. "JASPAR") and the table names used by
#' \code{get_data()} / \code{predict_target()} (e.g. "FIMO_JASPAR",
#' "PWMEnrich_JASPAR", "cor_COAD") are accepted.
#' @examples
#' \dontrun{
#' get_data_info("ENCODE")
#' get_data_info(c("FIMO_JASPAR","hTFtarget"))
#' }
#' @export
#'
get_data_info <- function(dataset = "ENCODE"){
  if (!is.character(dataset) || !length(dataset)){
    stop("dataset must be a non-empty character vector.")
  }
  # map prediction table names to the corresponding data_info entries
  alias <- c("FIMO_JASPAR" = "JASPAR",
             "PWMEnrich_JASPAR" = "JASPAR")
  ds <- unname(ifelse(dataset %in% names(alias), alias[dataset], dataset))
  if (all(ds %in% names(data_info))){
    for (d in unique(ds)){
      cat("----", d, "----\n")
      print(data_info[[d]])
    }
    invisible(data_info[unique(ds)])
  }else{
    cat("input error! try data(data_info) to see the available dataset names.\n")
    invisible(NULL)
  }
}
