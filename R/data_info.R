#' @title Get dataset information
#' @description
#' Get the data information of specific Transcription Factor dataset.
#' @import dplyr
#' @param dataset Transcription Factor datasets used in this package, e.g. "ENCODE".
#' @examples
#' \dontrun{
#' get_data_info("ENCODE")
#' }
#' @export
#'
get_data_info <- function(dataset = "ENCODE"){
  if (all(dataset %in% names(data_info))){
    print(data_info[dataset])
  }else{
    cat("input error! try data(data_info)")
  }
}
