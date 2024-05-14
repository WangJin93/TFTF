#' Run TFTF Shiny App
#' @import shiny
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' TFTF_app()
#' }
TFTF_app <- function() {
  shiny::shinyAppFile(system.file("shinyapp", "App.R", package = "TFTF"))
}
