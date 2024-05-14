#' Run TFTF Shiny App
#' @importFrom shiny validate
#' @import shiny
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' app_run()
#' }
TFTF_app <- function() {
  shiny::shinyAppFile(system.file("shinyapp", "App.R", package = "TFTF"))
}
