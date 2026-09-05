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
  app_dir <- system.file("shinyapp", package = "TFTF")
  app_file <- ""
  if (nzchar(app_dir)) {
    f <- file.path(app_dir, "App.R")
    if (!file.exists(f)) f <- file.path(app_dir, "app.R")  # fallback
    if (file.exists(f)) app_file <- f
  }
  if (!nzchar(app_file)) {
    stop("Could not locate the bundled Shiny application files (inst/shinyapp/). ",
         "Please reinstall the package, e.g. devtools::install_github(\"WangJin93/TFTF\").",
         call. = FALSE)
  }
  shiny::shinyAppFile(app_file)
}
