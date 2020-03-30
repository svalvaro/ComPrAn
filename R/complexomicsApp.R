#' Execute the complexomics Shiny app
#'
#' @export
complexomicsApp <- function() {
  appDir <- system.file("shinyApps", "complexomics", package = "complexomics")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `complexomics`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
