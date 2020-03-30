#' Execute the complexomics Shiny app
#'
#' @export
compranApp <- function() {
  appDir <- system.file("shinyApps", "ComPrAn", package = "ComPrAn")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `ComPrAn`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
