#' Execute the complexomics Shiny app
#'
#' @import scales VennDiagram RColorBrewer DT shinydashboard rio grid
#' @importFrom data.table fread
#'
#' @export
complexomicsApp <- function() {
  appDir <- system.file("shinyApps", "complexomics", package = "ComPrAn")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `ComPrAn`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
