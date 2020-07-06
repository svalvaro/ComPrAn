#' Execute the complexomics Shiny app
#'
#' @import scales VennDiagram RColorBrewer DT shinydashboard rio grid shinyjs
#' @importFrom data.table fread
#' @return Shiny app
#' @export
#' 
#' @examples
#' ##to start the shiny app associated with ComPrAn package run
#' compranApp()
compranApp <- function() {
    appDir <- system.file("shinyApps", "ComPrAn", package = "ComPrAn")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `ComPrAn`.", 
                call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
}
