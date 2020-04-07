###################
# ui.R
#
###################
# source('./app_components/header.R')
# source('./app_components/sidebar.R')
# source('./app_components/body.R')
#
# ui <- dashboardPage(header, sidebar, body, skin = "purple")

source(system.file("shinyApps/complexomics/app_components", "header.R", package = "ComPrAn"))
source(system.file("shinyApps/complexomics/app_components", "sidebar.R", package = "ComPrAn"))
source(system.file("shinyApps/complexomics/app_components", "body.R", package = "ComPrAn"))

ui <- dashboardPage(header, sidebar, body, skin = "purple")

