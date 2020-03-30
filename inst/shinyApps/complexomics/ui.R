###################
# ui.R
# 
###################
source('./app_components/header.R')
source('./app_components/sidebar.R')
source('./app_components/body.R')


ui <- dashboardPage(header, sidebar, body, skin = "purple")
