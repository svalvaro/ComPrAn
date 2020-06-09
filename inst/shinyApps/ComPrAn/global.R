##load libraries
library(ComPrAn)
library(shinydashboard)
library(scales)
library(VennDiagram)
library(RColorBrewer)
library(DT)
library(rio)
library(grid)
library(ggplot2)
library(readr)

##global options
options(shiny.maxRequestSize = 150*1024^2)

plotBkg <- "#ecf0f5"
rbCol <- RColorBrewer::brewer.pal(3, "Set1")[1:2]

# Read in data
#dataFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#peptides <- data.table::fread(dataFile)

