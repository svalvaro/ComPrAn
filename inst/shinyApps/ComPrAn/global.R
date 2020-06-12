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
library(shinyjs)

##global options
options(shiny.maxRequestSize = 150*1024^2)

plotBkg <- "#ecf0f5"
rbCol <- RColorBrewer::brewer.pal(3, "Set1")[1:2]

mtLSUProts <- c("Q9BYD6",    "Q5T653",    "P09001",    "Q9BYD3",    "Q9BYD2",
    "Q7Z7H8",    "Q9Y3B7",    "P52815",    "Q9BYD1",    "Q6P1L8",    "Q9P015",    "Q9NX20",    "Q9NRX2",
    "Q9H0U6",    "P49406",    "Q9BYC9",    "Q7Z2W9",    "Q9NWU5",    "Q16540",    "Q96A35",    "Q9P0M9",
    "Q13084",    "Q8TCC3",    "Q7Z7F7",    "Q9BYC8",    "O75394",    "Q9BQ48",    "Q9NZE8",    "Q9P0J6",
    "Q9BZE1",    "Q96DV4",    "Q9NYK5",    "Q9NQ50",    "Q8IXM3",    "Q9Y6G3",    "Q8N983",    "Q9H9J2",
    "Q9BRJ2",    "Q9H2W6",    "Q9HD33",    "Q96GC5",    "Q13405",    "Q8N5N7",    "Q4U2R6",    "Q86TS9",
    "Q96EL3",    "Q6P161",    "Q14197",    "Q9BQC6",    "Q8TAE8",    "Q9NP92",    "Q9NVS2")

mtSSUProts <- c("Q9Y2Q9",    "Q9Y399",    "Q96EL2",    "P82675",    "P82932",    "Q9Y2R9",    "P82933",    "P82664",
    "P82912",    "O15235",    "O60783",    "P82914",    "Q9Y3D3",    "Q9Y2R5",    "Q9Y3D5",    "P82921",
    "P82650",    "Q9Y3D9",    "P82663",    "Q9BYN8",    "Q92552",    "P51398",    "Q92665",    "Q9Y291",
    "P82930",    "P82673",    "Q96BP2",    "Q9NWT8",    "Q96EY7",    "Q9Y676")

# Read in data
#dataFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#peptides <- data.table::fread(dataFile)

