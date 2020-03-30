###################
# global.R
#
###################
library(data.table)
library(tidyverse)
library(scales)
library(VennDiagram)
library(RColorBrewer)
library(DT)
library(shinydashboard)
library(rio)
library(grid)
library(gridExtra)
library(complexomics)

# options(shiny.maxRequestSize = 30*1024^2)
options(shiny.maxRequestSize = 150*1024^2)

# global vars
plotBkg <- "#ecf0f5"
rbCol <- brewer.pal(3, "Set1")[1:2]

col_vector_peptides <- c('#ffc125', '#a020f0')
names(col_vector_peptides) <- c('TRUE', 'FALSE')

col_vector_proteins <- c('#ff9d2e', '#07b58a')
names(col_vector_proteins) <- c('TRUE', 'FALSE')

# This sources all files in the R directory
# for (i in list.files("R/", full.names = T)) {
#   source(i)
# }

# Read in data
# dataFile <- system.file("extdata", "data.txt", package = "complexomics")
# peptides <- data.table::fread(dataFile)
# names(peptides)
# max_frac <- max(peptides$Search.ID)
# max_frac <- max(peptides$Fraction)

# peptides <- peptides[peptides$`Protein Group Accessions` %in% c("Q9Y399", "P61604", "P27797", "Q8WW59", "O75439"),]

# peptides <- cleanData(peptides, fCol = "Search ID")

# peptides <- cleanData(data.table::fread(dataFile), fCol = "Search ID")
# names(peptides)


# splitModLab(peptides[1:5,])
# peptides <- splitModLab(peptides)

# set.seed(122)
# df <- data.frame(
#   n = rnorm(500),
#   m = rnorm(100))
