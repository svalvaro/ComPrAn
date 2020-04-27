#' Optional Filtering For Raw Data
#'
#' Filters only rows with specified values in columns Rank and Confidence Level
#' , specified as cl
#'
#' @param .data dataframe
#' @param rank integer
#' @param cl charater any combination of one or more of 'Low',
#'  'Middle', or 'High'
#'
#' @importFrom magrittr  %>%
#' @return a dataframe
#' @export
#'
#' @examples  
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#' peptides <- cleanData(data.table::fread(inputFile), fCol = "Search ID")
#' ##optional filtering based on rank and confidence level
#' peptides <- toFilter(peptides, rank = 1)
toFilter <- function(.data, rank = 1, cl = c('Low','Middle','High')) {
    .data %>%
        dplyr::filter(.data$Rank <= rank,
                        `Confidence Level` %in% cl)
    return(.data)
}