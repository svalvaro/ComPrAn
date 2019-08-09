#' Optional Filtering For Raw Data
#'
#' Filters only rows with specified values in columns Rank and Confidence Level, specified as cl
#'
#' @param .data dataframe
#' @param rank integer
#' @param cl charater any combination of one or more of 'Low', 'Middle', or 'High'
#'
#' @importFrom magrittr  %>%
#' @return a dataframe
#' @export
#'
toFilter <- function(.data, rank = 1, cl = c('Low','Middle','High')) {
  #This function filters only rows with specified values in columns Rank and`Confidence Level`
  # rank numeric
  # cl character vector of some or all of 'Low', 'Middle' and 'High' values
  .data %>%
    dplyr::filter(Rank <= rank,
           `Confidence Level` %in% cl)
}


