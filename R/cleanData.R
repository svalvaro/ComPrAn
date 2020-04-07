#' Clean raw complexomic data
#'
#' Perform initial, mandatory, cleaning of data
#' Check presence of reguired columns
#' Function to process raw input data into format required for subsequent analysis
#' .data is a data frame containing raw input data
#' This function checks (not neccessarily in this order):
#' 1) are all required columns present
#' 2) are these columns in correct format
#' 3) renames Sequence ID column to Fraction and converts values in this column from letters to numbers
#' 4) reorders Protein Group Accessions containing multiple proteins
#' 5) removes rows in which PSM Ambiguity == 'Rejected'
#' 6) removes rows in which # Protein Groups == 0
#' 7) removes rows in which Precursor Area is NA
#' 8) removes cols that are not used at all
#'
#' @param .data dataframe
#' @param fCol character The column containing the fractions, e.g. "Search ID" (default)
#'
#' @import dplyr forcats
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_detect
#'
#' @return dataframe
#' @export
#'
cleanData <- function(.data, fCol = "Search ID") {

  columnHeaders <- names(.data)

  #rename "Search ID" to something more meaningful e.g. "Fraction"
  if (fCol %in% columnHeaders) {
    names(.data)[grep(fCol, names(.data))] <- "Fraction"
  } else {
    stop(paste(fCol, "column not found"))
  }

  #all column names with "Search ID" renamed
  columnHeaders <- names(.data)

  # which columns we need
  colsToLookFor <- c("PSM Ambiguity", "Fraction", "Precursor Area",
                     "# Protein Groups", "Rank", "Confidence Level",
                     "Protein Group Accessions", "Protein Descriptions",
                     "Modifications", "Charge", "Sequence")

  #their expected classes
  classesOfCols <- c("factor", "integer", "numeric",
                     "integer", "integer", "factor",
                     "factor", "character",
                     "character", "integer", "character")


  # make a named vector:
  names(classesOfCols) <- colsToLookFor

  #are all columns we need in column names?
  if(sum(colsToLookFor %in% columnHeaders) < length(colsToLookFor)) {
    stop('Not all columns found') #write an ERROR message, do not continue
  }

  # convert Fraction (originaly "Search ID") to numbers instead of letters:
  .data %>%
    mutate(n = nchar(Fraction)) %>%
    arrange(n, Fraction) %>%
    mutate(Fraction = as_factor(Fraction),
                  Fraction = as.integer({Fraction})) %>%
    select(-n) -> .data

  # Convert to factor variables
  .data %>%
    mutate(`PSM Ambiguity` = as.factor(`PSM Ambiguity`),
                  `Confidence Level` = as.factor(`Confidence Level`),
                  `Protein Group Accessions` = as.factor(`Protein Group Accessions`)) -> .data

  #are all values in columns we need of correct type?
  # check to see that actual and desired match
  if(!identical(sapply(.data, class)[colsToLookFor], classesOfCols)) {
    stop('Not all columns are the correct type!') #write an ERROR message, do not continue
  }

  # Reorder multiple name accessions #
  strsplit(levels(.data$`Protein Group Accessions`), ";") %>%
    lapply(function (x) length(x)) %>%
    unlist() -> lengths

  # These are the names:
  strsplit(levels(.data$`Protein Group Accessions`)[lengths > 1], ";") %>%
    purrr::map(~ sort(.)) %>%
    purrr::map(~ paste(., collapse = ";")) %>%
    unlist() -> new_names

  # Replace levels
  levels(.data$`Protein Group Accessions`)[lengths > 1] <- new_names

  # Revert back to characte for future work
  .data$`Protein Group Accessions` <- as.character(.data$`Protein Group Accessions`)

  # Mandatory filtering
  .data %>%
    dplyr::select(colsToLookFor)  %>%
    dplyr::filter(`PSM Ambiguity` != 'Rejected',
                  `# Protein Groups` > 0) %>%
    tidyr::drop_na(`Precursor Area`) %>%
    dplyr::mutate(isLabel = ifelse(str_detect(Modifications, "\\(Label\\:"), TRUE, FALSE),
                  Sequence = toupper(Sequence)) -> data1

}
