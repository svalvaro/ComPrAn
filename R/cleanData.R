#' Clean raw peptide complexomics data
#'
#' Perform initial, mandatory, cleaning of data
#' Check presence of reguired columns
#' Function to process raw input data into format required for subsequent analysis
#' .data is a data frame containing raw input data
#' This function checks (not neccessarily in this order):
#' \itemize{
#' \item are all required columns present
#' \item are these columns in correct format
#' \item renames Sequence ID column to Fraction and converts values in this column from letters to numbers
#' \item reorders Protein Group Accessions containing multiple proteins
#' \item removes rows in which PSM Ambiguity == 'Rejected'
#' \item removes rows in which # Protein Groups == 0
#' \item removes rows in which Precursor Area is NA
#' \item removes cols that are not used at all
#'}
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
#' @author Petra Palenikova \email{pp451@@cam.ac.uk}
#' @author Rick Scavetta \email{office@@scavetta.academy}
#' 
#' @examples  
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#' peptides <- cleanData(data.table::fread(inputFile), fCol = "Search ID")
cleanData <- function(.data, fCol = "Search ID") {

  columnHeaders <- names(.data)

  #rename "Search ID" to something more meaningful e.g. "Fraction"
  #check if column specifying fraction number is present
  # if (fCol %in% columnHeaders) {
  #   names(.data)[grep(fCol, names(.data))] <- "Fraction"
  # } else {
  #   stop(paste(fCol, "column not found"))
  # }
  
  # check if column specifying fraction number is present
  if (!fCol %in% columnHeaders) {
      stop(paste0("\"",fCol, "\"column not found"))
      } 
  
  # convert Fraction (originaly "Search ID") to numbers instead of letters:
  ## This implementation causes issues when there are missing
  ## letters (= fractions).
  ## Let's say fractions present are A,B,C,F,G they are re-numbered 
  ## to 1,2,3,4,5 instead of 1,2,3,6,7
  # .data %>%
  #   dplyr::mutate(n = nchar(Fraction)) %>%
  #   dplyr::arrange(n, Fraction) %>%
  #   dplyr::mutate(Fraction = as.factor(Fraction),
  #          Fraction = as.integer(Fraction)) %>%
  #   dplyr::select(-n) -> .data
  
  ## Alternative implementation to account for issue described above
  # convert Fraction (originaly "Search ID") to numbers instead of letters:
  # 1) generate a data frame that we will use to map correce numers to letters
  #    let's make this work for up to 702 fractions (27x length of LETTERS)
  # 2) merge this data frame with .data to convert letters to numbers
  #    -> add "Fraction" column with correct numers, remove old fcol column
  letterSequence <- c(LETTERS, 
                      unlist(lapply(LETTERS, function(x) {paste0(x,LETTERS)})))
  
  letterToNumber <- data.frame("Letters" = letterSequence,
                               "Fraction" = seq_along(letterSequence),
                               stringsAsFactors = FALSE)
  
  names(letterToNumber)[grep("Letters", names(letterToNumber))] <- fCol
  nrowOld <- nrow(.data)
  .data <- dplyr::inner_join(.data,letterToNumber, by = fCol)
  
  # check if column specifying fraction number is present
  if (nrow(.data) < nrowOld) {
      stop(paste0("Some rows in \"", fCol, "\" column do not have a value"))
  } 
  

  #all column names with after adding "Fraction" column
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
