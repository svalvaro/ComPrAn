#' Clean raw peptide complexomics data
#'
#' Perform initial, mandatory, cleaning of data
#' Function to process raw input data into format required for
#' subsequent analysis. .data is a data frame containing raw input data.
#' This function checks (not neccessarily in this order):
#' \itemize{
#' \item renames Sequence ID column to Fraction and converts values in this
#'  column from letters to numbers
#' \item reorders Protein Group Accessions containing multiple proteins
#' \item removes rows in which PSM Ambiguity == 'Rejected'
#' \item removes rows in which # Protein Groups == 0
#' \item removes rows in which Precursor Area is NA
#' \item removes cols that are not used at all
#'}
#' @param .data dataframe
#' @param fCol character The column containing the fractions, 
#' e.g. "Search ID" (default)
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
#' inputFile <- system.file("extData", "data.txt", package = "ComPrAn")
#' peptides <- peptideImport(inputFile)
#' peptides <- cleanData(peptides, fCol = "Search ID")
cleanData <- function(.data, fCol = "Search ID") {
    columnHeaders <- names(.data)
    # check if column specifying fraction number is present
    if (!fCol %in% columnHeaders) {stop(paste0("\"",fCol,"\"column not found"))}
    letterSequence <- c(LETTERS, 
                        unlist(lapply(LETTERS, 
                                        function(x) {paste0(x,LETTERS)})))
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
    # Convert to factor variables
    .data %>%
        mutate( `PSM Ambiguity` = as.factor(`PSM Ambiguity`),
                `Confidence Level` = as.factor(`Confidence Level`),
                `Protein Group Accessions` = as.factor(
                    `Protein Group Accessions`)
                ) -> .data
    # Reorder multiple name accessions #
    strsplit(levels(.data$`Protein Group Accessions`), ";") %>%
        lapply(function (x) length(x)) %>% unlist() -> lengths
    # These are the names:
    strsplit(levels(.data$`Protein Group Accessions`)[lengths > 1], ";") %>%
        purrr::map(~ sort(.)) %>%
        purrr::map(~ paste(., collapse = ";")) %>%
        unlist() -> new_names
    # Replace levels
    levels(.data$`Protein Group Accessions`)[lengths > 1] <- new_names
    # Revert back to character for future work
    .data$`Protein Group Accessions` <- as.character(
        .data$`Protein Group Accessions`)
    colsToLookFor <- c( "PSM Ambiguity", "Fraction", "Precursor Area",
                        "# Protein Groups", "Rank", "Confidence Level",
                        "Protein Group Accessions", "Protein Descriptions",
                        "Modifications", "Charge", "Sequence")
    # Mandatory filtering
    .data %>%
        dplyr::select(colsToLookFor)  %>%
        dplyr::filter(`PSM Ambiguity` != 'Rejected',`# Protein Groups` > 0) %>%
        tidyr::drop_na(`Precursor Area`) %>%
        dplyr::mutate(isLabel = ifelse(str_detect(Modifications, "\\(Label\\:")
                                        ,TRUE, FALSE),
                        Sequence = toupper(Sequence)) -> .data
}