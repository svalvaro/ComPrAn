#' Import raw peptide complexomics data
#'
#' Check presence of required columns
#' inputFile is a character vector containing the location of peptide file
#' This function checks:
#' \itemize{
#' \item are all required columns present
#' \item are these columns in correct format
#'}
#' @param inputFile character
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
#' ##Use example peptide data set, read in data
#' inputFile <- system.file("extData", "data.txt", package = "ComPrAn")
#' peptides <- peptideImport(inputFile)
peptideImport <- function(inputFile) {
    .data <- data.table::fread(inputFile, stringsAsFactors = FALSE)
    columnHeaders <- names(.data)

    # which columns we need
    colsToLookFor <- c( "PSM Ambiguity", "Search ID", "Precursor Area",
                        "# Protein Groups", "Rank", "Confidence Level",
                        "Protein Group Accessions", "Protein Descriptions",
                        "Modifications", "Charge", "Sequence")
    #their expected classes
    classesOfCols <- c( "character", "character", "numeric",
                        "integer", "integer", "character",
                        "character", "character",
                        "character", "integer", "character")
    # make a named vector:
    names(classesOfCols) <- colsToLookFor
    #are all columns we need in column names?
    if(sum(colsToLookFor %in% columnHeaders) < length(colsToLookFor)) {
        stop(
        c("\n", 
            "Not all columns found!\nRequired columns and their classes are:\n"
            ,"\n",
            paste0(names(classesOfCols),": ",classesOfCols,"\n"))) 
    } #write an error message, do not continue

    #are all values in columns we need of correct type?
    # check to see that actual and desired match
    if(!identical(vapply(.data, class,"a")[colsToLookFor], classesOfCols)) {
        warning(c('\n Not all columns are the correct type!','\n ',
        'ComPrAn functions might not work properly with  this data set.',
        '\n Expected column types are:\n',
        paste0(names(classesOfCols),": ",classesOfCols,"\n")))
    }
    return(.data)
}