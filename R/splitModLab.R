#' Split Modification and Label tags
#'
#' Splits up the Modifications column into lists of vectors for 
#' modifications(Mods) and labels(Labels)
#' It adds two more columns to the data frame:
#' \itemize{
#' \item UniqueCombinedID_A: Unique combinations of Sequence, Mods and Charge 
#' for "scenario A".
#' \item UniqueCombinedID_B: Unique combinations of Sequence, Mods, Charge and 
#' Labels for "scenario B"
#'}
#' @param .data dataframe
#'
#' @importFrom stringr str_split
#'
#' @return dataframe
#' @export
#' @examples 
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#' peptides <- peptideImport(inputFile)
#' peptides <- cleanData(peptides, fCol = "Search ID")
#' ## separate chemical modifications and labelling into separate columns
#' peptides <- splitModLab(peptides) 
#'
splitModLab <- function(.data) {
    # split up all elements:
    .data$splitMod <- str_split(.data$Modifications, "; ") 
    # Give a cataloge of all the unique values in the data set:
    presentAll <- unique(unlist(.data$splitMod))
    presentLabels <- stringr::str_extract(presentAll, "^.*Label.*$")
    presentMods <- presentAll[is.na(presentLabels)]
    presentLabels <- presentLabels[!is.na(presentLabels)]
    # Make a column for mods and one for labels:
    .data$Mods <- .data$splitMod
    .data$Labels <- .data$splitMod
    #keep only ones belonging to Mods/Labels and sort them into same order
    .data$Mods <- .data$Mods %>%
        purrr::map(~ sort(.[. %in% presentMods]))
    .data$Labels <- .data$Labels %>%
        purrr::map(~ sort(.[. %in% presentLabels]))
    # To remove character(0)
    .data$Mods <- lapply(.data$Mods, 
                        function(x) if(identical(x, character(0))) "" else x)
    .data$Labels <- lapply(.data$Labels, 
                        function(x) if(identical(x, character(0))) "" else x)
    .data$Mods <- as.character(.data$Mods)
    .data$Labels <- as.character(.data$Labels)
    # Create unique IDs for scenario A or B:
    .data$UniqueCombinedID_A <- paste(  .data$Sequence,
                                        .data$Mods,
                                        .data$Charge, sep = "_")
    
    .data$UniqueCombinedID_B <- paste(  .data$Sequence,
                                        .data$Mods,
                                        .data$Charge,
                                        .data$Labels, sep = "_")
    .data %>%
        dplyr::select(-c(Modifications, .data$splitMod)) -> .data
    return(.data)
}