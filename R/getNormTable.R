#' Get normalised table for all proteins 
#'
#' Extracts values for representative peptides for each protein, 
#' for both scenario A and scenario B.
#' Results are combined into one data frame in a format either indended for
#' further analysis or for export.
#' @param .listDf list of data frames
#' @param purpose character, purpose of use of function output, values either
#'  "analysis" of "export"
#'
#' @importFrom purrr map_df
#' 
#' @return dataframe
#' @export
#' 
#' @examples 
#' 
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#' peptides <- cleanData(data.table::fread(inputFile), fCol = "Search ID")
#' ## separate chemical modifications and labelling into separate columns
#' peptides <- splitModLab(peptides) 
#' ## remove unneccessary columns, simplify rows
#' peptides <- simplifyProteins(peptides) 
#' ## Pick representative peptide for each protein for both scenarios
#' peptide_index <- pickPeptide(peptides)
#' ## extract table with normalised protein values for both scenarios
#' forAnalysis <- getNormTable(peptide_index,purpose = "analysis")
#' 
getNormTable <- function(.listDf, purpose = "analysis"){
    if (purpose != "analysis" & purpose != "export"){
        stop('Valid values for purpose are "analysis" or "export".')
    }
    names(.listDf) %>%
        map_df(~ extractRepPeps(
            .listDf[[.]], scenario = 'A', label = TRUE))  %>%
        normalizeTable() -> protNormLab
    
    names(.listDf) %>%
        map_df(~ extractRepPeps(
            .listDf[[.]], scenario = 'A', label = FALSE))  %>%
        normalizeTable() -> protNormUnlab
    
    names(.listDf) %>%
        map_df(~ extractRepPeps(
            .listDf[[.]], scenario = 'B')) %>%
        normalizeTable() -> protNormComb
    if(purpose == "analysis"){
        output <- normTableWideToLong(
            protNormLab, protNormUnlab, protNormComb)
    }
    else if(purpose == "export"){
        output <- normTableForExport(
            protNormLab, protNormUnlab, protNormComb)
    }
    return(output)
}