#' Simplify Raw Proteins file
#'
#' For rows: Keep only one row with highest Precursor Area in cases where for a single Protein Group Accession in
#' a single fraction there are multiple rows with the same combination of Sequence, Mods and Charge
#'
#' For cols: remove columns that are not neccessary any more
#'
#' @param .data a dataframe
#' @param direction character, rows, cols or both
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#' peptides <- cleanData(data.table::fread(inputFile), fCol = "Search ID")
#' ## separate chemical modifications and labelling into separate columns
#' peptides <- splitModLab(peptides) 
#' ## remove unneccessary columns, simplify rows
#' peptides <- simplifyProteins(peptides) 
#'
simplifyProteins <- function(.data, direction = c("rows", "cols")) {

  colsToSelect <- c("Fraction", "Precursor Area", "Protein Group Accessions",
                    "Protein Descriptions", "Mods", "Labels",
                    "Charge", "Sequence",
                    "isLabel", "UniqueCombinedID_A","UniqueCombinedID_B")

  if(!all(colsToSelect %in% names(.data))) {
    stop('Not all columns found, you need to run splitModLab function first.') #write an ERROR message, do not continue
  }

  if ("cols" %in% direction) {
    # Remove unnecessary columns to save memory
    .data %>%
      dplyr::select(colsToSelect) -> .data

  }

  if ("rows" %in% direction) {
    .data %>%
      dplyr::group_by(`Protein Group Accessions`, UniqueCombinedID_B, Fraction) %>%
      dplyr::filter(`Precursor Area` == max(`Precursor Area`),
             !duplicated(`Precursor Area`)) -> .data
  }

}

