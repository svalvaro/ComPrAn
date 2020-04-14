#' Create a Populated Environment
#'
#' split a large dataframe into an environment of dataframes
#'
#' @param .data peptides dataframe
#'
#' @return An environment
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
#' ##make environment
#' peptide_index <- makeEnv(peptides)
#' 
makeEnv <- function(.data) {
  # peptide_index <<- new.env(hash = TRUE, parent = emptyenv())
  # for (i in unique(.data$`Protein Group Accessions`)) {
  #   peptide_index[[i]] <<- .data[.data$`Protein Group Accessions` == i,]
  # }

  # Don't use <<-
  peptide_index <- new.env(hash = TRUE, parent = emptyenv())
  for (i in unique(.data$`Protein Group Accessions`)) {
    peptide_index[[i]] <- .data[.data$`Protein Group Accessions` == i,]
  }

  return(peptide_index)
}
