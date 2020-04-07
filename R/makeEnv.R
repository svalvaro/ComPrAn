#' Create a Populated Environment
#'
#' split a large dataframe into an environment of dataframes
#'
#' @param .data peptides dataframe
#'
#' @return An environment
#' @export
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
