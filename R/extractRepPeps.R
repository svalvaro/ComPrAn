#' Extract Only Data Belonging to Representative Peptide
#'
#' Incomplete labelling - there are cases when in peptides containing multiple Lys/Arg
#' not all of them are heavy in labelled samples. As in SILAC we assume that addition of
#' label does not affect peptide properties, we are taking a mean `Precursor Area` value
#' as the representative `Precursor Area` in such cases
#'
#' @param .data dataframe containing all peptides of one protein
#' @param scenario character A, or B
#' @param label character for output label
#'
#' @return dataframe containing only representative peptide
#' @export
extractRepPeps <- function(.data, scenario, label = 'Label neccessary for scenario A'){
  # Extract only data belonging to representative peptide
  # Input: .data, data frame containing all peptides of one protein
  # scenario: character, can have value of 'A' or 'B'
  # label: character, can have value of 'TRUE' or 'FALSE' required only if scenario == A
  # Output: data frame containing only representative peptide

  # Incomplete labelling - there are cases when in peptides containing multiple Lys/Arg
  # not all of them are heavy in labelled samples. As in SILAC we assume that addition of
  # label does not affect peptide properties, we are taking a mean `Precursor Area` value
  # as the representative `Precursor Area` in such cases
  #.data<- peptide_index[[chooseProtein]]
  if (scenario == 'B'){
    .data <- .data[.data$repPepB,]

  } else if(scenario == 'A' & label == T){
    .data <- .data[.data$repPepA & .data$isLabel,]

  } else if(scenario == 'A' & label == F){
    .data <- .data[.data$repPepA & !.data$isLabel,]
  } else {
    stop(paste('Invalid function use.',label, sep = ' '))
  }

  # summarize to end up with only neccessary columns
  # + handle incomplete labelling
  .data %>%
    group_by(isLabel,Fraction) %>%
    summarise(`Precursor Area` = mean(`Precursor Area`),
              `Protein Group Accessions` = `Protein Group Accessions`[1],
              `Protein Descriptions` = `Protein Descriptions`[1]) %>%
    ungroup()-> .data

}
