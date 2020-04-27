#' Extract Only Data Belonging to Representative Peptide
#'
#' Incomplete labelling - there are cases when in peptides containing multiple
#' Lys/Arg not all of them are heavy in labelled samples. As in SILAC we assume 
#' that addition of label does not affect peptide properties, we are taking 
#' a mean `Precursor Area` value as the representative `Precursor Area` 
#' in such cases.
#'
#' @param .data dataframe containing all peptides of one protein
#' @param scenario character "A", or "B"
#' @param label character, selects for which label state the representative 
#' peptides will be exported, can have value of "TRUE" or "FALSE",
#' required only for scenario "A" 
#'
#' @return dataframe containing only representative peptide
extractRepPeps <- function(.data, scenario, 
                            label = 'Label neccessary for scenario A'){
    if (scenario == 'B'){
        .data <- .data[.data$repPepB,]
        
    } else if(scenario == 'A' & label == TRUE){
        .data <- .data[.data$repPepA & .data$isLabel,]
        
    } else if(scenario == 'A' & label == FALSE){
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
